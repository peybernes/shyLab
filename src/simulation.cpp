/// Implementation of simulation.hpp

#include "simulation.hpp"

#include <algorithm>
#include <cstdlib>
#include <malloc.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI
#include <numeric>
#include <string>

#include <boost/foreach.hpp>
#include <boost/timer.hpp>
#include <sys/time.h>
#include <omp.h>

#include "numerics/kernels.hpp"
#include "numerics/kernel_tools.h"
#include "io/array_io.hpp"
#include "io/init.hpp"
#include "numerical_core.h"
#include "io/output.hpp"
#include "serialize.hpp"
#include "grid/structured_grid.hpp"
#include "variables/variable_metadata.hpp"
#include "euler_riemann_analytic_solver.hpp"

#include <ctime>

static void PrintTimings(const std::vector<RealType>& timings, const std::string& kernel_name) {

  const double mean = accumulate(timings.begin(), timings.end(), 0.0);

  std::cerr << std::setprecision(10) << kernel_name << " : "
	    << "min=" 
	    << *std::min_element(timings.begin(), timings.end()) << "ms, "
	    << "mean=" 
	    << mean / timings.size() << "ms, "
	    << "max=" 
	    << *std::max_element(timings.begin(), timings.end()) << "ms\n";

}

static double get_time() {

  struct timeval tv;
	
  gettimeofday(&tv, NULL);
	
  return (tv.tv_usec / 1000000.0) + tv.tv_sec;

}

#define GFLOP(m, n, s)      ((size_t) (s) * ((m) - 2) * ((n) - 2) * 9 * 2 / 1e9f)
#define MIN(a, b)           (((a) < (b))? (a): (b))

/// Default simulation constructor.
Simulation::Simulation():
  clock(Clock()) {}

void Simulation::CreateVariableDatabase() {

  InitVariableDatabase(&this->variables_database);

}

/// Load simulation from an XML file.
void Simulation::Load(const std::string &filename) {

  // Create an empty property tree object
  ptree empty_ptree;

  ptree pt_file;

  try {
    read_json(filename, pt_file);
  }

  catch (std::exception &e) {
    std::cout << "Error while reading json simulation file:\n\n"
	      << e.what() << "\n\n";
  }

  // Load simulation file.
  ptree pt_simulation;
  pt_simulation = pt_file.get_child("Simulation");

  // Load timetable.
  ptree pt_timetable;
  pt_timetable = pt_simulation.get_child("Timetable", empty_ptree);
  timetable.Load(pt_timetable);

  // Load clock.
  ptree pt_clock;
  pt_clock = pt_simulation.get_child("Clock", empty_ptree);
  clock.Load(pt_clock);

  // Load physical parameters.
  ptree pt_physical_params;
  pt_physical_params = pt_simulation.get_child("PhysicalParams", empty_ptree);
  physical_params.load(pt_physical_params);

  // Load numerical parameters.
  ptree pt_numerical_params;
  pt_numerical_params = pt_simulation.get_child("NumericalParams", empty_ptree);
  numerical_params.load(pt_numerical_params);

  ptree pt_grid;
  pt_grid = pt_simulation.get_child("Grid", empty_ptree);
  m_grid.Load(pt_grid);
  
  // Register variables statically defined in header.
  CreateVariableDatabase();

  ptree pt_variables;
  pt_variables = pt_simulation.get_child("Variables", empty_ptree);
  
  BOOST_FOREACH(ptree::value_type &v, pt_variables) {

    VariableEntry current_entry;
    
    std::string variable_name = v.first;
    ptree pt_variable = v.second;

    std::cerr << "new variable: \"" << variable_name << "\"\n";

    const std::string support_name = pt_variable.get<std::string>("support");

    VariableSupport support;

    if (support_name == "cell") {
     
      support = CELL;
 
    } else if (support_name == "facet") {

      support = FACET;

    } else if (support_name == "vertice") { 

      support = VERTICE;

    } else {

      std::cerr << "ERROR: " << "wrong support for variable \"" << variable_name 
		<< "\", expected \"cell\", \"facet\" or \"vertice\", got \"" 
		<< support_name << "\"\n";

      assert(0);
      
    }

    // Id de la variable dans le variable_store (trie par support).
    int current_variable_id = -1;
    
    switch (support) {

    case CELL:
      current_variable_id =  NumberOfCellVariables(variables_database);
      break;
      
    case VERTICE:
      current_variable_id =  NumberOfVerticeVariables(variables_database);
      break;

    case FACET:
      current_variable_id =  NumberOfFacetVariables(variables_database);
      break;
      

    default:
      assert(0);

    }

    std::cerr << "Current variable id=" << current_variable_id << "\n";

    const int default_attribute = PROTECTED | WRITTEN | COMMUNICATED | INITIALIZABLE;
    VariableRegion default_region = TOTAL;

    VariableEntry variable_entry(current_variable_id, support, default_attribute, default_region);

    AddVariableInDatabase(variable_name, variable_entry, &variables_database);

  }

  std::cerr << variables_database.size() << " variables created in database\n";

  math_parser.RegisterSimulationVariables(variables_database);

  math_parser.SetParameters(physical_params, numerical_params);

  // Load all mathematical expressions.
  ptree pt_math_expressions = 
    pt_simulation.get_child("Parser", empty_ptree);
  
  BOOST_FOREACH(ptree::value_type &v, pt_math_expressions) {

    if (v.first == "Variable")
      math_expressions.push_back(v.second.data());

    else if (v.first == "SecondaryValue") {

      ptree pt_sec_val = pt_math_expressions.get_child("SecondaryValue");
      std::string name = pt_sec_val.get<std::string>("name");
      std::string expr = pt_sec_val.get<std::string>("expression");

      secondary_value_expr[name] = expr;

    }
}
  
  // Load all events.
  ptree pt_events = pt_simulation.get_child("Events", empty_ptree);
  BOOST_FOREACH(ptree::value_type &v, pt_events) {

    // Will hold the (subclassed) event.
    boost::shared_ptr<Event> event_ptr;

    std::string event_name = v.first;
    ptree pt_event = v.second;

    // Any event has a timetable.
    ptree pt_timetable;
    pt_timetable = pt_event.get_child("Timetable", empty_ptree);
    Timetable event_timetable;
    event_timetable.Load(pt_timetable);

    // Shared pointer which will hold the current event.

    if (event_name == "InitVariable") {
      
      InitVariable* iv = new InitVariable(this, event_timetable);
      iv->Load(pt_event);

      event_ptr = boost::shared_ptr<Event>(iv);
      events.push_back(event_ptr);

    } else if (event_name == "OutputSimulation") {

      std::string stream_name = pt_event.get<std::string>("filename");
      const std::string io_format_string =
      pt_event.get<std::string>("format", "ascii");

      OutputSimulation *os = new OutputSimulation(this, stream_name,
						  io_format_string, event_timetable);

      event_ptr = boost::shared_ptr<Event>(os);
      events.push_back(event_ptr);

    } else if (event_name == "OutputMinMax") {

      std::string stream_name = pt_event.get<std::string>("filename");
      OutputMinMax *os = new OutputMinMax(this, stream_name, event_timetable);

      event_ptr = boost::shared_ptr<Event>(os);
      events.push_back(event_ptr);

    } else if (event_name == "OutputLocation") {

      std::string stream_name = pt_event.get<std::string>("filename");

      RealType x = pt_event.get<RealType>("x");
      RealType y = pt_event.get<RealType>("y");

      OutputLocation *os = 
	new OutputLocation(this, stream_name, event_timetable, x, y);

      event_ptr = boost::shared_ptr<Event>(os);
      events.push_back(event_ptr);

    } else if (event_name == "OutputNorms") {

      std::string stream_name = pt_event.get<std::string>("filename");

      const std::string name1 = pt_event.get<std::string>("variable1");
      const std::string name2 = pt_event.get<std::string>("variable2");

      OutputNorms *os = 
	new OutputNorms(this, stream_name, event_timetable, name1, name2);

      event_ptr = boost::shared_ptr<Event>(os);
      events.push_back(event_ptr);
    } else if (event_name == "OutputCartesianLine") {

      std::string stream_name = pt_event.get<std::string>("filename");

      const std::string direction = pt_event.get<std::string>("direction");
      const int indice = pt_event.get<int>("indice");

      OutputCartesianLine *os = 
	new OutputCartesianLine(this, stream_name, event_timetable, direction, indice);

      event_ptr = boost::shared_ptr<Event>(os);
      events.push_back(event_ptr);
    }
  }


  //assert(0);

}

/// Save simulation in XML format.
void Simulation::Save(const std::string &filename) {

  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree pt_simulation;

  // Save simulation clock.
  ptree pt_clock;
  clock.Save(pt_clock);
  pt_simulation.put_child("Clock", pt_clock);

  // Save simulation timetable.
  ptree pt_timetable;
  timetable.Save(pt_timetable);
  pt_simulation.put_child("Timetable", pt_timetable);

  // Save physical parameters.
  ptree pt_physical_params;
  physical_params.save(pt_physical_params);
  pt_simulation.put_child("PhysicalParams", pt_physical_params);

  // Save numerical parameters.
  ptree pt_numerical_params;
  numerical_params.save(pt_numerical_params);
  pt_simulation.put_child("NumericalParams", pt_numerical_params);

  // Save mesh.
  ptree pt_mesh;
  mesh.Save(pt_mesh);
  pt_simulation.put_child("Mesh", pt_mesh);

  // Save all events.
  ptree pt_events;
  BOOST_FOREACH(boost::shared_ptr<Event> &event, events) {
    ptree pt_event;
    event->Save(pt_event);
    pt_events.push_back(ptree::value_type(event->name(), pt_event));
  }
  pt_simulation.put_child("Events", pt_events);

  // Save simulation.
  ptree pt_dump;
  pt_dump.put_child("Simulation", pt_simulation);

  write_json(filename, pt_dump);
}

void Simulation::Init() {

  const RealType zero = 0.0;

  int nb_processes = 1;
  int process_rank = 0;
  
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif // HAVE_MPI

  // Timer for measuring for execution of the code.
  boost::timer simulation_timer = boost::timer();
  simulation_timer.restart();

  // Read mesh from file.
  if (!mesh.filename().empty()) {
    
    std::ifstream ifs(mesh.filename().c_str());

    mesh.ReadGmsh(nb_processes, process_rank, &ifs);

    ifs.close();

  }

//   for (int i = 0; i < nb_processes; ++i) {

//     mesh.mpi_boundary_conditions().at(i).SwitchToLocalNumbering(mesh);
//     mesh.mpi_boundary_conditions().at(i).AllocDataBuffers();

//   }

// #ifdef DEBUG_MPI
//   for (int i = 0; i < nb_processes; ++i) {
    
//     std::cerr << process_rank << "-->" << i << ":\n";
//     mesh.mpi_boundary_conditions().at(i).Print(&std::cerr);
//     std::cerr << "\n";

//   }
// #endif // DEBUG_MPI

  std::cerr << "reading from mesh file: " 
	    << simulation_timer.elapsed() << "\n";

  simulation_timer.restart();
  
  // mesh.ComputeFaceCellConnectivity();
  // mesh.ComputeFacesMap();
  
  simulation_timer.restart();

  const int nb_cells = m_grid.nb_cells();
  assert(0 <= nb_cells);

  const int nb_faces = m_grid.nb_faces();
  assert(0 <= nb_faces);

  const int nb_boundary_faces = 0;
  assert(0 <= nb_boundary_faces);
  
  const int nb_vertices = m_grid.nb_vertices();
  assert(0 <= nb_vertices);

  const int padding = 0;

  const int nb_cell_variables = NumberOfCellVariables(variables_database);
  cell_variables = VariableStore(nb_cell_variables, nb_cells, m_grid.nx(), m_grid.ny(), padding);

  const int nb_facet_variables = NumberOfFacetVariables(variables_database);
  face_variables = VariableStore(nb_facet_variables, nb_faces, m_grid.nx() + 1, m_grid.ny() + 1, padding);

  const int nb_vertice_variables = NumberOfVerticeVariables(variables_database);
  vertice_variables = VariableStore(nb_vertice_variables, nb_vertices, m_grid.nx() + 1, m_grid.ny() + 1, padding);

  std::cerr << "Allocation of simulation variables: " 
	    << simulation_timer.elapsed() << "\n";
  
  simulation_timer.restart();

  m_grid.ComputeGeometricQuantities(&cell_variables, 
				    &face_variables, 
				    &vertice_variables);
  
  std::cerr << "mesh geometrical quantities: "
	    << simulation_timer.elapsed() << "\n";

  simulation_timer.restart();

  // mesh.ReorderFaces(&face_variables);
  
  // mesh.ComputeCellFaceConnectivity();

  // std::ofstream ofs2("after.ppm");
  // mesh.WriteMeshBandwith(ofs2);
  // ofs2.close();

  // std::cerr << "Mesh face reordering and Cell->face connectivity: " 
  // 	    << simulation_timer.elapsed() << "\n";

  simulation_timer.restart();

  BOOST_FOREACH(std::string &expr, math_expressions) {
    
    math_parser.RegisterExpression(expr);
    math_parser.Eval();

  }

  std::cerr << "Registering all parser expressions: " 
	    << simulation_timer.elapsed() << "\n";

  simulation_timer.restart();

  // Execute all events.
  BOOST_FOREACH(boost::shared_ptr<Event> &event, events) {

    if (event->Timetable2().Happens(this->clock) && 
	(event->Type() == PRE)) {
      
      event->Execute();
      event->UpdateTimetableClock(0.0);

    }

  // Remove event if it is finished.
    if (event->Timetable2().IsFinished(this->clock) && 
	(event->Type() == PRE)) {
      
      events.erase(events.begin(), events.begin() + 1);

    }

  }

  BOOST_FOREACH(boost::shared_ptr<Event> &event, events) {
   
    if (event->Timetable2().Happens(this->clock) && 
	event->Type() == POST) {

      event->Execute();

    }
    
    // Remove event if it is finished.
    if (event->Timetable2().IsFinished(this->clock)) {

      events.erase(events.begin(), events.begin() + 1);

    }
  }

  std::cerr << "Evaluating all init events (initializing variables etc.): " 
	    << simulation_timer.elapsed() << "\n";

  // Goes from iteration 0 to iteration 1 ; simulation time remains 0;

  this->Save("dump.json");

}

void Simulation::Run() {
  
  int nb_processes = 1;
  int process_rank = 0;
  
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif // HAVE_MPI

  const int nb_cells = m_grid.nb_cells();
  const int nb_faces = m_grid.nb_faces();
  const int nb_boundary_faces = 0;
  const int nb_interior_faces = 0;

  const int nx = m_grid.nx();
  const int ny = m_grid.ny();

    
  RealType* cell_volumes = cell_variables.GetVariable(variables_database, "cell_volumes");

  // may need to move some variables in local part and be consistant
  // with json files -- here json variable name is ***_volume and
  // variable store in memory ishould be ***_volume_x , same for
  // lagrangian density
  RealType* directional_lagrangian_volume = cell_variables.GetVariable(variables_database, "directional_lagrangian_volume");
  RealType* directional_lagrangian_density = cell_variables.GetVariable(variables_database, "directional_lagrangian_density");
  RealType* directional_lagrangian_volume_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_volume_y");
  RealType* directional_lagrangian_density_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_y");
  
  RealType* in_rho = cell_variables.GetVariable(variables_database, "in_rho");
  RealType* out_rho = cell_variables.GetVariable(variables_database, "out_rho");

  RealType* in_cell_mass = cell_variables.GetVariable(variables_database, "in_cell_mass");
  RealType* out_cell_mass = cell_variables.GetVariable(variables_database, "out_cell_mass");

  RealType* in_e = cell_variables.GetVariable(variables_database, "in_e");
  RealType* out_e = cell_variables.GetVariable(variables_database, "out_e");

  RealType* in_p = cell_variables.GetVariable(variables_database, "in_pressure");

  RealType* rho_ref = cell_variables.GetVariable(variables_database, "rho_ref");
  RealType* p_ref = cell_variables.GetVariable(variables_database, "p_ref");

  RealType* u_lag = vertice_variables.GetVariable(variables_database, "lagrangian_u");
  RealType* v_lag = vertice_variables.GetVariable(variables_database, "lagrangian_v");

  RealType* in_u = vertice_variables.GetVariable(variables_database, "in_u");
  RealType* in_v = vertice_variables.GetVariable(variables_database, "in_v");

  RealType* predicted_u = vertice_variables.GetVariable(variables_database, "predicted_u");
  RealType* predicted_v = vertice_variables.GetVariable(variables_database, "predicted_v");

  RealType* out_u = vertice_variables.GetVariable(variables_database, "out_u");
  RealType* out_v = vertice_variables.GetVariable(variables_database, "out_v");

  RealType* u_ref = vertice_variables.GetVariable(variables_database, "u_ref");

  boost::timer simulation_timer = boost::timer();

  const RealType dx = (m_grid.xmax() - m_grid.xmin()) / nx;
  const RealType dy = (m_grid.ymax() - m_grid.ymin()) / ny;
  
  const RealType length_x = m_grid.xmax() - m_grid.xmin();
  const RealType length_y = m_grid.ymax() - m_grid.ymin();

  const int nb_nodes = (nx + 1) * (ny + 1);

  const int nb_faces_x = (nx + 1) * ny;

  const int nb_faces_y = nx * (ny + 1);

  // Face local variables.
  RealType* volume_fluxes_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* volume_fluxes_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  
  RealType* mass_flux_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* mass_flux_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* energy_flux_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* energy_flux_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  //  const int ALIGN_BYTES = 64;

  // Cell local variables.
  RealType* predicted_pressure = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* cell_pseudo_pressure = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* e_lag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));

  // Node local variables.
  RealType* gradient_v = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_u = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  // INIT
  RealType CFL = 0.3;

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = (nx * iy) + ix;

      const RealType x = ix * dx;
      const RealType y = iy * dy;
      cell_volumes[cell_ooo] = dx * dy;
      in_cell_mass[cell_ooo] = in_rho[cell_ooo] * cell_volumes[cell_ooo];
      in_p[cell_ooo] = 0.4 * in_rho[cell_ooo] * in_e[cell_ooo];
      cell_pseudo_pressure[cell_ooo] = 0.0;
      predicted_pressure[cell_ooo]=0.0;
      e_lag[cell_ooo]=0.0;
      
    }
  }

  for (int iy = 0; iy < ny + 1; ++iy) {
    for (int ix = 0; ix < nx + 1; ++ix) {

      const int node_ooo = ((nx + 1) * iy) +ix;

      u_lag[node_ooo] = 0.0;
      v_lag[node_ooo] = 0.0;
      gradient_u[node_ooo] = 0.0;
      gradient_v[node_ooo] = 0.0;
    }
  }

  std::vector<RealType> time_compressible_euler_physical_to_conservative_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_x_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_x_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_y_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_y_0;

  std::vector<RealType> time_time_step;
  
  std::vector<RealType> time_lagrange_pressure_predicted;
  std::vector<RealType> time_lagrange_velocity_predicted;
  std::vector<RealType> time_lagrange_correction;
  std::vector<RealType> time_lagrange_velocity_correction;
  std::vector<RealType> time_compute_volume_fluxes_X;
  std::vector<RealType> time_gradient_X;
  std::vector<RealType> time_mass_reconstruct_o2_X;
  std::vector<RealType> time_project_mass_X;
  std::vector<RealType> time_reconstruct_energy_o2_X;
  std::vector<RealType> time_project_energy_X;
  std::vector<RealType> time_gradient_nodal_X;
  std::vector<RealType> time_project_nodal_velocity_X;
  std::vector<RealType> time_compute_volume_fluxes_Y;
  std::vector<RealType> time_gradient_Y;
  std::vector<RealType> time_mass_reconstruct_o2_Y;
  std::vector<RealType> time_project_mass_Y;
  std::vector<RealType> time_reconstruct_energy_o2_Y;
  std::vector<RealType> time_project_energy_Y;
  std::vector<RealType> time_gradient_nodal_Y;
  std::vector<RealType> time_project_nodal_velocity_Y;

  std::vector<RealType> time_periodic_boundary;

  // Will hold the range of all simulation variables (including
  // temporaries for each time step).
  std::stringstream sstream;

  sstream << "minmax-pe" 
	  << std::setfill('0') << std::setw(4) 
	  << process_rank << ".txt";

  std::string filename = sstream.str();

  std::ofstream ofs(filename.c_str());

  double time = get_time();

  std::cerr << "Type of Boundary Conditions : " << numerical_params.BoundaryConditions << "\n";      
  std::cerr << "Type of Projection : " << numerical_params.TypeOfProjection << "\n";      

  while (!timetable.IsFinished(clock)) {
    
    if (process_rank == 0)
      std::cerr << "\n"
		<< "Time " << clock.time() << ", iteration " << clock.iter()
		<< "\n";
    

    // // Pack halo exchange data.
    // for (int i = 0; i < nb_processes; ++i)
    //   mesh.mpi_boundary_conditions().at(i).FillSendDataBuffer(cell_variables);

#ifdef HAVE_MPI
    // std::vector<MPI_Request> send_handles(nb_processes);
    // std::vector<MPI_Request> recv_handles(nb_processes);

    // // Start MPI communications.
    // CommunicateAllMpiBC(MPI_COMM_WORLD,
    // 			&(mesh.mpi_boundary_conditions()),
    // 			&send_handles, &recv_handles);

    // // Wait for MPI communications to finish.
    // FinishAllMpiBC(nb_processes, send_handles, recv_handles);

#endif // HAVE_MPI

    // PUT ALL NUMERICS HERE !!!

    const int nx = m_grid.nx();
    const int ny = m_grid.ny();
    const int nz = m_grid.nz();

    const int halo_width = 1;

    assert(halo_width < nx);
    assert(halo_width < ny);

    const RealType dx = (m_grid.xmax() - m_grid.xmin()) / nx;
    const RealType dy = (m_grid.ymax() - m_grid.ymin()) / ny;
    const RealType dz = (m_grid.zmax() - m_grid.zmin()) / nz;

    // Arbitrary. Should be computed using CFL constraint,
    RealType dt;

    const RealType kappa = 1.0 / 3.0;

    struct timespec time1;
    struct timespec time2;

    if (m_grid.dimension() == 2) {
      
      

      clock_gettime(CLOCK_REALTIME, &time1);
      dt = TimeStep(nx, ny, dx, dy, CFL, in_rho, in_p, in_u, in_v);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_time_step.push_back(diff(time1, time2)); 

      std::cerr << "CFL : " << numerical_params.CFL << "\n"; 
      std::cerr << "          dt = " << dt << "\n";
   
      // Lagrange 2D algorithm for compressible Euler equations.
      
      Lagrange2dDriver(//in
		       numerical_params.BoundaryConditions,
		       nx,
		       ny,
		       dx,
		       dy,
		       dt,
		       in_u,
		       in_v,
		       in_e,
		       in_cell_mass,
		       //out
		       in_p,
		       predicted_pressure,
		       cell_pseudo_pressure,
		       predicted_u,
		       predicted_v,
		       e_lag,
		       u_lag,
		       v_lag,
		       // timing
		       time_lagrange_pressure_predicted,
		       time_lagrange_velocity_predicted,
		       time_periodic_boundary,
		       time_lagrange_correction,
		       time_lagrange_velocity_correction);
  
      // Projection 2D algorithm for compressible Euler equations.
  
      if (numerical_params.TypeOfProjection == "AdProjection") {
	if (clock.iter()%2) {    // alternating X + Y and Y + X projection
	// Projection X     X then Y projection
	  AdProjection2dXDriver(//in
				numerical_params.BoundaryConditions,
				nx,
				ny,
				nb_faces_x,
				nb_faces_y,
				nb_cells,
				nb_nodes,
				dx,
				dy,
				dt,
				halo_width,
				predicted_u,
				e_lag,
				u_lag,
				v_lag,
				in_cell_mass, 
				// out
				out_u,
				out_v,
				out_e,
				out_cell_mass,
				directional_lagrangian_volume,
				directional_lagrangian_density,
				volume_fluxes_x,
				mass_flux_x,
				energy_flux_x,
				density_gradient,
				energy_gradient,
				gradient_u,
				gradient_v,
				//timing
				time_compute_volume_fluxes_X,
				time_gradient_X,
				time_mass_reconstruct_o2_X,
				//			      time_project_mass_X,
				time_reconstruct_energy_o2_X,
				time_project_energy_X,
				time_gradient_nodal_X,
				time_project_nodal_velocity_X,
				time_periodic_boundary);
	    
	  std::swap(in_cell_mass, out_cell_mass);
	  std::swap(u_lag, out_u);
	  std::swap(v_lag, out_v);
	  std::swap(e_lag, out_e);
	
	  // Projection Y
	  AdProjection2dYDriver(//in
				numerical_params.BoundaryConditions,
				nx,
				ny,
				nb_faces_x,
				nb_faces_y,
				nb_cells,
				nb_nodes,
				dx,
				dy,
				dt,
				halo_width,
				predicted_v,
				e_lag,
				u_lag,
				v_lag,
				in_cell_mass, 
				//out
				out_u,
				out_v,
				out_e,
				out_cell_mass,
				directional_lagrangian_volume_y,
				directional_lagrangian_density_y,
				volume_fluxes_y,
				mass_flux_y,
				energy_flux_y,
				density_gradient,
				energy_gradient,
				gradient_u,
				gradient_v,
				//timing
				time_compute_volume_fluxes_Y,
				time_gradient_Y,
				time_mass_reconstruct_o2_Y,
				//			      time_project_mass_Y,
				time_reconstruct_energy_o2_Y,
				time_project_energy_Y,
				time_gradient_nodal_Y,
				time_project_nodal_velocity_Y,
				time_periodic_boundary); 
	    
	  std::swap(in_cell_mass, out_cell_mass); 
	  std::swap(in_u, out_u); 
	  std::swap(in_v, out_v);
	  std::swap(in_e, out_e);
	  
	} else { // !clock.iter()%2  Y then X projection
 // Projection Y
	  AdProjection2dYDriver(//in
				numerical_params.BoundaryConditions,
				nx,
				ny,
				nb_faces_x,
				nb_faces_y,
				nb_cells,
				nb_nodes,
				dx,
				dy,
				dt,
				halo_width,
				predicted_v,
				e_lag,
				u_lag,
				v_lag,
				in_cell_mass, 
				//out
				out_u,
				out_v,
				out_e,
				out_cell_mass,
				directional_lagrangian_volume_y,
				directional_lagrangian_density_y,
				volume_fluxes_y,
				mass_flux_y,
				energy_flux_y,
				density_gradient,
				energy_gradient,
				gradient_u,
				gradient_v,
				//timing
				time_compute_volume_fluxes_Y,
				time_gradient_Y,
				time_mass_reconstruct_o2_Y,
				//			      time_project_mass_Y,
				time_reconstruct_energy_o2_Y,
				time_project_energy_Y,
				time_gradient_nodal_Y,
				time_project_nodal_velocity_Y,
				time_periodic_boundary); 

	  std::swap(in_cell_mass, out_cell_mass);
	  std::swap(u_lag, out_u);
	  std::swap(v_lag, out_v);
	  std::swap(e_lag, out_e);
	  
// Projection X
	  AdProjection2dXDriver(//in
				numerical_params.BoundaryConditions,
				nx,
				ny,
				nb_faces_x,
				nb_faces_y,
				nb_cells,
				nb_nodes,
				dx,
				dy,
				dt,
				halo_width,
				predicted_u,
				e_lag,
				u_lag,
				v_lag,
				in_cell_mass, 
				// out
				out_u,
				out_v,
				out_e,
				out_cell_mass,
				directional_lagrangian_volume,
				directional_lagrangian_density,
				volume_fluxes_x,
				mass_flux_x,
				energy_flux_x,
				density_gradient,
				energy_gradient,
				gradient_u,
				gradient_v,
				//timing
				time_compute_volume_fluxes_X,
				time_gradient_X,
				time_mass_reconstruct_o2_X,
				//			      time_project_mass_X,
				time_reconstruct_energy_o2_X,
				time_project_energy_X,
				time_gradient_nodal_X,
				time_project_nodal_velocity_X,
				time_periodic_boundary);

	  std::swap(in_cell_mass, out_cell_mass); 
	  std::swap(in_u, out_u); 
	  std::swap(in_v, out_v);
	  std::swap(in_e, out_e);
	}
	
      } else if (numerical_params.TypeOfProjection == "DirectProjection") {

	// Direct projection --- first order for now.
	DirectProjection2dDriver(//in
				 numerical_params.BoundaryConditions,
				 nx,
				 ny,
				 dx,
				 dy,
				 dt,
				 halo_width,
				 predicted_u,
				 predicted_v,
				 //out
				 e_lag,
				 u_lag,
				 v_lag,
				 in_u,
				 in_v,
				 in_e,
				 in_cell_mass,
				 out_u,
				 out_v,
				 out_e,
				 out_cell_mass,
				 directional_lagrangian_volume,
				 directional_lagrangian_density,
				 directional_lagrangian_volume_y,
				 directional_lagrangian_density_y,
				 volume_fluxes_x,
				 volume_fluxes_y,
				 mass_flux_x,
				 mass_flux_y,
				 energy_flux_x,
				 energy_flux_y,
				 //timing
				 time_compute_volume_fluxes_X,
				 time_gradient_X,
				 time_mass_reconstruct_o2_X,
				 time_project_mass_X,
				 time_reconstruct_energy_o2_X,
				 time_project_energy_X,
				 time_gradient_nodal_X,
				 time_project_nodal_velocity_X,
				 time_compute_volume_fluxes_Y,
				 time_gradient_Y,
				 time_mass_reconstruct_o2_Y,
				 time_project_mass_Y,
				 time_reconstruct_energy_o2_Y,
				 time_project_energy_Y,
				 time_gradient_nodal_Y,
				 time_project_nodal_velocity_Y,
				 time_periodic_boundary);
	
	std::swap(in_cell_mass, out_cell_mass); 
	std::swap(in_u, out_u); 
	std::swap(in_v, out_v);
	std::swap(in_e, out_e);
	
      }


      // Recompute the density from mass and volume. Only needed for output.

      for (int iy = 0; iy < ny; ++iy) {
	for (int ix = 0; ix < nx; ++ix) {
	  
	  const int cell_ooo = (nx * iy) + ix;
	  
	  const RealType x = ix * dx;
	  const RealType y = iy * dy;

	  in_rho[cell_ooo] =  in_cell_mass[cell_ooo] / cell_volumes[cell_ooo];
	}
      }

      // Compute analytical solution for the 1D Riemann problem.
      const RealType rho_left = 1.0;
      const RealType u_left = 0.0;
      const RealType p_left = 1.0;

      const RealType rho_right = 0.125;
      const RealType u_right = 0.0;
      const RealType p_right = 0.1;

      const RealType gamma = 1.4;

      const RealType x0 = 0.5;

      const RealType t = clock.time();
      
      RealType* cell_centers_x = cell_variables(CELL_CENTERS_X);

      RiemannAnalyticalSolver(rho_left, u_left, p_left,
      			      rho_right, u_right, p_right,
      			      gamma,
      			      t, x0, 
      			      m_grid.nb_cells(),
      			      cell_centers_x,
      			      rho_ref, u_ref, p_ref);
			      

      //assert(0);

    } else if (m_grid.dimension() == 3) {

      //AdvectionFvUw13d(nx, ny, nz, dt, dx, dy, dz, in_rho, u_face, v_face, w_face, out_rho);

    }

    // // Unpack halo exchange data.
    // for (int i = 0; i < nb_processes; ++i)
    //   mesh.mpi_boundary_conditions().at(i).UpdateVariables(&cell_variables);

    if (process_rank == 0)
      std::cerr << "-----------------------------\n"
		<< "timestep = " << dt << "\n";

    clock.Update(dt);

    std::cerr << clock.time() << "\n";

    // Execute all post events that actually happen.
    BOOST_FOREACH(boost::shared_ptr<Event> &event, events)
      event->UpdateTimetableClock(dt);
    
    BOOST_FOREACH(boost::shared_ptr<Event> &event, events) {

      if (event->Timetable2().Happens(this->clock)) {

	event->Execute();
	event->ResetTimetableClock();
    
      }

      // Remove event if finished.
      if (event->Timetable2().IsFinished(this->clock)) {
      
	events.erase(events.begin(), events.begin() + 1);

      }
    }
  }

  
  time = get_time() - time;

  std::cerr << "Elapsed time for simulation: " << time << "s\n";

  std::cerr << "\n-----------------------------------------------------------------------------------------\n\n";

  PrintTimings(time_time_step,                   "TimeStep                               ");
  std::cerr << "\n\n\n";
  PrintTimings(time_lagrange_pressure_predicted, "LagrangePressurePredicted              ");
  PrintTimings(time_lagrange_velocity_predicted, "LagrangeVelocityPredicted              ");
  PrintTimings(time_lagrange_correction,         "LagrangeCorrection                     ");
  PrintTimings(time_lagrange_velocity_correction,"LagrangeVelocityCorrection             ");
  std::cerr << "\n\n\n";
  PrintTimings(time_compute_volume_fluxes_X,     "ComputeDirectionalLagrangianQuantitiesX");
  PrintTimings(time_compute_volume_fluxes_Y,     "ComputeDirectionalLagrangianQuantitiesY");
  std::cerr << "\n";
  PrintTimings(time_gradient_X,                  "ReconstructGradientX                   ");
  PrintTimings(time_gradient_Y,                  "ReconstructGradientY                   ");
  std::cerr << "\n";
  PrintTimings(time_mass_reconstruct_o2_X,       "ReconstructMassFluxOrder2X             ");
  PrintTimings(time_mass_reconstruct_o2_Y,       "ReconstructMassFluxOrder2Y             ");
  std::cerr << "\n";
  //  PrintTimings(time_project_mass_X,              "ProjectMassX                           ");
  //  PrintTimings(time_project_mass_Y,              "ProjectMassY                           ");
  std::cerr << "\n";
  PrintTimings(time_reconstruct_energy_o2_X,     "ReconstructIntensiveVariableFluxOrder2X");
  PrintTimings(time_reconstruct_energy_o2_Y,     "ReconstructIntensiveVariableFluxOrder2Y");
  std::cerr << "\n";
  PrintTimings(time_project_energy_X,            "MassProjectIntensiveVariableX          ");
  PrintTimings(time_project_energy_Y,            "MassProjectIntensiveVariableY          ");
  std::cerr << "\n";
  PrintTimings(time_gradient_nodal_X,            "ReconstructNodalGradientX              ");
  PrintTimings(time_gradient_nodal_Y,            "ReconstructNodalGradientY              ");
  std::cerr << "\n";
  PrintTimings(time_project_nodal_velocity_X,    "ProjectNodalIntensiveVariableOrder2X   ");
  PrintTimings(time_project_nodal_velocity_Y,    "ProjectNodalIntensiveVariableOrder2Y   ");



  std::cerr << "\n";


  std::cerr << "\n";


  std::cerr << "\n";
  PrintTimings(time_periodic_boundary,           "**Aggregated** boundary functions      ");
  std::cerr << "\n-----------------------------------------------------------------------------------------\n\n";

  ofs.close();
}

void Simulation::Finalize() {

  cell_variables.DeAllocate();
  face_variables.DeAllocate();
  boundary_face_variables.DeAllocate();
  vertice_variables.DeAllocate();
  
}

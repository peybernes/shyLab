/// Implementation of simulation.hpp

#include "simulation.hpp"

#include <algorithm>
#include <cstdlib>
#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI
#include <numeric>
#include <string>

#include <boost/foreach.hpp>
#include <boost/timer.hpp>
#include <sys/time.h>
#include <omp.h>

#include "array_io.hpp"
#include "init.hpp"
#include "kernels.hpp"
#include "numerical_core.h"
#include "output.hpp"
#include "serialize.hpp"
#include "structured_grid.hpp"
#include "variable_metadata.hpp"

#include <ctime>

static double diff(timespec start, timespec end) {

  timespec result;

  if ((end.tv_nsec-start.tv_nsec) < 0) {

    result.tv_sec = end.tv_sec-start.tv_sec - 1;
    result.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;

  } else {

    result.tv_sec = end.tv_sec - start.tv_sec;
    result.tv_nsec = end.tv_nsec - start.tv_nsec;

  }

  const double  ms_time = 
    (1000.0 * result.tv_sec) + (result.tv_nsec / 1000000.0);

  return ms_time;
}


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

  CreateVariableDatabase();

  std::cerr << variables_database.size() << " variables\n";

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
  
  math_parser.RegisterSimulationVariables(variables_database);

  math_parser.SetParameters(physical_params, numerical_params);

  // Load the mesh.
  // ptree pt_mesh;
  // pt_mesh = pt_simulation.get_child("Mesh", empty_ptree);
  // mesh.Load(pt_mesh);
  // std::cerr << "Mesh filename: " << mesh.filename() << "\n";

  ptree pt_grid;
  pt_grid = pt_simulation.get_child("Grid", empty_ptree);
  m_grid.Load(pt_grid);

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

  // Write simulation in XML format.
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

  const int padding = 0;
  
  const int nb_cells = m_grid.nb_cells();
  assert(0 <= nb_cells);

  const int nb_faces = m_grid.nb_faces();
  assert(0 <= nb_faces);

  const int nb_boundary_faces = 0;
  assert(0 <= nb_boundary_faces);
  
  const int nb_vertices = m_grid.nb_vertices();
  assert(0 <= nb_vertices);

  //assert(0);

  cell_variables = VariableStore(NB_CELL_VALUES, nb_cells, m_grid.nx(), m_grid.ny(), padding);

  face_variables = VariableStore(NB_FACE_VALUES, nb_faces, m_grid.nx(), m_grid.ny(), padding);

  // boundary_face_variables = VariableStore(NB_BOUNDARY_FACE_VALUES, 
  // 					  nb_boundary_faces, m_grid.nx(), m_grid.ny(), padding);

  vertice_variables = VariableStore(NB_VERTICE_VALUES, nb_vertices, m_grid.nx(), m_grid.ny(), padding);

  std::cerr << "Allocation of simulation variables: " 
	    << simulation_timer.elapsed() << "\n";
  
  simulation_timer.restart();

  m_grid.ComputeGeometricQuantities(&cell_variables, 
				    &face_variables, 
				    &vertice_variables);
  
  // mesh.ComputeGeometricQuantities(&cell_variables, 
  // 				  &face_variables, 
  // 				  &vertice_variables);

  std::cerr << "mesh geometrical quantities: "
	    << simulation_timer.elapsed() << "\n";

  simulation_timer.restart();

  // mesh.ReorderFaces(&face_variables);
  
  // mesh.ComputeCellFaceConnectivity();

  // std::ofstream ofs2("after.ppm");
  // mesh.WriteMeshBandwith(ofs2);
  // ofs2.close();

  std::cerr << "Mesh face reordering and Cell->face connectivity: " 
	    << simulation_timer.elapsed() << "\n";

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
  
  RealType* cell_volumes = cell_variables(CELL_VOLUMES);

  RealType* local_timestep = cell_variables(OUT_TIMESTEP);


  RealType* in_rho = cell_variables(IN_RHO);
  RealType* in_u = cell_variables(IN_U);
  RealType* in_v = cell_variables(IN_V);
  RealType* in_w = cell_variables(IN_W);
  RealType* in_e = cell_variables(IN_E);

  RealType* predicted_rho = cell_variables(PREDICTED_RHO);
  RealType* predicted_u = cell_variables(PREDICTED_U);
  RealType* predicted_v = cell_variables(PREDICTED_V);
  RealType* predicted_w = cell_variables(PREDICTED_W);
  RealType* predicted_e = cell_variables(PREDICTED_E);

  RealType* out_rho = cell_variables(OUT_RHO);
  RealType* out_u = cell_variables(OUT_U);
  RealType* out_v = cell_variables(OUT_V);
  RealType* out_w = cell_variables(OUT_W);
  RealType* out_e = cell_variables(OUT_E);

  boost::timer simulation_timer = boost::timer();

  std::vector<RealType> time_compressible_euler_physical_to_conservative_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_x_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_x_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_y_0;
  std::vector<RealType> time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_y_0;
  std::vector<RealType> time_compressible_euler_conservative_to_physical_0;
  
  // Will hold the range of all simulation variables (including
  // temporaries for each time step).
  std::stringstream sstream;

  sstream << "minmax-pe" 
	  << std::setfill('0') << std::setw(4) 
	  << process_rank << ".txt";

  std::string filename = sstream.str();

  std::ofstream ofs(filename.c_str());

  double time = get_time();

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

    const int halo_width = 2;

    assert(halo_width < nx);
    assert(halo_width < ny);

    const RealType dx = (m_grid.xmax() - m_grid.xmin()) / nx;
    const RealType dy = (m_grid.ymax() - m_grid.ymin()) / ny;
    const RealType dz = (m_grid.zmax() - m_grid.zmin()) / nz;

    // Arbitrary. Should be computed using CFL constraint,
    const RealType dt = 0.1 * dx;

    const RealType kappa = 1.0 / 3.0;

    std::cerr << "dt=" << dt << ", dx=" << dx << ", dy=" << dy << "\n";

    struct timespec time1;
    struct timespec time2;

    if (m_grid.dimension() == 2) {

      // PREDICTION.
      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerPhysicalToConservative(nx, ny, in_rho, in_u, in_v, in_e,
					      predicted_rho, predicted_u, predicted_v, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_physical_to_conservative_0.push_back(diff(time1, time2));


      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerFvUwKappa2dX(nx, ny, halo_width, dt, dy, dx, kappa, in_rho, in_u, in_v, in_e, 
				    predicted_rho, predicted_u, predicted_v, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_fv_uw_kappa_2d_x_0.push_back(diff(time1, time2));


      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerFvUwKappa2dBoundaryConditionsX(nx, ny, halo_width, dt, dy, dx, kappa, in_rho, in_u, in_v, in_e, 
      						      predicted_rho, predicted_u, predicted_v, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_x_0.push_back(diff(time1,time2));


      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerFvUwKappa2dY(nx, ny, halo_width, dt, dy, dx, kappa, in_rho, in_v, in_u, in_e, 
				    predicted_rho, predicted_v, predicted_u, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_fv_uw_kappa_2d_y_0.push_back(diff(time1, time2));

      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerFvUwKappa2dBoundaryConditionsY(nx, ny, halo_width, dt, dy, dx, kappa, in_rho, in_v, in_u, in_e, 
      						      predicted_rho, predicted_v, predicted_u, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_y_0.push_back(diff(time1,time2));

      clock_gettime(CLOCK_REALTIME, &time1);
      CompressibleEulerConservativeToPhysical(nx, ny, predicted_rho, predicted_u, predicted_v, predicted_e);
      clock_gettime(CLOCK_REALTIME, &time2);
      time_compressible_euler_conservative_to_physical_0.push_back(diff(time1,time2));

      // CORRECTION.

      CompressibleEulerPhysicalToConservative(nx, ny, predicted_rho, predicted_u, predicted_v, predicted_e,
      					      out_rho, out_u, out_v, out_e);

      CompressibleEulerFvUwKappa2dX(nx, ny, halo_width, dt, dy, dx, kappa, predicted_rho, predicted_u, predicted_v, predicted_e, 
      				    out_rho, out_u, out_v, out_e);

      CompressibleEulerFvUwKappa2dBoundaryConditionsX(nx, ny, halo_width, dt, dy, dx, kappa, predicted_rho, predicted_u, predicted_v, predicted_e, 
      						      out_rho, out_u, out_v, out_e);

      CompressibleEulerFvUwKappa2dY(nx, ny, halo_width, dt, dy, dx, kappa, predicted_rho, predicted_v, predicted_u, predicted_e, 
      				    out_rho, out_v, out_u, out_e);

      CompressibleEulerFvUwKappa2dBoundaryConditionsY(nx, ny, halo_width, dt, dy, dx, kappa, predicted_rho, predicted_v, predicted_u, predicted_e, 
      						      out_rho, out_v, out_u, out_e);

      CompressibleEulerConservativeToPhysical(nx, ny, out_rho, out_u, out_v, out_e);
      
	
      for (int iy = 0; iy < ny; ++iy) {
	for (int ix = 0; ix < nx; ++ix) {

	  const int i = (nx * iy) + ix;

	  out_rho[i] = 0.5 * (in_rho[i] + predicted_rho[i]);
	  out_u[i] = 0.5 * (in_u[i] + predicted_u[i]);
	  out_v[i] = 0.5 * (in_v[i] + predicted_v[i]);
	  out_e[i] = 0.5 * (in_e[i] + predicted_e[i]);

	}
      }

      //assert(0);

    } else if (m_grid.dimension() == 3) {

      //AdvectionFvUw13d(nx, ny, nz, dt, dx, dy, dz, in_rho, u_face, v_face, w_face, out_rho);

    }

    std::swap(in_rho, out_rho);
    std::swap(in_u, out_u);
    std::swap(in_v, out_v);
    std::swap(in_e, out_e);

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

  PrintTimings(time_compressible_euler_physical_to_conservative_0, "CompressibleEulerPhysicalToConservative");
  PrintTimings(time_compressible_euler_conservative_to_physical_0, "CompressibleEulerConservativeToPhysical");

  std::cerr << "\n";

  PrintTimings(time_compressible_euler_fv_uw_kappa_2d_x_0, "CompressibleEulerFvUwKappa2dX");
  PrintTimings(time_compressible_euler_fv_uw_kappa_2d_y_0, "CompressibleEulerFvUwKappa2dY");

  std::cerr << "\n";

  PrintTimings(time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_x_0, "CompressibleEulerFvUwKappaBoundaryConditionsX");
  PrintTimings(time_compressible_euler_fv_uw_kappa_2d_boundary_conditions_y_0, "CompressibleEulerFvUwKappaBoundaryConditionsY");

  std::cerr << "\n-----------------------------------------------------------------------------------------\n\n";

  ofs.close();
}

void Simulation::Finalize() {

  cell_variables.DeAllocate();
  face_variables.DeAllocate();
  boundary_face_variables.DeAllocate();
  vertice_variables.DeAllocate();
  
}

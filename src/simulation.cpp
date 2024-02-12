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
#include <iostream>
#include <fstream>

using namespace std;

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


//#include <likwid.h>

#include <ctime>

#define PI 3.14159265358979323846

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

  int nb_mat = numerical_params.NumberOfMaterials;

  RealType speed_of_sound_mix[nb_cells];
  RealType gamma_mix[nb_cells];
  RealType pi_prime_mix[nb_cells];
  RealType gamma_k[nb_mat];
  RealType pi_prime_k   [nb_mat];
  gamma_k[0] = physical_params.gamma_1;
  if (nb_mat > 1)  gamma_k[1] = physical_params.gamma_2;
  if (nb_mat > 2)  gamma_k[2] = physical_params.gamma_3;
  if (nb_mat > 3)  gamma_k[3] = physical_params.gamma_4;

  pi_prime_k[0] = physical_params.pi_1 / gamma_k[0];
  if (nb_mat > 1)  pi_prime_k[1] = physical_params.pi_2 / gamma_k[1];
  if (nb_mat > 2)  pi_prime_k[2] = physical_params.pi_3 / gamma_k[2];
  if (nb_mat > 3)  pi_prime_k[3] = physical_params.pi_4 / gamma_k[3];
 
  RealType* cell_volumes = cell_variables.GetVariable(variables_database, "cell_volumes");

  // may need to move some variables in local part and be consistant
  // with json files -- here json variable name is ***_volume and
  // variable store in memory ishould be ***_volume_x , same for
  // lagrangian density
  RealType* directional_lagrangian_volume = cell_variables.GetVariable(variables_database, "directional_lagrangian_volume");
  RealType* directional_lagrangian_density = cell_variables.GetVariable(variables_database, "directional_lagrangian_density");
  RealType* directional_lagrangian_density_1 = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_1");
  RealType* directional_lagrangian_density_2 = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_2");
  RealType* directional_lagrangian_density_3 = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_3");
  RealType* directional_lagrangian_volume_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_volume_y");
  RealType* directional_lagrangian_density_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_y");
  RealType* directional_lagrangian_density_1_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_1_y");
  RealType* directional_lagrangian_density_2_y = cell_variables.GetVariable(variables_database, "directional_lagrangian_density_2_y");
  
  RealType* in_cell_mass = cell_variables.GetVariable(variables_database, "in_cell_mass");
  RealType* out_cell_mass = cell_variables.GetVariable(variables_database, "out_cell_mass");
  RealType* in_cell_mass_1 = cell_variables.GetVariable(variables_database, "in_cell_mass_1");
  RealType* out_cell_mass_1 = cell_variables.GetVariable(variables_database, "out_cell_mass_1");
  RealType* in_cell_mass_2 = cell_variables.GetVariable(variables_database, "in_cell_mass_2");
  RealType* out_cell_mass_2 = cell_variables.GetVariable(variables_database, "out_cell_mass_2");
  RealType* in_cell_mass_3 = cell_variables.GetVariable(variables_database, "in_cell_mass_3");
  RealType* out_cell_mass_3 = cell_variables.GetVariable(variables_database, "out_cell_mass_3");


  RealType* in_cell_volumic_fraction = cell_variables.GetVariable(variables_database, "in_cell_volumic_fraction");
  RealType* out_cell_volumic_fraction = cell_variables.GetVariable(variables_database, "out_cell_volumic_fraction");

  RealType* in_c_1   = cell_variables.GetVariable(variables_database, "in_c_1");
  RealType* out_c_1  = cell_variables.GetVariable(variables_database, "out_c_1");
  RealType* in_c_2   = cell_variables.GetVariable(variables_database, "in_c_2");
  RealType* out_c_2  = cell_variables.GetVariable(variables_database, "out_c_2");
  RealType* in_c_3   = cell_variables.GetVariable(variables_database, "in_c_3");
  RealType* out_c_3  = cell_variables.GetVariable(variables_database, "out_c_3");
  RealType** in_c_k  = new RealType*[nb_mat];
  RealType** out_c_k = new RealType*[nb_mat];
  in_c_k [0]  =  in_c_1;
  out_c_k[0]  =  out_c_1;
  if (nb_mat > 1) {
    in_c_k [1]  =  in_c_2;
    out_c_k[1]  =  out_c_2;
  }
  if (nb_mat > 2) {
    in_c_k [2]  = in_c_3;
    out_c_k[2]  = out_c_3;        
  }
  if (nb_mat > 3) {
    in_c_k [3]  = new RealType[nb_cells];
    out_c_k[3]  = new RealType[nb_cells];        
  }
  
  RealType* in_e = cell_variables.GetVariable(variables_database, "in_e");
  RealType* out_e = cell_variables.GetVariable(variables_database, "out_e");
  RealType* in_e_1 = cell_variables.GetVariable(variables_database, "in_e_1");
  RealType* out_e_1 = cell_variables.GetVariable(variables_database, "out_e_1");
  RealType* in_e_2 = cell_variables.GetVariable(variables_database, "in_e_2");
  RealType* out_e_2 = cell_variables.GetVariable(variables_database, "out_e_2");
  RealType* in_e_3 = cell_variables.GetVariable(variables_database, "in_e_3");
  RealType* out_e_3 = cell_variables.GetVariable(variables_database, "out_e_3");
  RealType** in_e_k  = new RealType*[nb_mat];
  RealType** out_e_k = new RealType*[nb_mat];
  in_e_k [0]  =  in_e_1;
  out_e_k[0]  =  out_e_1;
  if (nb_mat > 1) {
    in_e_k [1]  =  in_e_2;
    out_e_k[1]  =  out_e_2;
  }
  if (nb_mat > 2) {
    in_e_k [2]  = in_e_3;
    out_e_k[2]  = out_e_3;        
  }
  if (nb_mat > 3) {
    in_e_k [3]  = new RealType[nb_cells];
    out_e_k[3]  = new RealType[nb_cells];        
  }

  RealType* in_rho = cell_variables.GetVariable(variables_database, "in_rho");
  RealType* out_rho = cell_variables.GetVariable(variables_database, "out_rho");
  RealType* in_rho_1 = cell_variables.GetVariable(variables_database, "in_rho_1");
  RealType* out_rho_1 = cell_variables.GetVariable(variables_database, "out_rho_1");
  RealType* in_rho_2 = cell_variables.GetVariable(variables_database, "in_rho_2");
  RealType* out_rho_2 = cell_variables.GetVariable(variables_database, "out_rho_2");
  RealType* in_rho_3 = cell_variables.GetVariable(variables_database, "in_rho_3");
  RealType* out_rho_3 = cell_variables.GetVariable(variables_database, "out_rho_3");
  RealType** in_rho_k  = new RealType*[nb_mat];
  RealType** out_rho_k = new RealType*[nb_mat];
  in_rho_k [0]  =  in_rho_1;
  out_rho_k[0]  =  out_rho_1;
  if (nb_mat > 1) {
    in_rho_k [1]  =  in_rho_2;
    out_rho_k[1]  =  out_rho_2;
  }
  if (nb_mat > 2) {
    in_rho_k [2]  = in_rho_3;
    out_rho_k[2]  = out_rho_3;        
  }
  if (nb_mat > 3) {
    in_rho_k [3]  = new RealType[nb_cells];
    out_rho_k[3]  = new RealType[nb_cells];        
  }
 
  RealType* in_p = cell_variables.GetVariable(variables_database, "in_pressure");
  RealType* in_p_1 = cell_variables.GetVariable(variables_database, "in_pressure_1");
  RealType* in_p_2 = cell_variables.GetVariable(variables_database, "in_pressure_2");
  RealType* in_p_3 = cell_variables.GetVariable(variables_database, "in_pressure_3");
  RealType** in_p_k  = new RealType*[nb_mat];
  RealType** out_p_k = new RealType*[nb_mat];
  in_p_k [0]  =  in_p_1;
  if (nb_mat > 1) {
    in_p_k [1]  =  in_p_2;
  }
  if (nb_mat > 2) {
    in_p_k [2]  = in_p_3;
  }
  if (nb_mat > 3) {
    in_p_k [3]  = new RealType[nb_cells];
  }

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

  RealType* in_u_cell  = cell_variables.GetVariable(variables_database, "in_u_cell");
  RealType* out_u_cell = cell_variables.GetVariable(variables_database, "out_u_cell");
  RealType* in_v_cell  = cell_variables.GetVariable(variables_database, "in_v_cell");
  RealType* out_v_cell = cell_variables.GetVariable(variables_database, "out_v_cell");


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
  RealType* volume_fluxes_1_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* volume_fluxes_1_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));  
  RealType* volume_fluxes_2_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* volume_fluxes_2_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  
  RealType* mass_flux_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* mass_flux_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* mass_flux_1_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* mass_flux_1_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* mass_flux_2_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* mass_flux_2_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* energy_flux_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* energy_flux_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* energy_flux_1_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* energy_flux_1_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* energy_flux_2_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* energy_flux_2_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* reconstructed_density_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_density_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* reconstructed_density_1_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_density_1_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* reconstructed_density_2_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_density_2_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* reconstructed_energy_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_energy_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* reconstructed_energy_1_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_energy_1_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));
  RealType* reconstructed_energy_2_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_energy_2_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* reconstructed_concentration_faces_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* reconstructed_concentration_faces_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* bool_check_fluxes_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* bool_check_fluxes_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  RealType* concentration_flux_x = (RealType*) memalign(ALIGN_BYTES, nb_faces_x * sizeof(RealType));
  RealType* concentration_flux_y = (RealType*) memalign(ALIGN_BYTES, nb_faces_y * sizeof(RealType));

  //  const int ALIGN_BYTES = 64;

  // Cell local variables.
  RealType* predicted_pressure = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* predicted_pressure_1 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* predicted_pressure_2 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* cell_pseudo_pressure = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* e_lag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* e_1_lag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* e_2_lag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_1_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_2_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_1_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_2_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_1_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_2_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_1_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* density_2_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_1_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_2_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_1_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_2_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_1_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_2_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_1_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* energy_2_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* concentration_gradient = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* concentration_gradient_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* concentration_gradient_diag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* concentration_gradient_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* interface_normal_x = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* interface_normal_y = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_rho_0 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_cell_volumic_fraction_0 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_total_energy = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* out_total_energy = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_y_1 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_y_2 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* in_y_3 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* out_y_1 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* out_y_2 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* out_y_3 = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* rho_total_energy = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* rho_e = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* rho_U = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* rho_V = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* beta = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType));
  RealType* p_plus_pi_prime = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_plus_pi_prime_gradx = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_plus_pi_prime_grady = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_e_gradx_left  = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_e_gradx_right = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_e_grady_top   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_e_grady_bot   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_gradx_left  = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_gradx_right = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_grady_top   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* p_grady_bot   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_gradx_left  = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_gradx_right = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_grady_top   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  RealType* rho_grady_bot   = (RealType*) memalign(ALIGN_BYTES, nb_cells * sizeof(RealType)); 
  
  // arrays of pointers for multimat
  RealType** alphak_gradx_left   = new RealType*[nb_mat];
  RealType** alphak_gradx_right  = new RealType*[nb_mat];
  RealType** alphak_grady_bot    = new RealType*[nb_mat];
  RealType** alphak_grady_top    = new RealType*[nb_mat];
  RealType** alpha_beta_k        = new RealType*[nb_mat];
  RealType** alpha_beta_k_tmp    = new RealType*[nb_mat];
  RealType** masse_k             = new RealType*[nb_mat];
  RealType** masse_k_tmp         = new RealType*[nb_mat];
  for (int k = 0;k < nb_mat; k++) {
    alphak_gradx_left [k]  = new RealType[nb_cells];
    alphak_gradx_right[k]  = new RealType[nb_cells];
    alphak_grady_bot  [k]  = new RealType[nb_cells];
    alphak_grady_top  [k]  = new RealType[nb_cells];
    masse_k           [k]  = new RealType[nb_cells];
    masse_k_tmp       [k]  = new RealType[nb_cells];
    alpha_beta_k      [k]  = new RealType[nb_cells];
    alpha_beta_k_tmp  [k]  = new RealType[nb_cells];
  }  
  
  // Node local variables.
  RealType* gradient_v = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_v_y = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_v_diag = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_v_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_u = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_u_y = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_u_diag = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* gradient_u_antidiag = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* mass_corner_fluxes = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* mass_1_corner_fluxes = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* mass_2_corner_fluxes = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* volume_fluxes_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* volume_fluxes_1_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* volume_fluxes_2_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* energy_flux_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* energy_flux_1_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* energy_flux_2_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* concentration_flux_corner = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  index_t* sign_x_corner_fluxes = (index_t*) memalign(ALIGN_BYTES, nb_nodes * sizeof(index_t));
  index_t* sign_y_corner_fluxes = (index_t*) memalign(ALIGN_BYTES, nb_nodes * sizeof(index_t));
  RealType* in_X_x = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));
  RealType* in_X_y = (RealType*) memalign(ALIGN_BYTES, nb_nodes * sizeof(RealType));

  /*
   // FOR TEST HAAS

  //RealType xC = 57.0;
  //RealType yC = 4.45;
  //RealType r = 2.5;
  //If change, change in file_init --> INIT_MMX/  
  ifstream file_ini_HAAS("../INIT_MMX/fracs_HAAS.plt", ios::in);
  if (file_ini_HAAS) {
    string line;
    index_t ix, iy;
    RealType vol_fraction_1, vol_fraction_2;
    
    while (getline(file_ini_HAAS,line)) {
      file_ini_HAAS >> ix >> iy >> vol_fraction_2 >> vol_fraction_1;
      index_t cell_ooo = (iy - 1) * nx + (ix - 1);
      in_cell_volumic_fraction[cell_ooo] = vol_fraction_1;
    }
    file_ini_HAAS.close();
  } else {   
    cout << "Impossible to open INIT FILE !" << endl;
  }
  */


  /*  
  // FOR TEST IMPACT ED

  //RealType xC = 5.0;
  //RealType yC = 0.0;
  //RealType r = 2.0;
  //If change, change in file_init --> INIT_MMX/  
  ifstream file_ini_IMPACT("../INIT_MMX/fracs_IMPACT.plt", ios::in);
  if (file_ini_IMPACT) {
    string line;
    index_t ix, iy;
    RealType vol_fraction_1, vol_fraction_2, x;
    
    while (getline(file_ini_IMPACT,line)) {
      file_ini_IMPACT >> ix >> iy >> vol_fraction_2 >> vol_fraction_1;
      index_t cell_ooo = (iy - 1) * nx + (ix - 1);
      in_cell_volumic_fraction[cell_ooo] = vol_fraction_1;
      x = (ix - 1) * dx;
      if ((x < 1.5) || (x == 1.5)) {
	in_cell_volumic_fraction[cell_ooo] = 1.0;
      }
    }
    file_ini_IMPACT.close();
  } else {   
    cout << "Impossible to open INIT FILE !" << endl;
  }
  */



  /*
  // FOR TEST IMPACT ED SQUARE
  RealType xC = 5.0;
  RealType yC = 0.0;
  RealType r = 2.0;
  for (index_t ix = 0; ix < nx; ix++) {
    for (index_t iy = 0; iy < ny; iy++) {
      index_t cell_ooo = nx * iy + ix;
      RealType x = ix * dx;
      RealType y = iy * dy;
      RealType d = std::max(fabs(x - xC), fabs(y - yC));
      if ((d < r) || (d == r) || (x < 1.5) || (x == 1.5)) {
	in_cell_volumic_fraction[cell_ooo] = 1.0;
      }
    }
  }
  */

  /*
  // FOR ADVECTION OF A TILTED SQUARE TEST
  RealType xC = 0.4;
  RealType yC = 0.4;
  RealType r = 0.2;
  for (index_t ix = 0; ix < nx; ix++) {
    for (index_t iy = 0; iy < ny; iy++) {
      index_t cell_ooo = nx * iy + ix;
      RealType x = ix * dx;
      RealType y = iy * dy;
      RealType d = fabs(x - xC) + fabs(y - yC);
      if ((d < r) || (d ==r)) {
	in_rho[cell_ooo] = 10.;
	//in_cell_volumic_fraction[cell_ooo] = 1.;
      }
    }
  }
  */

  /*
  // FOR ADVECTION OF A CIRCLE TEST
  RealType xC = 1.4;
  RealType yC = 1.4;
  RealType r = 0.2;
  for (index_t ix = 0; ix < nx; ix++) {
    for (index_t iy = 0; iy < ny; iy++) {
      index_t cell_ooo = nx * iy + ix;
      RealType x = ix * dx;
      RealType y = iy * dy;
      RealType d2 = (x - xC) * (x - xC) + (y - yC) * (y - yC);
      RealType d = std::sqrt(d2);
      if ((d < r) || (d ==r)) {
	//in_rho[cell_ooo] = 10.;
	in_cell_volumic_fraction[cell_ooo] = 1.;
      }
    }
  }
  */


// INIT

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = (nx * iy) + ix;

      const RealType x = ix * dx;
      const RealType y = iy * dy;

      cell_volumes[cell_ooo] = dx * dy;

      if (numerical_params.TypeOfModel == "Monomaterial") {

	in_cell_mass[cell_ooo] = in_rho[cell_ooo] * cell_volumes[cell_ooo];
	out_cell_mass[cell_ooo] = in_cell_mass[cell_ooo];

	in_e[cell_ooo] = EnergyEOS(physical_params.gamma, in_rho[cell_ooo], in_p[cell_ooo], physical_params.pi);
	out_e[cell_ooo] = in_e[cell_ooo];

      } else if (numerical_params.TypeOfProjection == "LagrangeFluxes") {
	
	in_cell_mass[cell_ooo]         = in_rho[cell_ooo] * cell_volumes[cell_ooo];

	// [VM] Gamma*Pi from VM is equivalent to Pi here in Stiffened Gas 
	in_e_1[cell_ooo]               =  EnergyEOS(physical_params.gamma_1, in_rho_1[cell_ooo], in_p_1[cell_ooo], physical_params.pi_1);
	out_e_1[cell_ooo]              = in_e_1[cell_ooo];	
	in_e_2[cell_ooo]               =  EnergyEOS(physical_params.gamma_2, in_rho_2[cell_ooo], in_p_2[cell_ooo], physical_params.pi_2);
	out_e_2[cell_ooo]              = in_e_2[cell_ooo];
	if (nb_mat > 2) {
	  in_e_3[cell_ooo]             =  EnergyEOS(physical_params.gamma_3, in_rho_3[cell_ooo], in_p_3[cell_ooo], physical_params.pi_3);
	}
	out_e_3[cell_ooo]              = in_e_3[cell_ooo];
	
	in_e[cell_ooo]                 = in_e_1[cell_ooo] * in_c_1[cell_ooo] + in_e_2[cell_ooo] * in_c_2[cell_ooo] + in_e_3[cell_ooo] * in_c_3[cell_ooo];
	out_e[cell_ooo]                = in_e[cell_ooo];

	in_p[cell_ooo]                 = in_p_1[cell_ooo];

	// [VM] in VM mk = Alpha_k Rho_k with Alpha_k = in_c_1[cell_ooo] * cell_volumes[cell_ooo]
	in_cell_mass_1[cell_ooo]       = in_rho_1[cell_ooo] * in_c_1[cell_ooo] * cell_volumes[cell_ooo]; 
	out_cell_mass_1[cell_ooo]      = in_cell_mass_1[cell_ooo];
	in_cell_mass_2[cell_ooo]       = in_rho_2[cell_ooo] * in_c_2[cell_ooo] * cell_volumes[cell_ooo];
	out_cell_mass_2[cell_ooo]      = in_cell_mass_2[cell_ooo];
	in_cell_mass_3[cell_ooo]       = in_rho_3[cell_ooo] * in_c_3[cell_ooo] * cell_volumes[cell_ooo];
	out_cell_mass_3[cell_ooo]      = in_cell_mass_3[cell_ooo];

	in_y_1[cell_ooo]               = in_cell_mass_1[cell_ooo] / in_rho[cell_ooo];
	out_y_1[cell_ooo]              = in_y_1[cell_ooo];
	in_y_2[cell_ooo]               = in_cell_mass_2[cell_ooo] / in_rho[cell_ooo];
	out_y_2[cell_ooo]              = in_y_2[cell_ooo];
	in_y_3[cell_ooo]               = in_cell_mass_3[cell_ooo] / in_rho[cell_ooo];
	out_y_3[cell_ooo]              = in_y_3[cell_ooo];

	in_total_energy[cell_ooo]      = 0.5*in_u_cell[cell_ooo]*in_u_cell[cell_ooo] + 0.5*in_v_cell[cell_ooo]*in_v_cell[cell_ooo];
	in_total_energy[cell_ooo]      = in_total_energy[cell_ooo] + in_cell_mass_1[cell_ooo]*in_e_1[cell_ooo]/in_rho[cell_ooo] + in_cell_mass_2[cell_ooo]*in_e_2[cell_ooo]/in_rho[cell_ooo]+ in_cell_mass_3[cell_ooo]*in_e_3[cell_ooo]/in_rho[cell_ooo];
	out_total_energy[cell_ooo]     = in_total_energy[cell_ooo];

	rho_total_energy[cell_ooo]     = in_rho[cell_ooo] * in_total_energy[cell_ooo];

	rho_e[cell_ooo]                = rho_total_energy[cell_ooo] - 0.5*(in_u_cell[cell_ooo]*in_u_cell[cell_ooo] + in_v_cell[cell_ooo]*in_v_cell[cell_ooo])*in_rho[cell_ooo];

	beta[cell_ooo]                 = 1.;
	
	rho_U[cell_ooo]                = in_rho[cell_ooo] * in_u_cell[cell_ooo];               
	rho_V[cell_ooo]                = in_rho[cell_ooo] * in_v_cell[cell_ooo];               
	
	for (int imat = 0; imat < nb_mat; imat++) {
	  alpha_beta_k       [imat][cell_ooo] = in_c_k       [imat][cell_ooo];
	  alpha_beta_k_tmp   [imat][cell_ooo] = alpha_beta_k [imat][cell_ooo];
	  masse_k            [imat][cell_ooo] = in_c_k       [imat][cell_ooo] * in_rho_k[imat][cell_ooo];
	  masse_k_tmp        [imat][cell_ooo] = masse_k      [imat][cell_ooo];
	}
	
      } else if (numerical_params.TypeOfModel == "MultimaterialMix") {
	
	in_c_2[cell_ooo] = 1.0 - in_c_1[cell_ooo];
	out_c_2[cell_ooo] = in_c_2[cell_ooo];

	in_rho[cell_ooo] = in_rho_1[cell_ooo] * in_rho_2[cell_ooo] / (in_rho_1[cell_ooo] * in_c_2[cell_ooo] + in_rho_2[cell_ooo]*in_c_1[cell_ooo]);
	out_rho[cell_ooo] = in_rho[cell_ooo];
	
	in_cell_mass[cell_ooo] = in_rho[cell_ooo] * cell_volumes[cell_ooo];
	
	in_e_1[cell_ooo] =  EnergyEOS(physical_params.gamma_1, in_rho_1[cell_ooo], in_p_1[cell_ooo], physical_params.pi_1);
	out_e_1[cell_ooo] = in_e_1[cell_ooo];	
	in_e_2[cell_ooo] =  EnergyEOS(physical_params.gamma_2, in_rho_2[cell_ooo], in_p_2[cell_ooo], physical_params.pi_2);
	out_e_2[cell_ooo] = in_e_2[cell_ooo];
	
	in_e[cell_ooo] = in_e_1[cell_ooo] * in_c_1[cell_ooo] + in_e_2[cell_ooo] * in_c_2[cell_ooo];
	out_e[cell_ooo] = in_e[cell_ooo];

	in_p[cell_ooo] = in_p_1[cell_ooo];

	
      } else if (numerical_params.TypeOfModel == "MultimaterialInterface") {
	
	const RealType vol_fraction1 = in_cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction2 = 1.0 - vol_fraction1;
	const RealType vol_fraction_1 = vol_fraction1 / (vol_fraction1 + vol_fraction2);
	const RealType vol_fraction_2 = vol_fraction2 / (vol_fraction1 + vol_fraction2);
	
	in_cell_mass_1[cell_ooo] = in_rho_1[cell_ooo] * cell_volumes[cell_ooo] * vol_fraction_1;
	out_cell_mass_1[cell_ooo] = in_cell_mass_1[cell_ooo];
	in_cell_mass_2[cell_ooo] = in_rho_2[cell_ooo] * cell_volumes[cell_ooo] * vol_fraction_2;
	out_cell_mass_2[cell_ooo] = in_cell_mass_2[cell_ooo];
	
	in_cell_mass[cell_ooo] = in_cell_mass_1[cell_ooo] + in_cell_mass_2[cell_ooo];
	out_cell_mass[cell_ooo] = in_cell_mass[cell_ooo];
	
	in_rho[cell_ooo] = in_cell_mass[cell_ooo] / cell_volumes[cell_ooo];
	out_rho[cell_ooo] = in_rho[cell_ooo];
	
	in_e_1[cell_ooo] =  EnergyEOS(physical_params.gamma_1, in_rho_1[cell_ooo], in_p_1[cell_ooo], physical_params.pi_1);
	out_e_1[cell_ooo] = in_e_1[cell_ooo];	
	in_e_2[cell_ooo] =  EnergyEOS(physical_params.gamma_2, in_rho_2[cell_ooo], in_p_2[cell_ooo], physical_params.pi_2);
	out_e_2[cell_ooo] = in_e_2[cell_ooo];

	if (vol_fraction_1 == 0.0) {
	  in_rho_1[cell_ooo] = 0.0;
	  in_e_1[cell_ooo] = 0.0;
	  in_p_1[cell_ooo] = 0.0;
	}

	if (vol_fraction_2 == 0.0) {
	  in_rho_2[cell_ooo] = 0.0;
	  in_e_2[cell_ooo] = 0.0;
	  in_p_2[cell_ooo] = 0.0;
	}

	in_e[cell_ooo] = ( in_e_1[cell_ooo] * in_cell_mass_1[cell_ooo] + in_e_2[cell_ooo] * in_cell_mass_2[cell_ooo] ) / in_cell_mass[cell_ooo];
	out_e[cell_ooo] = in_e[cell_ooo];

	in_p[cell_ooo] = in_p_1[cell_ooo] * vol_fraction_1 + in_p_2[cell_ooo] * vol_fraction_2;
	
      } else {

	in_cell_mass[cell_ooo] = in_rho[cell_ooo] * cell_volumes[cell_ooo];
	out_cell_mass[cell_ooo] = in_cell_mass[cell_ooo];

	in_e[cell_ooo] = EnergyEOS(physical_params.gamma, in_rho[cell_ooo], in_p[cell_ooo], physical_params.pi);
	out_e[cell_ooo] = in_e[cell_ooo];

      }	
            
      cell_pseudo_pressure[cell_ooo] = 0.0;
      predicted_pressure[cell_ooo] = 0.0;
      predicted_pressure_1[cell_ooo] = 0.0;
      predicted_pressure_2[cell_ooo] = 0.0;
      e_lag[cell_ooo] = 0.0;
      e_1_lag[cell_ooo] = 0.0;
      e_2_lag[cell_ooo] = 0.0;
      
    }
  }
 
  for (int iy = 0; iy < ny + 1; ++iy) {
    for (int ix = 0; ix < nx + 1; ++ix) {

      const int node_ooo = ((nx + 1) * iy) +ix;

      u_lag[node_ooo] = 0.0;
      v_lag[node_ooo] = 0.0;
      gradient_u[node_ooo] = 0.0;
      gradient_v[node_ooo] = 0.0;
      in_X_x[node_ooo] = m_grid.xmin() + ix * dx ; 
      in_X_y[node_ooo] = m_grid.xmin() + iy * dy ;
    }
  }

  //For computation of the L2 error
#pragma omp parallel for
  for (index_t cell = 0; cell < nx * ny; cell++) {
    in_rho_0[cell] = in_rho[cell];
    in_cell_volumic_fraction_0[cell] = in_cell_volumic_fraction[cell];
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
  std::cerr << "Type of Model : " << numerical_params.TypeOfModel << "\n";      

  //likwid_markerInit(); //init
  //likwid_markerThreadInit(); //init for deamon access and/or multithread profiling
    index_t ret = 0;

  while (!timetable.IsFinished(clock)) {
    
    if (process_rank == 0)
      std::cerr << "\n"
		<< "Time " << clock.time() << ", iteration " << clock.iter()
		<< "\n";


    /*   
    //FOR TEST ADVECTION RETURN ONLY
    if (clock.time() < 0.2399) {
      for (index_t node = 0; node < nb_nodes; node++) {
	in_u[node] = 5.;
	in_v[node] = 5.;
      }
    } else {
      for (index_t node = 0; node < nb_nodes; node++) {
	in_u[node] = -5.;
	in_v[node] = -5.;
      }
    }      

    for (index_t cell = 0; cell < nb_cells; cell++) {
      in_p[cell] = 1.;
      in_e[cell] = EnergyEOS(physical_params.gamma, in_rho[cell], in_p[cell], physical_params.pi);
    }

    if ((clock.time() > 0.2399) && (ret == 0)) {
      for (index_t node = 0; node < nb_nodes; node++) {
	in_u[node] = - in_u[node];
	in_v[node] = - in_v[node];
	ret = 1;
      }
    }      
   */ 

    
    /*    
    //FOR TEST SOLID ROTATION ONLY
    RealType x0 = 2.0;
    RealType y0 = 2.0;
    RealType omega0 = 2 * PI * 1000 * 100;
    for (index_t ix = 0; ix < nx + 1; ix++) {
      for (index_t iy = 0; iy < ny + 1; iy++) {
	index_t node_ooo = (nx + 1) * iy + ix;
	RealType x = ix * dx;
	RealType y = iy * dy;
	in_u[node_ooo] = - (y - y0) * omega0;
	in_v[node_ooo] = (x - x0) * omega0;
      }
    }

    for (index_t ix = 0; ix < nx; ix++) {
      { index_t iy = 0;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t iy = 1;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t iy = ny - 1;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t iy = ny - 2;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
    }
    for (index_t iy = 0; iy < ny; iy++) {
      { index_t ix = 0;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t ix = 1;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t ix = nx - 1;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
      { index_t ix = nx - 2;
	index_t cell_ooo = nx * iy + ix;
	in_rho[cell_ooo] = 0.1;
      }
    }
    in_u[0] = 0.;
    in_v[0] = 0.;
    //in_rho[0] = 0.1;
    in_u[nx] = 0.;
    in_v[nx] = 0.;
    //in_rho[nx - 1] = 0.1;
    in_u[(nx + 1) * ny] = 0.;
    in_v[(nx + 1) * ny] = 0.;
    //in_rho[nx * (ny - 1)] = 0.1;
    in_u[(nx + 1) * ny + nx] = 0.;
    in_v[(nx + 1) * ny + nx] = 0.;
    //in_rho[nx * ny - 1] = 0.1;

    for (index_t cell = 0; cell < nb_cells; cell++) {
      in_p_2[cell] = 1.0e5*10;
      in_e_2[cell] = EnergyEOS(physical_params.gamma_2, in_rho_2[cell], in_p_2[cell], physical_params.pi_2);
    }
    */       

    
    /*
    // FOR TEST IMPACT ED DEBUG
    RealType M_1_tot = 0.0;
    RealType M_2_tot = 0.0;
    RealType M_tot = 0.0;   
 
    for (index_t cell = 0; cell < nb_cells; cell++) {
      M_1_tot = M_1_tot + in_cell_mass_1[cell];
      M_2_tot = M_2_tot + in_cell_mass_2[cell];
      M_tot = M_tot + in_cell_mass[cell];
    }
    printf("\nM1 = %f,   M2 = %f,   M_tot = %f \n\n",M_1_tot,M_2_tot,M_tot); 
    */

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
	
    struct timespec time1;
    struct timespec time2;

    if (m_grid.dimension() == 2) {

      if (numerical_params.TypeOfProjection == "LagrangeFluxes") {

	dt = TimeStepLFMix(nx, ny, dx, dy, nb_mat, numerical_params.CFL, gamma_mix, pi_prime_mix, speed_of_sound_mix, gamma_k,pi_prime_k, in_rho, in_p, in_c_k, in_u_cell, in_v_cell);
	std::cout << "dt = " << dt << std::endl; 
	// Direct projection 2nd order
	LagrangeFluxes2dDriver(//in
				 numerical_params.BoundaryConditions,
				 numerical_params.TypeOfModel,
				 numerical_params.NumberOfMaterials,
				 nx,
				 ny,
				 nb_cells,
				 dx,
				 dy,
				 dt,
				 halo_width,
				 gamma_mix,
				 speed_of_sound_mix,
				 gamma_k,
				 pi_prime_k,
				 in_c_k,
				 in_rho_k,
				 in_rho,
				 in_cell_volumic_fraction,
				 cell_volumes,
				 // out
				 in_u_cell,
				 in_v_cell,
				 out_u,
				 out_v,
				 out_e,
				 out_e_1,
				 out_e_2,
				 out_rho,
				 out_rho_1,
				 out_rho_2,
				 out_c_k,
				 out_cell_volumic_fraction,
				 density_gradient,
				 density_gradient_y,
				 density_1_gradient,
				 density_1_gradient_y,
				 density_2_gradient,
				 density_2_gradient_y,
				 energy_gradient,
				 energy_gradient_y,
				 energy_1_gradient,
				 energy_1_gradient_y,
				 energy_2_gradient,
				 energy_2_gradient_y,
				 concentration_gradient,
				 concentration_gradient_y,
				 gradient_u,
				 gradient_u_y,
				 gradient_v,
				 gradient_v_y,
			         rho_total_energy,
			         rho_e,
			         rho_U,
			         rho_V,
				 in_p,
				 beta,
				 pi_prime_mix,
				 p_plus_pi_prime,
				 p_plus_pi_prime_gradx,
				 p_plus_pi_prime_grady,
				 rho_e_gradx_left,
				 rho_e_gradx_right,
				 rho_e_grady_top,
				 rho_e_grady_bot,
				 p_gradx_left,
				 p_gradx_right,
				 p_grady_top,
				 p_grady_bot,
				 rho_gradx_left,
				 rho_gradx_right,
				 rho_grady_top,
				 rho_grady_bot,
				 alpha_beta_k,
				 alpha_beta_k_tmp,
				 masse_k,
				 masse_k_tmp,
				 alphak_gradx_left,
				 alphak_gradx_right,
				 alphak_grady_bot,
				 alphak_grady_top);

	std::cout << "Start Lagrange Flux ok ! " << gamma_mix[0] << std::endl;
	exit(0);
	
	std::swap(in_cell_mass, out_cell_mass); 
	std::swap(in_cell_mass_1, out_cell_mass_1);
	std::swap(in_cell_mass_2, out_cell_mass_2);
	std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	std::swap(in_u, out_u); 
	std::swap(in_v, out_v);
	std::swap(in_e, out_e);
	std::swap(in_c_1, out_c_1);
	std::swap(in_c_2, out_c_2);
	std::swap(in_e_1, out_e_1);
	std::swap(in_e_2, out_e_2);
	std::swap(in_rho, out_rho);
	std::swap(in_rho_1, out_rho_1);
	std::swap(in_rho_2, out_rho_2);

      }
      else {

	clock_gettime(CLOCK_REALTIME, &time1);
	
	if ((numerical_params.TypeOfModel == "MultimaterialMix") || (numerical_params.TypeOfModel == "MultimaterialInterface")) {	
	  dt = TimeStepSGPCMultimat(nx, ny, dx, dy, numerical_params.CFL, physical_params.gamma_1, physical_params.gamma_1, physical_params.pi_1, physical_params.pi_1, in_rho_1, in_rho_2, in_p_1, in_p_2, in_cell_volumic_fraction, in_u, in_v);
	  
	} else {
	  dt = TimeStepSGPC(nx, ny, dx, dy, numerical_params.CFL, physical_params.gamma, physical_params.pi, in_rho, in_p, in_u, in_v);
	}
	
	clock_gettime(CLOCK_REALTIME, &time2);
	time_time_step.push_back(diff(time1, time2)); 
	
	std::cerr << "CFL : " << numerical_params.CFL << "\n"; 
	std::cerr << "          dt = " << dt << "\n";
	
	// Lagrange 2D algorithm for compressible Euler equations.
	
	Lagrange2dDriver(//in
			 numerical_params.BoundaryConditions,
			 numerical_params.TypeOfModel,
			 nx,
			 ny,
			 dx,
			 dy,
			 dt,
			 physical_params.gamma,
			 physical_params.gamma_1,
			 physical_params.gamma_2,
			 physical_params.pi,
			 physical_params.pi_1,
			 physical_params.pi_2,
			 in_u,
			 in_v,
			 in_e,
			 in_e_1,
			 in_e_2,
			 in_cell_mass,
			 in_cell_mass_1,
			 in_cell_mass_2,
			 in_rho_1,
			 in_rho_2,
			 in_cell_volumic_fraction,
			 //out
			 in_p,
			 in_p_1,
			 in_p_2,
			 predicted_pressure,
			 predicted_pressure_1,
			 predicted_pressure_2,
			 cell_pseudo_pressure,
			 predicted_u,
			 predicted_v,
			 e_lag,
			 e_1_lag,
			 e_2_lag,
			 u_lag,
			 v_lag,
			 in_X_x,
			 in_X_y,
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
				  numerical_params.TypeOfModel,
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
				  physical_params.gamma_1,
				  physical_params.gamma_2,
				  physical_params.pi_1,
				  physical_params.pi_2,
				  predicted_u,
				  e_lag,
				  e_1_lag,
				  e_2_lag,
				  u_lag,
				  v_lag,
				  in_rho_1,
				  in_rho_2,
				  in_cell_mass,
				  in_cell_mass_1,
				  in_cell_mass_2,
				  in_c_1,
				  in_c_2,
				  in_cell_volumic_fraction,
				  cell_volumes,
				  // out
				  out_u,
				  out_v,
				  out_e,
				  out_cell_mass,
				  out_cell_mass_1,
				  out_cell_mass_2,
				  out_e_1,
				  out_e_2,
				  out_rho,
				  out_rho_1,
				  out_rho_2,
				  out_c_1,
				  out_c_2,
				  out_cell_volumic_fraction,
				  directional_lagrangian_volume,
				  directional_lagrangian_density,
				  directional_lagrangian_density_1,
				  directional_lagrangian_density_2,
				  volume_fluxes_x,
				  volume_fluxes_1_x,
				  volume_fluxes_2_x,
				  mass_flux_x,
				  mass_flux_1_x,
				  mass_flux_2_x,
				  energy_flux_x,
				  energy_flux_1_x,
				  energy_flux_2_x,
				  concentration_flux_x,
				  density_gradient,
				  density_1_gradient,
				  density_2_gradient,
				  energy_gradient,
				  energy_1_gradient,
				  energy_2_gradient,
				  concentration_gradient,
				  gradient_u,
				  gradient_v,
				  interface_normal_x,
				  interface_normal_y,
				  //timing
				  time_compute_volume_fluxes_X,
				  time_gradient_X,
				  time_mass_reconstruct_o2_X,
				  time_project_mass_X,
				  time_reconstruct_energy_o2_X,
				  time_project_energy_X,
				  time_gradient_nodal_X,
				  time_project_nodal_velocity_X,
				  time_periodic_boundary);	
	    
	    std::swap(in_cell_mass, out_cell_mass);
	    std::swap(in_cell_mass_1, out_cell_mass_1);
	    std::swap(in_cell_mass_2, out_cell_mass_2);
	    std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	    std::swap(u_lag, out_u);
	    std::swap(v_lag, out_v);
	    std::swap(e_lag, out_e);
	    std::swap(e_1_lag, out_e_1);
	    std::swap(e_2_lag, out_e_2);
	    std::swap(in_c_1, out_c_1);
	    std::swap(in_c_2, out_c_2);
	    std::swap(in_rho_1, out_rho_1);
	    std::swap(in_rho_2, out_rho_2);
	    
	    // Projection Y
	    
	    AdProjection2dYDriver(//in
				  numerical_params.BoundaryConditions,
				  numerical_params.TypeOfModel,
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
				  physical_params.gamma_1,
				  physical_params.gamma_2,
				  physical_params.pi_1,
				  physical_params.pi_2,
				  predicted_v,
				  e_lag,
				  e_1_lag,
				  e_2_lag,
				  u_lag,
				  v_lag,
				  in_rho_1,
				  in_rho_2,
				  in_cell_mass, 
				  in_cell_mass_1,
				  in_cell_mass_2,
				  in_c_1,
				  in_c_2,
				  in_cell_volumic_fraction,
				  cell_volumes,
				  //out
				  out_u,
				  out_v,
				  out_e,
				  out_cell_mass,
				  out_cell_mass_1,
				  out_cell_mass_2,
				  out_e_1,
				  out_e_2,
				  out_rho,
				  out_rho_1,
				  out_rho_2,
				  out_c_1,
				  out_c_2,
				  out_cell_volumic_fraction,
				  directional_lagrangian_volume_y,
				  directional_lagrangian_density_y,
				  directional_lagrangian_density_1_y,
				  directional_lagrangian_density_2_y,
				  volume_fluxes_y,
				  volume_fluxes_1_y,
				  volume_fluxes_2_y,
				  mass_flux_y,
				  mass_flux_1_y,
				  mass_flux_2_y,
				  energy_flux_y,
				  energy_flux_1_y,
				  energy_flux_2_y,
				  concentration_flux_y,
				  density_gradient_y,
				  density_1_gradient_y,
				  density_2_gradient_y,
				  energy_gradient_y,
				  energy_1_gradient,
				  energy_2_gradient_y,
				  concentration_gradient_y,
				  gradient_u_y,
				  gradient_v_y,
				  interface_normal_x,
				  interface_normal_y,
				  //timing
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
	    std::swap(in_cell_mass_1, out_cell_mass_1);
	    std::swap(in_cell_mass_2, out_cell_mass_2);
	    std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	    std::swap(in_u, out_u); 
	    std::swap(in_v, out_v);
	    std::swap(in_e, out_e);
	    std::swap(in_c_1, out_c_1);
	    std::swap(in_c_2, out_c_2);
	    std::swap(in_e_1, out_e_1);
	    std::swap(in_e_2, out_e_2);
	    std::swap(in_rho, out_rho);
	    std::swap(in_rho_1, out_rho_1);
	    std::swap(in_rho_2, out_rho_2);
	    
	  } else { // !clock.iter()%2  Y then X projection
	    // Projection Y
	    
	    AdProjection2dYDriver(//in
				  numerical_params.BoundaryConditions,
				  numerical_params.TypeOfModel,
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
				  physical_params.gamma_1,
				  physical_params.gamma_2,
				  physical_params.pi_1,
				  physical_params.pi_2,
				  predicted_v,
				  e_lag,
				  e_1_lag,
				  e_2_lag,
				  u_lag,
				  v_lag,
				  in_rho_1,
				  in_rho_2,
				  in_cell_mass, 
				  in_cell_mass_1,
				  in_cell_mass_2,
				  in_c_1,
				  in_c_2,
				  in_cell_volumic_fraction,
				  cell_volumes,
				  //out
				  out_u,
				  out_v,
				  out_e,
				  out_cell_mass,
				  out_cell_mass_1,
				  out_cell_mass_2,
				  out_e_1,
				  out_e_2,
				  out_rho,
				  out_rho_1,
				  out_rho_2,
				  out_c_1,
				  out_c_2,
				  out_cell_volumic_fraction,
				  directional_lagrangian_volume_y,
				  directional_lagrangian_density_y,
				  directional_lagrangian_density_1_y,
				  directional_lagrangian_density_2_y,
				  volume_fluxes_y,
				  volume_fluxes_1_y,
				  volume_fluxes_2_y,
				  mass_flux_y,
				  mass_flux_1_y,
				  mass_flux_2_y,
				  energy_flux_y,
				  energy_flux_1_y,
				  energy_flux_2_y,
				  concentration_flux_y,
				  density_gradient_y,
				  density_1_gradient_y,
				  density_2_gradient_y,
				  energy_gradient_y,
				  energy_1_gradient,
				  energy_2_gradient_y,
				  concentration_gradient_y,
				  gradient_u_y,
				  gradient_v_y,
				  interface_normal_x,
				  interface_normal_y,
				  //timing
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
	    std::swap(in_cell_mass_1, out_cell_mass_1);
	    std::swap(in_cell_mass_2, out_cell_mass_2);
	    std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	    std::swap(u_lag, out_u);
	    std::swap(v_lag, out_v);
	    std::swap(e_lag, out_e);
	    std::swap(e_1_lag, out_e_1);
	    std::swap(e_2_lag, out_e_2);
	    std::swap(in_c_1, out_c_1);
	    std::swap(in_c_2, out_c_2);
	    std::swap(in_rho_1, out_rho_1);
	    std::swap(in_rho_2, out_rho_2);
	    
	    // Projection X
	    
	    AdProjection2dXDriver(//in
				  numerical_params.BoundaryConditions,
				  numerical_params.TypeOfModel,
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
				  physical_params.gamma_1,
				  physical_params.gamma_2,
				  physical_params.pi_1,
				  physical_params.pi_2,
				  predicted_u,
				  e_lag,
				  e_1_lag,
				  e_2_lag,
				  u_lag,
				  v_lag,
				  in_rho_1,
				  in_rho_2,
				  in_cell_mass,
				  in_cell_mass_1,
				  in_cell_mass_2,
				  in_c_1,
				  in_c_2,
				  in_cell_volumic_fraction,
				  cell_volumes,
				  // out
				  out_u,
				  out_v,
				  out_e,
				  out_cell_mass,
				  out_cell_mass_1,
				  out_cell_mass_2,
				  out_e_1,
				  out_e_2,
				  out_rho,
				  out_rho_1,
				  out_rho_2,
				  out_c_1,
				  out_c_2,
				  out_cell_volumic_fraction,
				  directional_lagrangian_volume,
				  directional_lagrangian_density,
				  directional_lagrangian_density_1,
				  directional_lagrangian_density_2,
				  volume_fluxes_x,
				  volume_fluxes_1_x,
				  volume_fluxes_2_x,
				  mass_flux_x,
				  mass_flux_1_x,
				  mass_flux_2_x,
				  energy_flux_x,
				  energy_flux_1_x,
				  energy_flux_2_x,
				  concentration_flux_x,
				  density_gradient,
				  density_1_gradient,
				  density_2_gradient,
				  energy_gradient,
				  energy_1_gradient,
				  energy_2_gradient,
				  concentration_gradient,
				  gradient_u,
				  gradient_v,
				  interface_normal_x,
				  interface_normal_y,
				  //timing
				  time_compute_volume_fluxes_X,
				  time_gradient_X,
				  time_mass_reconstruct_o2_X,
				  time_project_mass_X,
				  time_reconstruct_energy_o2_X,
				  time_project_energy_X,
				  time_gradient_nodal_X,
				  time_project_nodal_velocity_X,
				  time_periodic_boundary);	
	    
	    std::swap(in_cell_mass, out_cell_mass); 
	    std::swap(in_cell_mass_1, out_cell_mass_1);
	    std::swap(in_cell_mass_2, out_cell_mass_2);
	    std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	    std::swap(in_u, out_u); 
	    std::swap(in_v, out_v);
	    std::swap(in_e, out_e);
	    std::swap(in_c_1, out_c_1);
	    std::swap(in_c_2, out_c_2);
	    std::swap(in_e_1, out_e_1);
	    std::swap(in_e_2, out_e_2);
	    std::swap(in_rho, out_rho);
	    std::swap(in_rho_1, out_rho_1);
	    std::swap(in_rho_2, out_rho_2);
	    
	  }
	
	} else if (numerical_params.TypeOfProjection == "DirectProjection") {
	  
	  // Direct projection 2nd order
	  DirectProjection2dDriver(//in
				   numerical_params.BoundaryConditions,
				   numerical_params.TypeOfModel,
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
				   physical_params.gamma_1,
				   physical_params.gamma_2,
				   physical_params.pi_1,
				   physical_params.pi_2,
				   predicted_u,
				   predicted_v,
				   e_lag,
				   e_1_lag,
				   e_2_lag,
				   u_lag,
				   v_lag,
				   in_rho_1,
				   in_rho_2,
				   in_cell_mass,
				   in_cell_mass_1,
				   in_cell_mass_2,
				   in_c_1,
				   in_c_2,
				   in_cell_volumic_fraction,
				   cell_volumes,
				   // out
				   out_u,
				   out_v,
				   out_e,
				   out_cell_mass,
				   out_cell_mass_1,
				   out_cell_mass_2,
				   out_e_1,
				   out_e_2,
				   out_rho,
				   out_rho_1,
				   out_rho_2,
				   out_c_1,
				   out_c_2,
				   out_cell_volumic_fraction,
				   directional_lagrangian_volume,
				   directional_lagrangian_volume_y,
				   directional_lagrangian_density,
				   directional_lagrangian_density_y,
				   directional_lagrangian_density_1,
				   directional_lagrangian_density_1_y,
				   directional_lagrangian_density_2,
				   directional_lagrangian_density_2_y,
				   volume_fluxes_x,
				   volume_fluxes_y,
				   volume_fluxes_1_x,
				   volume_fluxes_1_y,
				   volume_fluxes_2_x,
				   volume_fluxes_2_y,
				   mass_flux_x,
				   mass_flux_y,
				   mass_flux_1_x,
				   mass_flux_1_y,
				   mass_flux_2_x,
				   mass_flux_2_y,
				   energy_flux_x,
				   energy_flux_y,
				   energy_flux_1_x,
				   energy_flux_1_y,
				   energy_flux_2_x,
				   energy_flux_2_y,
				   concentration_flux_x,
				   concentration_flux_y,
				   bool_check_fluxes_x,
				   bool_check_fluxes_y,
				   density_gradient,
				   density_gradient_y,
				   density_1_gradient,
				   density_1_gradient_y,
				   density_2_gradient,
				   density_2_gradient_y,
				   energy_gradient,
				   energy_gradient_y,
				   energy_1_gradient,
				   energy_1_gradient_y,
				   energy_2_gradient,
				   energy_2_gradient_y,
				   concentration_gradient,
				   concentration_gradient_y,
				   gradient_u,
				   gradient_u_y,
				   gradient_v,
				   gradient_v_y,
				   interface_normal_x,
				   interface_normal_y,
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
	  std::swap(in_cell_mass_1, out_cell_mass_1);
	  std::swap(in_cell_mass_2, out_cell_mass_2);
	  std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	  std::swap(in_u, out_u); 
	  std::swap(in_v, out_v);
	  std::swap(in_e, out_e);
	  std::swap(in_c_1, out_c_1);
	  std::swap(in_c_2, out_c_2);
	  std::swap(in_e_1, out_e_1);
	  std::swap(in_e_2, out_e_2);
	  std::swap(in_rho, out_rho);
	  std::swap(in_rho_1, out_rho_1);
	  std::swap(in_rho_2, out_rho_2);
	  
	} else if (numerical_params.TypeOfProjection == "DirectProjectionCornerFluxes") {

	  // Direct projection 2nd order with corner fluxes
	  DirectProjectionCornerFluxes2dDriver(//in
					       numerical_params.BoundaryConditions,
					       numerical_params.TypeOfModel,
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
					       physical_params.gamma,
					       physical_params.gamma_1,
					       physical_params.gamma_2,
					       physical_params.pi,
					       physical_params.pi_1,
					       physical_params.pi_2,
					       predicted_u,
					       predicted_v,
					       e_lag,
					       e_1_lag,
					       e_2_lag,
					       u_lag,
					       v_lag,
					       in_rho_1,
					       in_rho_2,
					       in_cell_mass,
					       in_cell_mass_1,
					       in_cell_mass_2,
					       in_p,
					       in_p_1,
					       in_p_2,
					       in_c_1,
					       in_c_2,
					       in_cell_volumic_fraction,
					       cell_volumes,
					       // out
					       out_u,
					       out_v,
					       out_e,
					       out_cell_mass,
					       out_cell_mass_1,
					       out_cell_mass_2,
					       out_e_1,
					       out_e_2,
					       out_rho,
					       out_rho_1,
					       out_rho_2,
					       out_c_1,
					       out_c_2,
					       out_cell_volumic_fraction,
					       directional_lagrangian_volume,
					       directional_lagrangian_volume_y,
					       directional_lagrangian_density,
					       directional_lagrangian_density_y,
					       directional_lagrangian_density_1,
					       directional_lagrangian_density_1_y,
					       directional_lagrangian_density_2,
					       directional_lagrangian_density_2_y,
					       volume_fluxes_x,
					       volume_fluxes_y,
					       volume_fluxes_corner,
					       volume_fluxes_1_x,
					       volume_fluxes_1_y,
					       volume_fluxes_1_corner,
					       volume_fluxes_2_x,
					       volume_fluxes_2_y,
					       volume_fluxes_2_corner,
					       mass_flux_x,
					       mass_flux_y,
					       mass_flux_1_x,
					       mass_flux_1_y,
					       mass_flux_2_x,
					       mass_flux_2_y,
					       mass_corner_fluxes,
					       mass_1_corner_fluxes,
					       mass_2_corner_fluxes,
					       sign_x_corner_fluxes,
					       sign_y_corner_fluxes,
					       energy_flux_x,
					       energy_flux_y,
					       energy_flux_corner,
					       energy_flux_1_x,
					       energy_flux_1_y,
					       energy_flux_1_corner,
					       energy_flux_2_x,
					       energy_flux_2_y,
					       energy_flux_2_corner,
					       concentration_flux_x,
					       concentration_flux_y,
					       concentration_flux_corner,
					       bool_check_fluxes_x,
					       bool_check_fluxes_y,
					       density_gradient,
					       density_gradient_y,
					       density_gradient_diag,
					       density_gradient_antidiag,
					       density_1_gradient,
					       density_1_gradient_y,
					       density_1_gradient_diag,
					       density_1_gradient_antidiag,
					       density_2_gradient,
					       density_2_gradient_y,
					       density_2_gradient_diag,
					       density_2_gradient_antidiag,
					       energy_gradient,
					       energy_gradient_y,
					       energy_gradient_diag,
					       energy_gradient_antidiag,
					       energy_1_gradient,
					       energy_1_gradient_y,
					       energy_1_gradient_diag,
					       energy_1_gradient_antidiag,
					       energy_2_gradient,
					       energy_2_gradient_y,
					       energy_2_gradient_diag,
					       energy_2_gradient_antidiag,
					       concentration_gradient,
					       concentration_gradient_y,
					       concentration_gradient_diag,
					       concentration_gradient_antidiag,
					       reconstructed_density_faces_x,
					       reconstructed_density_faces_y,
					       reconstructed_density_1_faces_x,
					       reconstructed_density_1_faces_y,
					       reconstructed_density_2_faces_x,
					       reconstructed_density_2_faces_y,
					       reconstructed_energy_faces_x,
					       reconstructed_energy_faces_y,
					       reconstructed_energy_1_faces_x,
					       reconstructed_energy_1_faces_y,
					       reconstructed_energy_2_faces_x,
					       reconstructed_energy_2_faces_y,
					       reconstructed_concentration_faces_x,
					       reconstructed_concentration_faces_y,
					       gradient_u,
					       gradient_u_y,
					       gradient_u_diag,
					       gradient_u_antidiag,
					       gradient_v,
					       gradient_v_y,
					       gradient_v_diag,
					       gradient_v_antidiag,
					       interface_normal_x,
					       interface_normal_y,
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
	  std::swap(in_cell_mass_1, out_cell_mass_1);
	  std::swap(in_cell_mass_2, out_cell_mass_2);
	  std::swap(in_cell_volumic_fraction, out_cell_volumic_fraction);
	  std::swap(in_u, out_u); 
	  std::swap(in_v, out_v);
	  std::swap(in_e, out_e);
	  std::swap(in_c_1, out_c_1);
	  std::swap(in_c_2, out_c_2);
	  std::swap(in_e_1, out_e_1);
	  std::swap(in_e_2, out_e_2);
	  std::swap(in_rho, out_rho);
	  std::swap(in_rho_1, out_rho_1);
	  std::swap(in_rho_2, out_rho_2);
	  
	}
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

      const RealType x0 = 0.5;

      const RealType t = clock.time();
      
      RealType* cell_centers_x = cell_variables(CELL_CENTERS_X);

      RiemannAnalyticalSolver(rho_left, u_left, p_left,
      			      rho_right, u_right, p_right,
      			      physical_params.gamma,
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


  // Computation of the L2 error
  RealType error_L2 = 0.0;
  RealType norm_rho_exact = 0.0;		 
  if (numerical_params.TypeOfModel == "Monomaterial") {
    for (index_t ix = 0; ix < nx; ix++) {
      for (index_t iy = 0; iy < ny; iy++) {      
	index_t cell_ooo = nx * iy + ix;
	error_L2 = error_L2 + (in_rho[cell_ooo] - in_rho_0[cell_ooo]) * (in_rho[cell_ooo] - in_rho_0[cell_ooo]);
	norm_rho_exact = norm_rho_exact + in_rho_0[cell_ooo] * in_rho_0[cell_ooo];      
      }
    }
  } else if (numerical_params.TypeOfModel == "MultimaterialInterface") {
    for (index_t ix = 0; ix < nx; ix++) {
      for (index_t iy = 0; iy < ny; iy++) {      
	index_t cell_ooo = nx * iy + ix;     
	error_L2 = error_L2 + (in_cell_volumic_fraction[cell_ooo] - in_cell_volumic_fraction_0[cell_ooo]) * (in_cell_volumic_fraction[cell_ooo] - in_cell_volumic_fraction_0[cell_ooo]);
	norm_rho_exact = norm_rho_exact + in_cell_volumic_fraction_0[cell_ooo] * in_cell_volumic_fraction_0[cell_ooo];
      }
    }
  }
  error_L2 = error_L2 / norm_rho_exact;
  error_L2 = std::sqrt(error_L2);
  // end Computation of the L2 error


  //likwid_markerClose(); // stop likwid 
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

  std::cerr << "Error L2: " << error_L2 << "\n";

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

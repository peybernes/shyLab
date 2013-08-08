#ifndef SIMULATION_HPP
#define SIMULATION_HP

#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "clock.hpp"
#include "config.hpp"
#include "math_parser.hpp"
#include "mesh.hpp"
#include "numerical_params.hpp"
#include "physical_params.hpp"
#include "structured_grid.hpp"
#include "timetable.hpp"
#include "variable_database.hpp"
#include "variable_store.hpp"

class Event;

typedef std::vector< boost::shared_ptr<Event> > EventCollection;
typedef std::vector<std::string> ExpressionCollection;
typedef std::map<std::string, std::string> ValueExpr_t;

class Simulation {
public:
  // Constructor.
  Simulation();
  /// Initialization. Allocate and initialize variables, read mesh and
  /// other simulation state.
  void Init();
  /// Execution of the main timeloop. Variables updating and I/O.
  void Run();
  /// Cleaning up. Deallocation, etc.
  void Finalize();
  /// Load simulation from parameter file.
  void Load(const std::string&);
  /// Save simulation to parameter file.
  void Save(const std::string&);
  /// Simulation (unstructured triangular) mesh.
  void CreateVariableDatabase();
  Mesh mesh;
  MathParser math_parser;
  Clock clock;  
  PhysicalParams physical_params;
  NumericalParams numerical_params;
  VariableStore cell_variables;
  VariableStore face_variables;
  VariableStore vertice_variables;
  VariableStore boundary_face_variables;
  VariableDatabase variables_database;
  StructuredGrid m_grid;
private:
  ExpressionCollection math_expressions;
  ValueExpr_t secondary_value_expr;
  EventCollection events;
  Timetable timetable;
  // Boundary conditions timeseries for Catalina2 test case.
#ifdef CATALINA2
  std::vector<RealType> bc_time;
  std::vector<RealType> bc_height;
#endif // CATALINA2
  DISALLOW_COPY_AND_ASSIGN(Simulation);
};

#endif // SIMULATION_HPP

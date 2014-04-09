#include "init.hpp"

#include <cassert>
#include <fstream>
#include <omp.h>

#include "array_io.hpp"
#include "serialize.hpp"
#include "simulation.hpp"

InitVariable::InitVariable(Simulation *sim_ptr, const Timetable &timetable):
  Event(sim_ptr, timetable, PRE) {
 
  m_name = "InitVariable";
  m_layout = ROW_MAJOR;

}

int InitVariable::id() {
  return m_variable_entry.id();
}

std::string InitVariable::expression() {
  return m_expression;
}

const std::string InitVariable::expression() const {
  return m_expression;
}

void InitVariable::Execute() {

  // We can not initialise a variable from both a file 
  // And a mathematical expression.
  assert(m_expression.empty() || m_stream_name.empty());

  const int id_variable = m_variable_entry.id();

  std::cerr << "InitVariable::Execute, id=" << id_variable << "\n";

  if (!m_expression.empty()) {
    
    std::cerr << "Initialising variable " << m_var_name 
	      << " with mathematical expression: \""
	      << m_expression << "\"\n";

    // Register the expression with the simulation math. parser.
    m_sim_ptr->math_parser.SetExpr(m_expression);
    
    const VariableSupport support = m_variable_entry.support();
    
    // "x", "y" and "z" variables have different meaning depending if
    // we consider a cell, facet or vertice variable. So we need to
    // reinitialize it.
    m_sim_ptr->math_parser.RemoveVar("x");
    m_sim_ptr->math_parser.RemoveVar("y");
    m_sim_ptr->math_parser.RemoveVar("z");

    if (support == CELL) {

      m_sim_ptr->math_parser.DefineVar("x", &(m_sim_ptr->math_parser.m_local_variables["cell_centers_x"]));
      m_sim_ptr->math_parser.DefineVar("y", &(m_sim_ptr->math_parser.m_local_variables["cell_centers_y"]));
      m_sim_ptr->math_parser.DefineVar("z", &(m_sim_ptr->math_parser.m_local_variables["cell_centers_z"]));

      // Evaluate the expression for each cell.
      for (int i = 0; i < m_sim_ptr->cell_variables.nb_elements(); ++i) {

	m_sim_ptr->math_parser.UpdateReservedVariables(m_sim_ptr, CELL, i);
	m_sim_ptr->cell_variables(id_variable)[i] = m_sim_ptr->math_parser.Eval();

      }

    } else if (support == VERTICE) {

      m_sim_ptr->math_parser.DefineVar("x", &(m_sim_ptr->math_parser.m_local_variables["vertices_x"]));
      m_sim_ptr->math_parser.DefineVar("y", &(m_sim_ptr->math_parser.m_local_variables["vertices_y"]));
      m_sim_ptr->math_parser.DefineVar("z", &(m_sim_ptr->math_parser.m_local_variables["vertices_z"]));

      // Evaluate the expression for each cell.
      for (int i = 0; i < m_sim_ptr->vertice_variables.nb_elements(); ++i) {

	m_sim_ptr->math_parser.UpdateReservedVariables(m_sim_ptr, VERTICE, i);
	m_sim_ptr->vertice_variables(id_variable)[i] = m_sim_ptr->math_parser.Eval();

      }


    } else if (support == FACET) {

      m_sim_ptr->math_parser.DefineVar("x", &(m_sim_ptr->math_parser.m_local_variables["face_centers_x"]));
      m_sim_ptr->math_parser.DefineVar("y", &(m_sim_ptr->math_parser.m_local_variables["face_centers_y"]));
      m_sim_ptr->math_parser.DefineVar("z", &(m_sim_ptr->math_parser.m_local_variables["face_centers_z"]));

      // Evaluate the expression for each cell.
      for (int i = 0; i < m_sim_ptr->face_variables.nb_elements(); ++i) {

	m_sim_ptr->math_parser.UpdateReservedVariables(m_sim_ptr, FACET, i);
	m_sim_ptr->face_variables(id_variable)[i] = m_sim_ptr->math_parser.Eval();

      }

    } else {

      assert(0);

    }
      
  }

  // Read values of a scalar variable from file.
  else if (!m_stream_name.empty()) {

    std::ifstream ifs(m_stream_name.c_str());

    const int nb_cells = m_sim_ptr->cell_variables.nb_elements();
    
    ReadTxtAsciiScalar(nb_cells, &ifs, m_sim_ptr->cell_variables(id_variable));

    if (m_layout == COLUMN_MAJOR) {

      std::cerr << "variable_id=" << id_variable << "\n";
      m_sim_ptr->cell_variables.Transpose(id_variable);

    }

    ifs.close();
  }
}

void InitVariable::Save(ptree &pt) {

  assert(m_expression.empty() || m_stream_name.empty());
  
  ptree pt_timetable;
  m_timetable.Save(pt_timetable);
  
  pt.put("name", "Height"); // FIXME
  pt.put_child("Timetable", pt_timetable);
  pt.put("expression", m_expression);
  pt.put("filename", m_stream_name);

  if (m_layout == ROW_MAJOR) {

    pt.put("layout", "row_major");

  } else if (m_layout == COLUMN_MAJOR) {

    pt.put("layout", "column_major");

  } else {

    assert(0);

  }
  
}

void InitVariable::Load(ptree &pt) {
  
  std::string var_name = pt.get<std::string>("name"); // mandatory

  const bool exists_in_database = 
    (m_sim_ptr->variables_database.count(var_name) > 0);

  if (exists_in_database == false) {

    std::cerr << "ERROR: the variable " << var_name 
	      << " has not been declared in the simulation.\n";
    assert(0);

  }

  m_expression = pt.get<std::string>("expression", "");
  m_stream_name = pt.get<std::string>("filename", "");
  m_var_name = var_name;

  m_variable_entry = m_sim_ptr->variables_database[var_name];
  std::cerr << "InitVariable::Load, variable " << var_name << ", id=" << m_variable_entry.id() << "\n";
  
  assert((!m_expression.empty()) || (!m_stream_name.empty()));

  const std::string layout_name = pt.get<std::string>("layout", "");

  if (layout_name == "row_major") {

    m_layout = ROW_MAJOR;

  } else if (layout_name == "column_major") {

    m_layout = COLUMN_MAJOR;

  } else {

    if (!m_stream_name.empty()) {

      std::cerr << "ERROR: in InitVariable, expected \"row_major\" or \"column_major\", got \"" 
		<< layout_name << "\"\n";
      assert(0);

    }

  }

}



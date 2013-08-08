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

  if (!m_expression.empty()) {

    // Register the expression with the simulation math. parser.
    m_sim_ptr->math_parser.SetExpr(m_expression);
    
    const VariableSupport support = m_variable_entry.support();
    
    if (support == CELL) {

      // Evaluate the expression for each cell.
      for (int i = 0; i < m_sim_ptr->cell_variables.nb_elements(); ++i) {

	m_sim_ptr->math_parser.UpdateReservedVariables(m_sim_ptr, CELL, i);
	m_sim_ptr->cell_variables(id_variable)[i] = m_sim_ptr->math_parser.Eval();

      }

    } else if (support == VERTICE) {

      // Evaluate the expression for each cell.
      for (int i = 0; i < m_sim_ptr->vertice_variables.nb_elements(); ++i) {

	m_sim_ptr->math_parser.UpdateReservedVariables(m_sim_ptr, VERTICE, i);
	m_sim_ptr->vertice_variables(id_variable)[i] = m_sim_ptr->math_parser.Eval();

      }


    } else if (support == FACET) {
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
  
}

void InitVariable::Load(ptree &pt) {
  
  std::string var_name = pt.get<std::string>("name"); // mandatory

  const bool exists_in_database = 
    (m_sim_ptr->variables_database.count(var_name) > 0);

  if (exists_in_database == false) {

    std::cerr << "ERROR: the variable " << var_name 
	      << " has not been declared in the simulation.\n";

  }

  m_expression = pt.get<std::string>("expression", "");
  m_stream_name = pt.get<std::string>("filename", "");
  m_var_name = var_name;
  m_variable_entry = m_sim_ptr->variables_database[var_name];
  
  assert((!m_expression.empty()) || (!m_stream_name.empty()));

}

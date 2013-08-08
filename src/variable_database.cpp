#include "variable_database.hpp"

#include <cassert>
#include <iostream>

#include "variable_metadata.hpp"

VariableEntry::VariableEntry():
  m_id(-1), m_support(NONE),
  m_attribute(PROTECTED |  WRITTEN | COMMUNICATED |  INITIALIZABLE), 
  m_region(TOTAL) {}

VariableEntry::VariableEntry(int id, VariableSupport support,
			     int attribute, VariableRegion region):
  m_id(id), m_support(support),
  m_attribute(attribute), m_region(region) {};

void InitVariableDatabase(VariableDatabase* database_ptr) {

  // FIX ME. We should replace static dispatching by a more
  // dynamical approch; should be done by defining all variables in
  // the parameter file.

  for (int id_variable = 0; 
       id_variable < NB_CELL_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, CELL, 
					cell_variable_attributes[id_variable],
					TOTAL);
  
    (*database_ptr)[cell_variable_names[id_variable]] = entry;
  
    std::cerr << "Stored cell variable " 
	      << cell_variable_names[id_variable]
	      << " in database\n";
    
  }

  for (int id_variable = 0; 
       id_variable < NB_FACE_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, FACET, 
					face_variable_attributes[id_variable],
					TOTAL);

    (*database_ptr)[face_variable_names[id_variable]] = entry;

    std::cerr << "Stored face variable " 
	      << face_variable_names[id_variable]
	      << " in database\n";
  }

  for (int id_variable = 0; 
       id_variable < NB_VERTICE_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, VERTICE, 
					vertice_variable_attributes[id_variable],
					TOTAL);

    std::cerr << "Stored vertice variable " 
	      << vertice_variable_names[id_variable]
	      << " in database\n";

  (*database_ptr)[vertice_variable_names[id_variable]] = entry;
    
  }

}

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

void AddVariableInDatabase(const std::string& variable_name, 
			   const VariableEntry& variable_entry, 
			   VariableDatabase* database_ptr) {

  (*database_ptr)[variable_name] = variable_entry;

}

void InitVariableDatabase(VariableDatabase* database_ptr) {

  for (int id_variable = 0; id_variable < NB_CELL_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, CELL, 
					cell_variable_attributes[id_variable],
					TOTAL);
 
    AddVariableInDatabase(cell_variable_names[id_variable], entry, database_ptr);
 
    // std::cerr << "Stored cell variable " 
    // 	      << cell_variable_names[id_variable]
    // 	      << " in database\n";
    
  }

  for (int id_variable = 0; id_variable < NB_FACE_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, FACET, 
					face_variable_attributes[id_variable],
					TOTAL);

    AddVariableInDatabase(face_variable_names[id_variable], entry, database_ptr);

    // std::cerr << "Stored face variable " 
    // 	      << face_variable_names[id_variable]
    // 	      << " in database\n";
  }

  for (int id_variable = 0; id_variable < NB_VERTICE_VALUES; ++id_variable) {
    
    VariableEntry entry = VariableEntry(id_variable, VERTICE, 
					vertice_variable_attributes[id_variable],
					TOTAL);

    AddVariableInDatabase(vertice_variable_names[id_variable], entry, database_ptr);

    // std::cerr << "Stored vertice variable " 
    // 	      << vertice_variable_names[id_variable]
    // 	      << " in database\n";
    
  }

}

int NumberOfVariableOfGivenSupport(VariableSupport given_support, const VariableDatabase& database) {

  int result = 0;

  for (VariableDatabase::const_iterator it = database.begin(); it != database.end(); ++it) {

    const VariableSupport current_support = it->second.support();

    if (current_support == given_support)
      result += 1;

  }

  return result;
  
}

int NumberOfCellVariables(const VariableDatabase& database) {

  return NumberOfVariableOfGivenSupport(CELL, database);

}

int NumberOfVerticeVariables(const VariableDatabase& database) {

  return NumberOfVariableOfGivenSupport(VERTICE, database);

}

int NumberOfFacetVariables(const VariableDatabase& database) {

  return NumberOfVariableOfGivenSupport(FACET, database);

}

#include "variable_store.hpp"

#include <cassert>
#include <cstdlib>
#include <malloc.h>
#include <cstring> // For memcpy.
#include <iostream>

VariableStore::VariableStore():
  m_nb_variables(0), m_nb_elements(0), m_nx(0), m_ny(0), m_padding(0), m_data(NULL) {};

VariableStore::VariableStore(int nb_variables, int nb_elements, int nx, int ny):
  m_nb_variables(nb_variables), 
  m_nb_elements(nb_elements), 
  m_nx(nx), m_ny(ny),
  m_padding(0), m_data(NULL) {

  Allocate(m_nx, m_ny);

};

VariableStore::VariableStore(int nb_variables, int nb_elements, int nx, int ny, int padding):
  m_nb_variables(nb_variables), 
  m_nb_elements(nb_elements), 
  m_nx(nx), m_ny(ny),
  m_padding(padding), 
  m_data(NULL) {
  
  Allocate(m_nx, m_ny);

};

void VariableStore::Allocate(int nx, int ny) {

  assert(m_padding == 0);

  const int data_size = m_nb_variables * (m_nb_elements + m_padding);

  // alignment : 64 for MIC, 32 for Xeon.
  const int success = posix_memalign((void**)&m_data, 64, data_size * sizeof(RealType));

  assert(success == 0);

  std::cerr << "allocated " << data_size * sizeof(RealType) << " bytes\n";
  
  for (int var = 0; var < m_nb_variables; ++var) {
    //#pragma omp parallel for
    for (int j = 2; j < ny - 2; ++j) {
      for (int i = 0; i < nx - 2; ++i) {

	m_data[(var * (m_nb_elements + m_padding)) + (nx * j) + i] = 0.0;

      }
    }
  }
    
}

void VariableStore::DeAllocate() {

  // if (m_data != NULL)
  //   free(m_data);

  // m_data = NULL;

}

void VariableStore::CopyFrom(const VariableStore& other) {

  DeAllocate();

  m_nb_variables = other.nb_variables();
  m_nb_elements = other.nb_elements();
  m_padding = other.padding();

  Allocate(m_nx, m_ny);

  memcpy(m_data, other.m_data,
	 (m_nb_variables * (m_nb_elements + m_padding)) * sizeof(RealType));


}

RealType* VariableStore::operator()(int id_variable) {

  assert(m_data != NULL);

  RealType* result = NULL;

  if ((0 <= id_variable) && (id_variable < m_nb_variables)) {

    result = m_data + (id_variable * (m_nb_elements + m_padding));

  }

  // If index is invalid, we do not crash; instead, the data pointer
  // is set to NULL.
  else {

    std::cerr << "WARNING: index of variable (" << id_variable
	      << ") greater than number of variable for this store("
	      << m_nb_variables << ")\n";

  }

  assert(result != NULL);

  return result;
}

RealType const * VariableStore::operator()(int id_variable) const {

  assert(m_data != NULL);

  const RealType* result = NULL;

  if ((0 <= id_variable) && (id_variable < m_nb_variables)) {

    result = m_data + (id_variable * (m_nb_elements + m_padding));

  }

  // If index is invalid, we do not crash; instead, the data pointer
  // is set to NULL.
  else {

    std::cerr << "WARNING: index of variable (" << id_variable
	      << ") greater than number of variable for this store("
	      << m_nb_variables << ")\n";

  }

  assert(result != NULL);

  return result;
}

RealType* VariableStore::GetVariable(VariableDatabase& db,
				     std::string& variable_name) {

  const int variable_id = db[variable_name].id();

  return (*this)(variable_id);

}


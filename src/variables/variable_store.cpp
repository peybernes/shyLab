#include "variable_store.hpp"

#include <cassert>
#include <cstdlib>
#include <malloc.h>
#include <cstring> // For memcpy.
#include <iostream>

void VariableStore::Validate() {

  if ((m_nx * m_ny) != m_nb_elements) {

    std::cerr << "ERROR: in this data store, nb_elements=" << m_nb_elements
	      << ", but nx=" << m_nx << " and ny=" << m_ny << "\n";

    assert((m_nx * m_ny) == m_nb_elements);

  }
  
}

void VariableStore::Transpose(int variable_id) {

  assert(m_nb_elements == m_ny * m_nx);

  RealType* temporary_variable = (RealType*)malloc((m_nb_elements + m_padding) * sizeof(RealType));

#pragma omp parallel for
  for (int j = 0; j < m_ny; ++j) {
    for (int i = 0; i < m_nx; ++i) {
    
      temporary_variable[(m_nx * j) + i] = 
      	m_data[(variable_id * (m_nb_elements + m_padding)) + (m_ny * i) + j];

    }
  }

  // A remplacer par un memcpy ??
#pragma omp parallel for
  for (int j = 0; j < m_ny; ++j) {
    for (int i = 0; i < m_nx; ++i) {

      m_data[(variable_id * (m_nb_elements + m_padding)) + (m_nx * j) + i] = 
	temporary_variable[(m_nx * j) + i];

    }
  }

}

VariableStore::VariableStore():
  m_nb_variables(0), m_nb_elements(0), m_nx(0), m_ny(0), m_padding(0), m_data(NULL) {};

VariableStore::VariableStore(int nb_variables, int nb_elements, int nx, int ny):
  m_nb_variables(nb_variables), 
  m_nb_elements(nb_elements), 
  m_nx(nx), m_ny(ny),
  m_padding(0), m_data(NULL) {

  Validate();

  Allocate(m_nx, m_ny);

};

VariableStore::VariableStore(int nb_variables, int nb_elements, int nx, int ny, int padding):
  m_nb_variables(nb_variables), 
  m_nb_elements(nb_elements), 
  m_nx(nx), m_ny(ny),
  m_padding(padding), 
  m_data(NULL) {
  
  Validate();

  Allocate(m_nx, m_ny);

};

void VariableStore::Allocate(int nx, int ny) {

  assert(m_padding == 0);

  const int data_size = m_nb_variables * (m_nb_elements + m_padding);

  // alignment : 64 for MIC, 32 for Xeon.
  //const int success = posix_memalign((void**)&m_data, 32, data_size * sizeof(RealType));
  m_data = (RealType*)malloc(data_size * sizeof(RealType));

  assert(m_data != NULL);

  std::cerr << "Number of elements for variables of this store: " << m_nb_elements 
	    << ", allocated " << data_size * sizeof(RealType) << " bytes\n";
  
  for (int var = 0; var < m_nb_variables; ++var) {
#pragma omp parallel for
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

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

    std::cerr << "ERROR: invalid index (" << id_variable
	      << ") for variable\n";

    assert(0);

  }

  assert(result != NULL);

  return result;
}

RealType const * VariableStore::operator()(int id_variable) const {

  assert(m_data != NULL);
  assert(m_padding == 0);

  const RealType* result = NULL;

  if ((0 <= id_variable) && (id_variable < m_nb_variables)) {

    result = m_data + (id_variable * (m_nb_elements + m_padding));

  } else {

    std::cerr << "ERROR: invalid index (" << id_variable
	      << ") for variable\n";

    assert(0);

  }

  assert(result != NULL);

  return result;
}

RealType* VariableStore::GetVariable(VariableDatabase& db,
				     const std::string& variable_name) {

  std::cerr << "Accessing variable " << variable_name << "\n";
  const int variable_id = db[variable_name].id();

  return (*this)(variable_id);

}


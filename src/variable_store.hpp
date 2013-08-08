#ifndef VARIABLE_STORE_HPP_
#define VARIABLE_STORE_HPP_

#include <string>

//#include <Eigen/Core>

#include "variable_database.hpp"

//typedef Eigen::Array<RealType, Eigen::Dynamic, Eigen::Dynamic> Container2D;
//typedef Eigen::Array<RealType, 1, Eigen::Dynamic> ScalarVariable;

//typedef RealType __attribute__ ((aligned(16)))* Container2D;

/// A very thin wrapper around a linear array representing mesh values
/// stored in a "struct of arrays" fashion.
/// 
class VariableStore {
public:
  VariableStore();
  VariableStore(int nb_variables, int nb_elements, int nx, int ny);
  VariableStore(int nb_variables, int nb_elements, int nx, int ny, int padding);
  void CopyFrom(const VariableStore& other);
  RealType* operator()(int id_variable);
  RealType const * operator()(int id_variable) const;
  inline int nb_variables() {return m_nb_variables;}
  inline int nb_variables() const {return m_nb_variables;}
  inline int nb_elements() {return m_nb_elements;}
  inline int nb_elements() const {return m_nb_elements;}
  inline int padding() {return m_padding;}
  inline int padding() const {return m_padding;}
  void DeAllocate();
  RealType* GetVariable(VariableDatabase& db,
			std::string& variable_name);
protected:
  RealType* m_data;
private:
  void Allocate(int nx, int ny);
  int m_nb_variables;
  int m_nb_elements;
  int m_nx;
  int m_ny;
  int m_padding;
};

#endif // VARIABLE_STORE_HPP_

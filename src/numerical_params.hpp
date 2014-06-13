#ifndef NUMERICALPARAMS_HPP
#define NUMERICALPARAMS_HPP

#include "serialize_fwd.hpp"

// template defining a small value (used for testing for dry cells) at
// compile time. User can overwrite this value.
template<typename T>
struct EPSILON {
  static T value() {return static_cast<T>(0);}
};
// EPSILON specialization for floats
template<>
struct EPSILON<float> {
  static float value() {return 1e-5;}
};
// EPSILON specialization for doubles
template<>
struct EPSILON<double> {
  static double value() {return 1e-10;}
};

RealType const EPS = EPSILON<RealType>::value();

struct NumericalParams {
  //! CFL stability factor.
  RealType CFL;
  //! Artificial viscosity for Roe type numerical fluxes.
  RealType beta;
  //! Value of water height for dry cells.
  RealType dry_height;
  //! Type of Boundary Conditions
  std::string BoundaryConditions;
  //! Type of Projection
  std::string TypeOfProjection;
  NumericalParams();
  void save(ptree &pt);
  void load(ptree &pt);
};

#endif // NUMERICALPARAMS_HPP

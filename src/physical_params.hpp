#ifndef PHYSICALPARAMS_HPP
#define PHYSICALPARAMS_HPP

#include "config.hpp"
#include "serialize_fwd.hpp"

struct PhysicalParams {
  RealType g;
  RealType gamma;
  RealType gamma_1;
  RealType gamma_2;
  RealType pi;
  RealType pi_1;
  RealType pi_2;
  PhysicalParams();
  void save(ptree &pt);
  void load(ptree &pt);
};

#endif // PHYSICALPARAMS_HPP

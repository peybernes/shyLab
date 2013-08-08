#ifndef PHYSICALPARAMS_HPP
#define PHYSICALPARAMS_HPP

#include "config.hpp"
#include "serialize_fwd.hpp"

struct PhysicalParams {
  RealType g;
  PhysicalParams();
  void save(ptree &pt);
  void load(ptree &pt);
};

#endif // PHYSICALPARAMS_HPP

#include "physical_params.hpp"
#include "serialize.hpp"

PhysicalParams::PhysicalParams():
      g(9.81) {}

void PhysicalParams::save(ptree &pt) {
  pt.put("gravity", g);
}

void PhysicalParams::load(ptree &pt) {
  g = pt.get<RealType>("gravity", 9.81);
}

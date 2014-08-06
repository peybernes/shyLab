#include <cassert>

#include "numerical_params.hpp"
#include "serialize.hpp"

NumericalParams::NumericalParams():
  CFL(0.34), beta(0.0), dry_height(EPS), BoundaryConditions("Wall"), TypeOfProjection("AdProjection"), TypeOfModel("Monomaterial") {}

void NumericalParams::save(ptree &pt) {
  pt.put("CFL", CFL);
  pt.put("beta", beta);
  pt.put("dryHeight", dry_height);
  pt.put("BoundaryConditions", "Wall");
  pt.put("TypeOfProjection", "AdProjection");
  pt.put("TypeOfModel", "Monomaterial");
}

void NumericalParams::load(ptree &pt) {
  
  CFL = pt.get<RealType>("CFL", 0.34);
  beta = pt.get<RealType>("beta", 0.0);
  dry_height = pt.get<RealType>("dryHeight", EPS);
  BoundaryConditions = pt.get<std::string>("BoundaryConditions", "Wall");
  TypeOfProjection = pt.get<std::string>("TypeOfProjection", "AdProjection");
  TypeOfModel = pt.get<std::string>("TypeOfModel", "Monomaterial");
}


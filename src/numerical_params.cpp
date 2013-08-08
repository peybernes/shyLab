#include "numerical_params.hpp"
#include "serialize.hpp"

NumericalParams::NumericalParams():
  CFL(0.40), beta(0.0), dry_height(EPS) {}

void NumericalParams::save(ptree &pt) {
  pt.put("CFL", CFL);
  pt.put("beta", beta);
  pt.put("dryHeight", dry_height);
}

void NumericalParams::load(ptree &pt) {
  pt.get<RealType>("CFL", 0.49);
  pt.get<RealType>("beta", 0.0);
  pt.get<RealType>("dryHeight", EPS);
}


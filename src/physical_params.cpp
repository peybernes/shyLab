#include "physical_params.hpp"
#include "serialize.hpp"

PhysicalParams::PhysicalParams():
  g(9.81), gamma(1.4), gamma_1(1.4), gamma_2(1.4), gamma_3(1.4), gamma_4(1.4), pi(0.0), pi_1(0.0), pi_2(0.0), pi_3(0.0), pi_4(0.0) {}
void PhysicalParams::save(ptree &pt) {
  pt.put("gravity", g);
  pt.put("EOS_gamma", gamma);
  pt.put("EOS_gamma_1", gamma_1);
  pt.put("EOS_gamma_2", gamma_2);
  pt.put("EOS_gamma_3", gamma_3);
  pt.put("EOS_gamma_4", gamma_4);
  pt.put("EOS_pi", pi);
  pt.put("EOS_pi_1", pi_1);
  pt.put("EOS_pi_2", pi_2);
  pt.put("EOS_pi_3", pi_3);
  pt.put("EOS_pi_4", pi_4);
}

void PhysicalParams::load(ptree &pt) {
  g = pt.get<RealType>("gravity", 9.81);
  gamma = pt.get<RealType>("EOS_gamma", 1.4);
  gamma_1 = pt.get<RealType>("EOS_gamma_1", 1.4);
  gamma_2 = pt.get<RealType>("EOS_gamma_2", 1.4);
  gamma_3 = pt.get<RealType>("EOS_gamma_3", 1.4);
  gamma_4 = pt.get<RealType>("EOS_gamma_4", 1.4);
  pi = pt.get<RealType>("EOS_pi", 0.0);
  pi_1 = pt.get<RealType>("EOS_pi_1", 0.0);
  pi_2 = pt.get<RealType>("EOS_pi_2", 0.0);
  pi_3 = pt.get<RealType>("EOS_pi_3", 0.0);
  pi_4 = pt.get<RealType>("EOS_pi_4", 0.0);
}

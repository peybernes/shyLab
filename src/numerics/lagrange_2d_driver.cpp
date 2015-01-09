// -*- c++ -*- (for emacs users)
#include "lagrange_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "kernel_tools.h"

void Lagrange2dDriver(const std::string BoundaryConditions,
		      const std::string TypeOfModel,
		      const int nx,
		      const int ny,
		      const RealType dx,
		      const RealType dy,
		      const RealType dt,
		      const RealType gamma,
		      const RealType gamma_1,
		      const RealType gamma_2,
		      const RealType pi,
		      const RealType pi_1,
		      const RealType pi_2,
		      const RealType* RESTRICT in_u,
		      const RealType* RESTRICT in_v,
		      const RealType* RESTRICT in_e,
		      const RealType* RESTRICT in_e_1,
		      const RealType* RESTRICT in_e_2,
		      const RealType* RESTRICT in_cell_mass,
		      const RealType* RESTRICT in_cell_mass_1,
		      const RealType* RESTRICT in_cell_mass_2,
		      const RealType* RESTRICT in_rho_1,
		      const RealType* RESTRICT in_rho_2,
		      const RealType* RESTRICT in_cell_volumic_fraction,
		      RealType* RESTRICT in_p,
		      RealType* RESTRICT in_p_1,
		      RealType* RESTRICT in_p_2,
		      RealType* RESTRICT predicted_pressure,
		      RealType* RESTRICT predicted_pressure_1,
		      RealType* RESTRICT predicted_pressure_2,
		      RealType* RESTRICT cell_pseudo_pressure,
		      RealType* RESTRICT predicted_u,
		      RealType* RESTRICT predicted_v,
		      RealType* RESTRICT e_lag,
		      RealType* RESTRICT e_1_lag,
		      RealType* RESTRICT e_2_lag,
		      RealType* RESTRICT u_lag,
		      RealType* RESTRICT v_lag,
		      RealType* RESTRICT in_X_x,
		      RealType* RESTRICT in_X_y,
		      //timing
		      std::vector<RealType>& time_lagrange_pressure_predicted,
		      std::vector<RealType>& time_lagrange_velocity_predicted,
		      std::vector<RealType>& time_periodic_boundary,
		      std::vector<RealType>& time_lagrange_correction,
		      std::vector<RealType>& time_lagrange_velocity_correction) {
  

  if (TypeOfModel == "MultimaterialInterface") {

    Lagrange2dMultimatInterfaceDriver(BoundaryConditions, nx, ny, dx, dy, dt, gamma_1, gamma_2, pi_1, pi_2, in_u, in_v, in_e, in_e_1, in_e_2, in_cell_mass, in_cell_mass_1, in_cell_mass_2, in_rho_1, in_rho_2, in_cell_volumic_fraction, in_p, in_p_1, in_p_2, predicted_pressure, predicted_pressure_1, predicted_pressure_2, cell_pseudo_pressure, predicted_u, predicted_v, e_lag, e_1_lag, e_2_lag, u_lag, v_lag, in_X_x, in_X_y, time_lagrange_pressure_predicted, time_lagrange_velocity_predicted, time_periodic_boundary, time_lagrange_correction, time_lagrange_velocity_correction);

  } else {

    Lagrange2dDefaultDriver(BoundaryConditions, nx, ny, dx, dy, dt, gamma, pi, in_u, in_v, in_e, in_cell_mass, in_p, predicted_pressure, cell_pseudo_pressure, predicted_u, predicted_v, e_lag, u_lag, v_lag, in_X_x, in_X_y, time_lagrange_pressure_predicted, time_lagrange_velocity_predicted, time_periodic_boundary, time_lagrange_correction, time_lagrange_velocity_correction);

  }

}
    

#ifndef LAGRANGE_2D_INTERFACE_H
#define LAGRANGE_2D_INTERFACE_H
#include <vector>
#include <string>

#include "lagrange_2d_default_driver.hpp"
#include "lagrange_2d_multimat_interface_driver.hpp"


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
		      const RealType* in_u,
		      const RealType* in_v,
		      const RealType* in_e,
		      const RealType* in_e_1,
		      const RealType* in_e_2,
		      const RealType* in_cell_mass,
		      const RealType* in_cell_mass_1,
		      const RealType* in_cell_mass_2,
		      const RealType* in_rho_1,
		      const RealType* in_rho_2,
		      const RealType* in_cell_volumic_fraction,
		      RealType* in_p,
		      RealType* in_p_1,
		      RealType* in_p_2,
		      RealType* predicted_pressure,
		      RealType* predicted_pressure_1,
		      RealType* predicted_pressure_2,
		      RealType* cell_pseudo_pressure,
		      RealType* predicted_u,
		      RealType* predicted_v,
		      RealType* e_lag,
		      RealType* e_1_lag,
		      RealType* e_2_lag,
		      RealType* u_lag,
		      RealType* v_lag,
		      RealType* in_X_x,
		      RealType* in_X_y,
		      //timing
		      std::vector<RealType>& time_lagrange_pressure_predicted,
		      std::vector<RealType>& time_lagrange_velocity_predicted,
		      std::vector<RealType>& time_periodic_boundary,
		      std::vector<RealType>& time_lagrange_correction,
		      std::vector<RealType>& time_lagrange_velocity_correction);



#endif 

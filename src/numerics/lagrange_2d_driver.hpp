#ifndef LAGRANGE_2D_INTERFACE_H
#define LAGRANGE_2D_INTERFACE_H
#include <vector>
#include <string>

#include "kernel_lagrange_2d.h"
#include "kernel_ad_projection_2d_boundary_interface.hpp"

void Lagrange2dDriver(const std::string BoundaryConditions,
		      const int nx,
		      const int ny,
		      const RealType dx,
		      const RealType dy,
		      const RealType dt,
		      const RealType* in_u,
		      const RealType* in_v,
		      const RealType* in_e,
		      const RealType* in_cell_mass,
		      RealType* in_p,
		      RealType* predicted_pressure,
		      RealType* cell_pseudo_pressure,
		      RealType* predicted_u,
		      RealType* predicted_v,
		      RealType* e_lag,
		      RealType* u_lag,
		      RealType* v_lag,
		      //timing
		      std::vector<RealType>& time_lagrange_pressure_predicted,
		      std::vector<RealType>& time_lagrange_velocity_predicted,
		      std::vector<RealType>& time_periodic_boundary,
		      std::vector<RealType>& time_lagrange_correction,
		      std::vector<RealType>& time_lagrange_velocity_correction);

#endif 

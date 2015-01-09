// -*- c++ -*- (for emacs users)
#include "lagrange_2d_multimat_interface_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdio.h>

#include "kernel_tools.h"


void Lagrange2dMultimatInterfaceDriver(const std::string BoundaryConditions,
				       const int nx,
				       const int ny,
				       const RealType dx,
				       const RealType dy,
				       const RealType dt,
				       const RealType gamma_1,
				       const RealType gamma_2,
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

  //for timing
  struct timespec time1;
  struct timespec time2;
  //

  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangePressurePredictedOptimisedMultimat(nx, ny, dt, dx, dy, gamma_1, gamma_2, pi_1, pi_2, in_cell_mass_1, in_cell_mass_2, in_rho_1, in_rho_2, in_e_1, in_e_2, in_cell_volumic_fraction, in_u, in_v,
					     in_p, in_p_1, in_p_2, predicted_pressure, predicted_pressure_1, predicted_pressure_2, cell_pseudo_pressure);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_pressure_predicted.push_back(diff(time1, time2)); 

  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeVelocityPredicted(nx, ny, dt, dx, dy, in_cell_mass, predicted_pressure, cell_pseudo_pressure, in_u, in_v,
			    predicted_u , predicted_v);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_velocity_predicted.push_back(diff(time1, time2)); 
    
  clock_gettime(CLOCK_REALTIME, &time1);
  BoundaryVelocityPrediction(BoundaryConditions, nx, ny, dt, dx, dy, in_cell_mass, predicted_pressure,cell_pseudo_pressure, in_u, in_v,
			     predicted_u , predicted_v);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_periodic_boundary.push_back(diff(time1, time2));   
    
  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeCorrectionOptimisedMultimat(nx, ny, dt, dx, dy, in_cell_mass_1, in_cell_mass_2, in_e_1, in_e_2, predicted_pressure_1, predicted_pressure_2, cell_pseudo_pressure, in_cell_volumic_fraction, predicted_u, predicted_v,
				      e_lag, e_1_lag, e_2_lag);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_correction.push_back(diff(time1, time2)); 
    
  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeVelocityCorrection(nx, ny, in_u, in_v, predicted_u, predicted_v,
			     u_lag, v_lag);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_velocity_correction.push_back(diff(time1, time2)); 


  
}

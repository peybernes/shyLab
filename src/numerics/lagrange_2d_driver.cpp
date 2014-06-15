// -*- c++ -*- (for emacs users)
#include "lagrange_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "kernel_tools.h"

void Lagrange2dDriver(const std::string BoundaryConditions,
		      const int nx,
		      const int ny,
		      const RealType dx,
		      const RealType dy,
		      const RealType dt,
		      const RealType* RESTRICT in_u,
		      const RealType* RESTRICT in_v,
		      const RealType* RESTRICT in_e,
		      const RealType* RESTRICT in_cell_mass,
		      RealType* RESTRICT in_p,
		      RealType* RESTRICT predicted_pressure,
		      RealType* RESTRICT cell_pseudo_pressure,
		      RealType* RESTRICT predicted_u,
		      RealType* RESTRICT predicted_v,
		      RealType* RESTRICT e_lag,
		      RealType* RESTRICT u_lag,
		      RealType* RESTRICT v_lag,
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
  LagrangePressurePredictedOptimised(nx, ny, dt, dx, dy, in_cell_mass, in_e, in_u, in_v,
				     in_p , predicted_pressure, cell_pseudo_pressure);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_pressure_predicted.push_back(diff(time1, time2)); 
    
  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeVelocityPredicted(nx, ny, dt, dx, dy, in_cell_mass, predicted_pressure,cell_pseudo_pressure, in_u, in_v,
			    predicted_u , predicted_v);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_velocity_predicted.push_back(diff(time1, time2)); 
    
  clock_gettime(CLOCK_REALTIME, &time1);
  BoundaryVelocityPrediction(BoundaryConditions, nx, ny, dt, dx, dy, in_cell_mass, predicted_pressure,cell_pseudo_pressure, in_u, in_v,
			     predicted_u , predicted_v);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_periodic_boundary.push_back(diff(time1, time2));   
    
  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeCorrection(nx, ny, dt, dx, dy, in_cell_mass, in_e, predicted_pressure,cell_pseudo_pressure, predicted_u, predicted_v,
		     e_lag);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_correction.push_back(diff(time1, time2)); 
    
  clock_gettime(CLOCK_REALTIME, &time1);
  LagrangeVelocityCorrection(nx, ny, in_u, in_v, predicted_u, predicted_v,
			     u_lag, v_lag);
  clock_gettime(CLOCK_REALTIME, &time2);
  time_lagrange_velocity_correction.push_back(diff(time1, time2)); 
  
}

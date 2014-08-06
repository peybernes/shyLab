// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_2d_boundary_interface.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"



void ProjectNodalIntensiveVariableOrder2BoundaryDirect(std::string BoundaryConditions,
						  std::string variable_projected,
						  index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable_x,
						  const RealType* gradient_variable_y,
						  const RealType* mass_flux_x,
						  const RealType* mass_flux_y,
						  RealType* out_variable) {

  
  if (BoundaryConditions == "periodic") {
    ProjectNodalIntensiveVariableOrder2PeriodicBoundaryDirect(nx, ny, halo_width, dx, dy, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, v_velocity_pred, in_variable, gradient_variable_x, gradient_variable_y, mass_flux_x, mass_flux_y, out_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ProjectNodalIntensiveVariableOrder2BoundaryRtDirect(nx, ny, halo_width, dx, dy, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, v_velocity_pred, in_variable, gradient_variable_x, gradient_variable_y, mass_flux_x, mass_flux_y, out_variable);
  }
  else if (BoundaryConditions == "Wall") {
    if (variable_projected == "project_ux") {
      ProjectNodalIntensiveVariableUxWallBoundaryDirect(nx, ny, halo_width, lag_cell_mass, out_cell_mass, in_variable, mass_flux_x, mass_flux_y, out_variable);
    }
    else if (variable_projected == "project_uy") {
      ProjectNodalIntensiveVariableUyWallBoundaryDirect(nx, ny, halo_width, lag_cell_mass, out_cell_mass, in_variable, mass_flux_x, mass_flux_y,  out_variable);
    }
  }
  else{
    ProjectNodalIntensiveVariableOrder2PeriodicBoundaryDirect(nx, ny, halo_width, dx, dy, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, v_velocity_pred, in_variable, gradient_variable_x, gradient_variable_y, mass_flux_x, mass_flux_y, out_variable); 
  }


}

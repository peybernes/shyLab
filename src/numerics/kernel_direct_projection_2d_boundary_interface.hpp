#ifndef KERNEL_DIRECT_PROJECTION_2D_BOUNDARY_INTERFACE_H
#define KERNEL_DIRECT_PROJECTION_2D_BOUNDARY_INTERFACE_H
#include <string>
typedef int index_t;

#include "kernel_direct_projection_2d_boundary_conditions.h"
#include "kernel_lagrange_2d.h"


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
						  RealType* out_variable);



#endif 

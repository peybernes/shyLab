#ifndef DIRECT_PROJECTION_2D_MONOMAT_INTERFACE_H
#define DIRECT_PROJECTION_2D_MONOMAT_INTERFACE_H
#include <vector>
#include <string>

#include "kernels.hpp"

void DirectProjection2dMonomaterialDriver(//in
					  const std::string BoundaryConditions,
					  const int nx,
					  const int ny,
					  const RealType dx,
					  const RealType dy,
					  const RealType dt,
					  const RealType halo_width,
					  const RealType* predicted_u,
					  const RealType* predicted_v,
					  const RealType* cell_volumes,
					  const RealType* e_lag,
					  const RealType* u_lag,
					  const RealType* v_lag,
					  const RealType* in_cell_mass,
					  //out
					  RealType* out_u,
					  RealType* out_v,
					  RealType* out_e,
					  RealType* out_cell_mass,
					  RealType* out_rho,
					  RealType* directional_lagrangian_volume,
					  RealType* directional_lagrangian_volume_y,
					  RealType* directional_lagrangian_density,
					  RealType* directional_lagrangian_density_y,
					  RealType* volume_fluxes_x,
					  RealType* volume_fluxes_y,
					  RealType* mass_flux_x,
					  RealType* mass_flux_y,
					  RealType* energy_flux_x,
					  RealType* energy_flux_y,
					  RealType* density_gradient,
					  RealType* density_gradient_y,
					  RealType* energy_gradient,
					  RealType* energy_gradient_y,
					  RealType* gradient_u,
					  RealType* gradient_u_y,
					  RealType* gradient_v,
					  RealType* gradient_v_y,
					  //timing
					  std::vector<RealType>& time_compute_volume_fluxes_X,
					  std::vector<RealType>& time_gradient_X,
					  std::vector<RealType>& time_mass_reconstruct_o2_X,
					  std::vector<RealType>& time_project_mass_X,
					  std::vector<RealType>& time_reconstruct_energy_o2_X,
					  std::vector<RealType>& time_project_energy_X,
					  std::vector<RealType>& time_gradient_nodal_X,
					  std::vector<RealType>& time_project_nodal_velocity_X,
					  std::vector<RealType>& time_compute_volume_fluxes_Y,
					  std::vector<RealType>& time_gradient_Y,
					  std::vector<RealType>& time_mass_reconstruct_o2_Y,
					  std::vector<RealType>& time_project_mass_Y,
					  std::vector<RealType>& time_reconstruct_energy_o2_Y,
					  std::vector<RealType>& time_project_energy_Y,
					  std::vector<RealType>& time_gradient_nodal_Y,
					  std::vector<RealType>& time_project_nodal_velocity_Y,
					  std::vector<RealType> time_periodic_boundary);



#endif 

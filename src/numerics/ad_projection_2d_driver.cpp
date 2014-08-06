// -*- c++ -*- (for emacs users)
#include "ad_projection_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"



void AdProjection2dXDriver(//in
			   const std::string BoundaryConditions,
			   const std::string TypeOfModel,
			   const int nx,
			   const int ny,
			   const int nb_faces_x,
			   const int nb_faces_y,
			   const int nb_cells,
			   const int nb_nodes,
			   const RealType dx,
			   const RealType dy,
			   const RealType dt,
			   const RealType halo_width,
			   const RealType gamma_1,
			   const RealType gamma_2,
			   const RealType pi_1,
			   const RealType pi_2,
			   const RealType* RESTRICT predicted_u,
			   const RealType* RESTRICT e_lag,
			   const RealType* RESTRICT e_1_lag,
			   const RealType* RESTRICT e_2_lag,
			   const RealType* RESTRICT u_lag,
			   const RealType* RESTRICT v_lag,
			   const RealType* RESTRICT in_rho_1,
			   const RealType* RESTRICT in_rho_2,
			   const RealType* RESTRICT in_cell_mass,
			   const RealType* RESTRICT in_cell_mass_1,
			   const RealType* RESTRICT in_cell_mass_2,
			   const RealType* RESTRICT in_c_1,
			   const RealType* RESTRICT in_c_2,
			   const RealType* RESTRICT in_cell_volumic_fraction,
			   const RealType* RESTRICT cell_volumes,
			   // out
			   RealType* RESTRICT out_u,
			   RealType* RESTRICT out_v,
			   RealType* RESTRICT out_e,
			   RealType* RESTRICT out_cell_mass,
			   RealType* RESTRICT out_cell_mass_1,
			   RealType* RESTRICT out_cell_mass_2,
			   RealType* RESTRICT out_e_1,
			   RealType* RESTRICT out_e_2,
			   RealType* RESTRICT out_rho,
			   RealType* RESTRICT out_rho_1,
			   RealType* RESTRICT out_rho_2,
			   RealType* RESTRICT out_c_1,
			   RealType* RESTRICT out_c_2,
			   RealType* RESTRICT out_cell_volumic_fraction,
			   RealType* RESTRICT directional_lagrangian_volume,
			   RealType* RESTRICT directional_lagrangian_density,
			   RealType* RESTRICT directional_lagrangian_density_1,
			   RealType* RESTRICT directional_lagrangian_density_2,
			   RealType* RESTRICT volume_fluxes_x,
			   RealType* RESTRICT volume_fluxes_1_x,
			   RealType* RESTRICT volume_fluxes_2_x,
			   RealType* RESTRICT mass_flux_x,
			   RealType* RESTRICT mass_flux_1_x,
			   RealType* RESTRICT mass_flux_2_x,
			   RealType* RESTRICT energy_flux_x,
			   RealType* RESTRICT energy_flux_1_x,
			   RealType* RESTRICT energy_flux_2_x,
			   RealType* RESTRICT concentration_flux_x,
			   RealType* RESTRICT density_gradient,
			   RealType* RESTRICT density_1_gradient,
			   RealType* RESTRICT density_2_gradient,
			   RealType* RESTRICT energy_gradient,
			   RealType* RESTRICT energy_1_gradient,
			   RealType* RESTRICT energy_2_gradient,
			   RealType* RESTRICT concentration_gradient,
			   RealType* RESTRICT gradient_u,
			   RealType* RESTRICT gradient_v,
			   RealType* RESTRICT interface_normal_x,
			   RealType* RESTRICT interface_normal_y,
			   //timing
			   std::vector<RealType>& time_compute_volume_fluxes_X,
			   std::vector<RealType>& time_gradient_X,
			   std::vector<RealType>& time_mass_reconstruct_o2_X,
			   std::vector<RealType>& time_project_mass_X,
			   std::vector<RealType>& time_reconstruct_energy_o2_X,
			   std::vector<RealType>& time_project_energy_X,
			   std::vector<RealType>& time_gradient_nodal_X,
			   std::vector<RealType>& time_project_nodal_velocity_X,
			   std::vector<RealType> time_periodic_boundary) {


  if (TypeOfModel == "Momomaterial") {
    
    AdProjection2dMonomaterialXDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, predicted_u, e_lag, u_lag, v_lag, in_cell_mass, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, mass_flux_x, energy_flux_x, density_gradient, energy_gradient, gradient_u, gradient_v, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_periodic_boundary);

  } else if (TypeOfModel == "MultimaterialMix") {
    
    AdProjection2dMultimaterialMixXDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, gamma_1, gamma_2, pi_1, pi_2, predicted_u, e_lag, u_lag, v_lag, in_cell_mass, in_c_1, in_c_2, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_e_1,out_e_2, out_rho, out_rho_1, out_rho_2, out_c_1, out_c_2, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, mass_flux_x, energy_flux_x, concentration_flux_x, density_gradient, energy_gradient, concentration_gradient, gradient_u, gradient_v,time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_periodic_boundary);

  } else if (TypeOfModel == "MultimaterialInterface") {

    AdProjection2dMultimaterialInterfaceXDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, gamma_1, gamma_2, pi_1, pi_2, predicted_u, e_lag,e_1_lag, e_2_lag, u_lag, v_lag, in_rho_1, in_rho_2, in_cell_mass, in_cell_mass_1,in_cell_mass_2, in_cell_volumic_fraction, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_cell_mass_1, out_cell_mass_2, out_e_1, out_e_2, out_rho, out_rho_1, out_rho_2, out_cell_volumic_fraction, directional_lagrangian_volume, directional_lagrangian_density, directional_lagrangian_density_1, directional_lagrangian_density_2, volume_fluxes_x, volume_fluxes_1_x, volume_fluxes_2_x, mass_flux_x, mass_flux_1_x, mass_flux_2_x, energy_flux_x, energy_flux_1_x, energy_flux_2_x, density_gradient, density_1_gradient, density_2_gradient, energy_gradient, energy_1_gradient, energy_2_gradient, gradient_u, gradient_v, interface_normal_x, interface_normal_y, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X,time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_periodic_boundary);

  } else {
    
    AdProjection2dMonomaterialXDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, predicted_u, e_lag, u_lag, v_lag, in_cell_mass, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, mass_flux_x, energy_flux_x, density_gradient, energy_gradient, gradient_u, gradient_v, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_periodic_boundary);

  }
  
}



void AdProjection2dYDriver(//in
			   const std::string BoundaryConditions,
			   const std::string TypeOfModel,
			   const int nx,
			   const int ny,
			   const int nb_faces_x,
			   const int nb_faces_y,
			   const int nb_cells,
			   const int nb_nodes,
			   const RealType dx,
			   const RealType dy,
			   const RealType dt,
			   const RealType halo_width,
			   const RealType gamma_1,
			   const RealType gamma_2,
			   const RealType pi_1,
			   const RealType pi_2,
			   const RealType* RESTRICT predicted_v,
			   const RealType* RESTRICT e_lag,
			   const RealType* RESTRICT e_1_lag,
			   const RealType* RESTRICT e_2_lag,
			   const RealType* RESTRICT u_lag,
			   const RealType* RESTRICT v_lag,
			   const RealType* RESTRICT in_rho_1,
			   const RealType* RESTRICT in_rho_2,
			   const RealType* RESTRICT in_cell_mass,
			   const RealType* RESTRICT in_cell_mass_1,
			   const RealType* RESTRICT in_cell_mass_2,
			   const RealType* RESTRICT in_c_1,
			   const RealType* RESTRICT in_c_2,
			   const RealType* RESTRICT in_cell_volumic_fraction,
			   const RealType* RESTRICT cell_volumes,
			   // out
			   RealType* RESTRICT out_u,
			   RealType* RESTRICT out_v,
			   RealType* RESTRICT out_e,
			   RealType* RESTRICT out_cell_mass,
			   RealType* RESTRICT out_cell_mass_1,
			   RealType* RESTRICT out_cell_mass_2,
			   RealType* RESTRICT out_e_1,
			   RealType* RESTRICT out_e_2,
			   RealType* RESTRICT out_rho,
			   RealType* RESTRICT out_rho_1,
			   RealType* RESTRICT out_rho_2,
			   RealType* RESTRICT out_c_1,
			   RealType* RESTRICT out_c_2,
			   RealType* RESTRICT out_cell_volumic_fraction,
			   RealType* RESTRICT directional_lagrangian_volume_y,
			   RealType* RESTRICT directional_lagrangian_density_y,
			   RealType* RESTRICT directional_lagrangian_density_1_y,
			   RealType* RESTRICT directional_lagrangian_density_2_y,
			   RealType* RESTRICT volume_fluxes_y,
			   RealType* RESTRICT volume_fluxes_1_y,
			   RealType* RESTRICT volume_fluxes_2_y,
			   RealType* RESTRICT mass_flux_y,
			   RealType* RESTRICT mass_flux_1_y,
			   RealType* RESTRICT mass_flux_2_y,
			   RealType* RESTRICT energy_flux_y,
			   RealType* RESTRICT energy_flux_1_y,
			   RealType* RESTRICT energy_flux_2_y,
			   RealType* RESTRICT concentration_flux_y,
			   RealType* RESTRICT density_gradient_y,
			   RealType* RESTRICT density_1_gradient_y,
			   RealType* RESTRICT density_2_gradient_y,
			   RealType* RESTRICT energy_gradient_y,
			   RealType* RESTRICT energy_1_gradient_y,
			   RealType* RESTRICT energy_2_gradient_y,
			   RealType* RESTRICT concentration_gradient_y,
			   RealType* RESTRICT gradient_u_y,
			   RealType* RESTRICT gradient_v_y,
			   RealType* RESTRICT interface_normal_x,
			   RealType* RESTRICT interface_normal_y,
			   //timing
			   std::vector<RealType>& time_compute_volume_fluxes_Y,
			   std::vector<RealType>& time_gradient_Y,
			   std::vector<RealType>& time_mass_reconstruct_o2_Y,
			   std::vector<RealType>& time_project_mass_Y,
			   std::vector<RealType>& time_reconstruct_energy_o2_Y,
			   std::vector<RealType>& time_project_energy_Y,
			   std::vector<RealType>& time_gradient_nodal_Y,
			   std::vector<RealType>& time_project_nodal_velocity_Y,
			   std::vector<RealType> time_periodic_boundary) {


  if (TypeOfModel == "Momomaterial") {
    
    AdProjection2dMonomaterialYDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, predicted_v, e_lag, u_lag, v_lag, in_cell_mass, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume_y, directional_lagrangian_density_y, volume_fluxes_y, mass_flux_y, energy_flux_y, density_gradient_y, energy_gradient_y, gradient_u_y, gradient_v_y, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);

  } else if (TypeOfModel == "MultimaterialMix") {
    
    AdProjection2dMultimaterialMixYDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, gamma_1, gamma_2, pi_1, pi_2, predicted_v, e_lag, u_lag, v_lag, in_cell_mass, in_c_1, in_c_2, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_e_1, out_e_2, out_rho, out_rho_1, out_rho_2, out_c_1, out_c_2, directional_lagrangian_volume_y, directional_lagrangian_density_y, volume_fluxes_y, mass_flux_y, energy_flux_y, concentration_flux_y, density_gradient_y, energy_gradient_y, concentration_gradient_y, gradient_u_y, gradient_v_y, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);

  } else if (TypeOfModel == "MultimaterialInterface") {

    AdProjection2dMultimaterialInterfaceYDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, gamma_1, gamma_2, pi_1, pi_2, predicted_v, e_lag, e_1_lag, e_2_lag, u_lag, v_lag, in_rho_1, in_rho_2, in_cell_mass, in_cell_mass_1, in_cell_mass_2, in_cell_volumic_fraction, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_cell_mass_1, out_cell_mass_2, out_e_1, out_e_2, out_rho, out_rho_1, out_rho_2, out_cell_volumic_fraction, directional_lagrangian_volume_y, directional_lagrangian_density_y, directional_lagrangian_density_1_y, directional_lagrangian_density_2_y, volume_fluxes_y, volume_fluxes_1_y, volume_fluxes_2_y, mass_flux_y, mass_flux_1_y, mass_flux_2_y, energy_flux_y, energy_flux_1_y, energy_flux_2_y, density_gradient_y, density_1_gradient_y, density_2_gradient_y, energy_gradient_y, energy_1_gradient_y, energy_2_gradient_y, gradient_u_y, gradient_v_y, interface_normal_x, interface_normal_y, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y,time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);

  } else {
    
    AdProjection2dMonomaterialYDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, predicted_v, e_lag, u_lag, v_lag, in_cell_mass, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume_y, directional_lagrangian_density_y, volume_fluxes_y, mass_flux_y, energy_flux_y, density_gradient_y, energy_gradient_y, gradient_u_y, gradient_v_y, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);

  }
  
}

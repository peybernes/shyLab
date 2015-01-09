// -*- c++ -*- (for emacs users)
#include "direct_projection_corner_fluxes_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"


void DirectProjectionCornerFluxes2dDriver(//in
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
					  const RealType gamma,
					  const RealType gamma_1,
					  const RealType gamma_2,
					  const RealType pi,
					  const RealType pi_1,
					  const RealType pi_2,
					  const RealType* RESTRICT predicted_u,
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
					  const RealType* RESTRICT in_pressure,
					  const RealType* RESTRICT in_p_1,
					  const RealType* RESTRICT in_p_2,
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
					  RealType* RESTRICT directional_lagrangian_volume_y,
					  RealType* RESTRICT directional_lagrangian_density,
					  RealType* RESTRICT directional_lagrangian_density_y,
					  RealType* RESTRICT directional_lagrangian_density_1,
					  RealType* RESTRICT directional_lagrangian_density_1_y,
					  RealType* RESTRICT directional_lagrangian_density_2,
					  RealType* RESTRICT directional_lagrangian_density_2_y,
					  RealType* RESTRICT volume_fluxes_x,
					  RealType* RESTRICT volume_fluxes_y,
					  RealType* RESTRICT volume_fluxes_corner,
					  RealType* RESTRICT volume_fluxes_1_x,
					  RealType* RESTRICT volume_fluxes_1_y,
					  RealType* RESTRICT volume_fluxes_1_corner,
					  RealType* RESTRICT volume_fluxes_2_x,
					  RealType* RESTRICT volume_fluxes_2_y,
					  RealType* RESTRICT volume_fluxes_2_corner,
					  RealType* RESTRICT mass_flux_x,
					  RealType* RESTRICT mass_flux_y,
					  RealType* RESTRICT mass_flux_1_x,
					  RealType* RESTRICT mass_flux_1_y,
					  RealType* RESTRICT mass_flux_2_x,
					  RealType* RESTRICT mass_flux_2_y,
					  RealType* RESTRICT mass_corner_fluxes,
					  RealType* RESTRICT mass_1_corner_fluxes,
					  RealType* RESTRICT mass_2_corner_fluxes,
					  index_t* RESTRICT sign_x_corner_fluxes,
					  index_t* RESTRICT sign_y_corner_fluxes,
					  RealType* RESTRICT energy_flux_x,
					  RealType* RESTRICT energy_flux_y,
					  RealType* RESTRICT energy_flux_corner,
					  RealType* RESTRICT energy_flux_1_x,
					  RealType* RESTRICT energy_flux_1_y,
					  RealType* RESTRICT energy_flux_1_corner,
					  RealType* RESTRICT energy_flux_2_x,
					  RealType* RESTRICT energy_flux_2_y,
					  RealType* RESTRICT energy_flux_2_corner,
					  RealType* RESTRICT concentration_flux_x,
					  RealType* RESTRICT concentration_flux_y,
					  RealType* RESTRICT concentration_flux_corner,
					  RealType* RESTRICT bool_check_fluxes_x,
					  RealType* RESTRICT bool_check_fluxes_y,
					  RealType* RESTRICT density_gradient,
					  RealType* RESTRICT density_gradient_y,
					  RealType* RESTRICT density_gradient_diag,
					  RealType* RESTRICT density_gradient_antidiag,
					  RealType* RESTRICT density_1_gradient,
					  RealType* RESTRICT density_1_gradient_y,
					  RealType* RESTRICT density_1_gradient_diag,
					  RealType* RESTRICT density_1_gradient_antidiag,
					  RealType* RESTRICT density_2_gradient,
					  RealType* RESTRICT density_2_gradient_y,
					  RealType* RESTRICT density_2_gradient_diag,
					  RealType* RESTRICT density_2_gradient_antidiag,
					  RealType* RESTRICT energy_gradient,
					  RealType* RESTRICT energy_gradient_y,
					  RealType* RESTRICT energy_gradient_diag,
					  RealType* RESTRICT energy_gradient_antidiag,
					  RealType* RESTRICT energy_1_gradient,
					  RealType* RESTRICT energy_1_gradient_y,
					  RealType* RESTRICT energy_1_gradient_diag,
					  RealType* RESTRICT energy_1_gradient_antidiag,
					  RealType* RESTRICT energy_2_gradient,
					  RealType* RESTRICT energy_2_gradient_y,
					  RealType* RESTRICT energy_2_gradient_diag,
					  RealType* RESTRICT energy_2_gradient_antidiag,
					  RealType* RESTRICT concentration_gradient,
					  RealType* RESTRICT concentration_gradient_y,
					  RealType* RESTRICT concentration_gradient_diag,
					  RealType* RESTRICT concentration_gradient_antidiag,
					  RealType* RESTRICT reconstructed_density_faces_x,
					  RealType* RESTRICT reconstructed_density_faces_y,
					  RealType* RESTRICT reconstructed_density_1_faces_x,
					  RealType* RESTRICT reconstructed_density_1_faces_y,
					  RealType* RESTRICT reconstructed_density_2_faces_x,
					  RealType* RESTRICT reconstructed_density_2_faces_y,
					  RealType* RESTRICT reconstructed_energy_faces_x,
					  RealType* RESTRICT reconstructed_energy_faces_y,
					  RealType* RESTRICT reconstructed_energy_1_faces_x,
					  RealType* RESTRICT reconstructed_energy_1_faces_y,
					  RealType* RESTRICT reconstructed_energy_2_faces_x,
					  RealType* RESTRICT reconstructed_energy_2_faces_y,
					  RealType* RESTRICT reconstructed_concentration_faces_x,
					  RealType* RESTRICT reconstructed_concentration_faces_y,
					  RealType* RESTRICT gradient_u,
					  RealType* RESTRICT gradient_u_y,
					  RealType* RESTRICT gradient_u_diag,
					  RealType* RESTRICT gradient_u_antidiag,
					  RealType* RESTRICT gradient_v,
					  RealType* RESTRICT gradient_v_y,
					  RealType* RESTRICT gradient_v_diag,
					  RealType* RESTRICT gradient_v_antidiag,
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
    
    DirectProjectionCornerFluxes2dMonomaterialDriver(BoundaryConditions, nx, ny, dx, dy, dt, gamma, pi, halo_width, predicted_u, predicted_v, cell_volumes, e_lag, u_lag, v_lag, in_cell_mass, in_pressure, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, volume_fluxes_y, volume_fluxes_corner, mass_flux_x, mass_flux_y, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes, energy_flux_x, energy_flux_y, energy_flux_corner, density_gradient, density_gradient_y, density_gradient_diag, density_gradient_antidiag, energy_gradient, energy_gradient_y, energy_gradient_diag, energy_gradient_antidiag, reconstructed_density_faces_x, reconstructed_density_faces_y, reconstructed_energy_faces_x, reconstructed_energy_faces_y, gradient_u, gradient_u_y, gradient_u_diag, gradient_u_antidiag, gradient_v, gradient_v_y, gradient_v_diag, gradient_v_antidiag, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);


  } else if (TypeOfModel == "MultimaterialMix") {
    
    DirectProjectionCornerFluxes2dMultimaterialMixDriver(BoundaryConditions, nx, ny, dx, dy, dt, halo_width, gamma_1, gamma_2, pi_1, pi_2, predicted_u, predicted_v, cell_volumes, e_lag, u_lag, v_lag, in_cell_mass, in_pressure, in_c_1, in_c_2, out_u, out_v, out_e, out_cell_mass, out_e_1, out_e_2, out_rho, out_rho_1, out_rho_2, out_c_1, out_c_2, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, volume_fluxes_y, volume_fluxes_corner, mass_flux_x, mass_flux_y, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes, energy_flux_x, energy_flux_y, energy_flux_corner, concentration_flux_x, concentration_flux_y, concentration_flux_corner, density_gradient, density_gradient_y, density_gradient_diag, density_gradient_antidiag, energy_gradient, energy_gradient_y, energy_gradient_diag, energy_gradient_antidiag, concentration_gradient, concentration_gradient_y, concentration_gradient_diag, concentration_gradient_antidiag, reconstructed_density_faces_x, reconstructed_density_faces_y, reconstructed_energy_faces_x, reconstructed_energy_faces_y, reconstructed_concentration_faces_x, reconstructed_concentration_faces_y,gradient_u, gradient_u_y, gradient_u_diag, gradient_u_antidiag, gradient_v, gradient_v_y, gradient_v_diag, gradient_v_antidiag, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);  


  } else if (TypeOfModel == "MultimaterialInterface") {

    DirectProjectionCornerFluxes2dMultimaterialInterfaceDriver(BoundaryConditions, nx, ny, nb_faces_x, nb_faces_y, nb_cells, nb_nodes, dx, dy, dt, halo_width, gamma_1, gamma_2,pi_1, pi_2, predicted_u, predicted_v, e_lag, e_1_lag, e_2_lag, u_lag, v_lag, in_p_1, in_p_2, in_rho_1, in_rho_2, in_cell_mass, in_cell_mass_1, in_cell_mass_2, in_cell_volumic_fraction, cell_volumes, out_u, out_v, out_e, out_cell_mass, out_cell_mass_1, out_cell_mass_2, out_e_1, out_e_2, out_rho, out_rho_1, out_rho_2, out_cell_volumic_fraction, directional_lagrangian_volume, directional_lagrangian_density, directional_lagrangian_density_1, directional_lagrangian_density_2,volume_fluxes_x, volume_fluxes_y, volume_fluxes_1_x, volume_fluxes_1_y, volume_fluxes_1_corner, volume_fluxes_2_x, volume_fluxes_2_y, volume_fluxes_2_corner, mass_flux_x, mass_flux_y, mass_flux_1_x, mass_flux_1_y, mass_flux_2_x, mass_flux_2_y, mass_corner_fluxes, mass_1_corner_fluxes, mass_2_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes, energy_flux_1_x, energy_flux_1_y, energy_flux_1_corner, energy_flux_2_x, energy_flux_2_y, energy_flux_2_corner, bool_check_fluxes_x, bool_check_fluxes_y, density_gradient, density_gradient_y, density_gradient_diag, density_gradient_antidiag, density_1_gradient, density_1_gradient_y, density_1_gradient_diag, density_1_gradient_antidiag, density_2_gradient, density_2_gradient_y, density_2_gradient_diag, density_2_gradient_antidiag, energy_gradient, energy_gradient_y, energy_gradient_diag, energy_gradient_antidiag, energy_1_gradient, energy_1_gradient_y, energy_1_gradient_diag, energy_1_gradient_antidiag, energy_2_gradient, energy_2_gradient_y, energy_2_gradient_diag, energy_2_gradient_antidiag, reconstructed_density_1_faces_x, reconstructed_density_1_faces_y, reconstructed_density_2_faces_x, reconstructed_density_2_faces_y, reconstructed_energy_1_faces_x, reconstructed_energy_1_faces_y, reconstructed_energy_2_faces_x, reconstructed_energy_2_faces_y, gradient_u, gradient_u_y, gradient_u_diag, gradient_u_antidiag, gradient_v, gradient_v_y, gradient_v_diag, gradient_v_antidiag, interface_normal_x, interface_normal_y, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);


  } else {
    
    DirectProjectionCornerFluxes2dMonomaterialDriver(BoundaryConditions, nx, ny, dx, dy, dt, gamma, pi, halo_width, predicted_u, predicted_v, cell_volumes, e_lag, u_lag, v_lag, in_cell_mass, in_pressure, out_u, out_v, out_e, out_cell_mass, out_rho, directional_lagrangian_volume, directional_lagrangian_density, volume_fluxes_x, volume_fluxes_y, volume_fluxes_corner, mass_flux_x, mass_flux_y, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes, energy_flux_x, energy_flux_y, energy_flux_corner, density_gradient, density_gradient_y, density_gradient_diag, density_gradient_antidiag, energy_gradient, energy_gradient_y, energy_gradient_diag, energy_gradient_antidiag, reconstructed_density_faces_x, reconstructed_density_faces_y, reconstructed_energy_faces_x, reconstructed_energy_faces_y, gradient_u, gradient_u_y, gradient_u_diag, gradient_u_antidiag, gradient_v, gradient_v_y, gradient_v_diag, gradient_v_antidiag, time_compute_volume_fluxes_X, time_gradient_X, time_mass_reconstruct_o2_X, time_project_mass_X, time_reconstruct_energy_o2_X, time_project_energy_X, time_gradient_nodal_X, time_project_nodal_velocity_X, time_compute_volume_fluxes_Y, time_gradient_Y, time_mass_reconstruct_o2_Y, time_project_mass_Y, time_reconstruct_energy_o2_Y, time_project_energy_Y, time_gradient_nodal_Y, time_project_nodal_velocity_Y, time_periodic_boundary);
    

  }
  
}

  

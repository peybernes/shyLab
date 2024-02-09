// -*- c++ -*- (for emacs users)
#include "lagrange_fluxes_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>     // std::cout, std::end
#include <stdlib.h>    

#include "kernel_tools.h"
#include "kernels.hpp"



void LagrangeFluxes2dDriver(//in
			      const std::string BoundaryConditions,
			      const std::string TypeOfModel,
			      const int nb_mat,
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
			      RealType* RESTRICT volume_fluxes_x,
			      RealType* RESTRICT volume_fluxes_y,
			      RealType* RESTRICT volume_fluxes_1_x,
			      RealType* RESTRICT volume_fluxes_1_y,
			      RealType* RESTRICT volume_fluxes_2_x,
			      RealType* RESTRICT volume_fluxes_2_y,
			      RealType* RESTRICT mass_flux_x,
			      RealType* RESTRICT mass_flux_y,
			      RealType* RESTRICT mass_flux_1_x,
			      RealType* RESTRICT mass_flux_1_y,
			      RealType* RESTRICT mass_flux_2_x,
			      RealType* RESTRICT mass_flux_2_y,
			      RealType* RESTRICT energy_flux_x,
			      RealType* RESTRICT energy_flux_y,
			      RealType* RESTRICT energy_flux_1_x,
			      RealType* RESTRICT energy_flux_1_y,
			      RealType* RESTRICT energy_flux_2_x,
			      RealType* RESTRICT energy_flux_2_y,
			      RealType* RESTRICT concentration_flux_x,
			      RealType* RESTRICT concentration_flux_y,
			      RealType* RESTRICT bool_check_fluxes_x,
			      RealType* RESTRICT bool_check_fluxes_y,
			      RealType* RESTRICT density_gradient,
			      RealType* RESTRICT density_gradient_y,
			      RealType* RESTRICT density_1_gradient,
			      RealType* RESTRICT density_1_gradient_y,
			      RealType* RESTRICT density_2_gradient,
			      RealType* RESTRICT density_2_gradient_y,
			      RealType* RESTRICT energy_gradient,
			      RealType* RESTRICT energy_gradient_y,
			      RealType* RESTRICT energy_1_gradient,
			      RealType* RESTRICT energy_1_gradient_y,
			      RealType* RESTRICT energy_2_gradient,
			      RealType* RESTRICT energy_2_gradient_y,
			      RealType* RESTRICT concentration_gradient,
			      RealType* RESTRICT concentration_gradient_y,
			      RealType* RESTRICT gradient_u,
			      RealType* RESTRICT gradient_u_y,
			      RealType* RESTRICT gradient_v,
			      RealType* RESTRICT gradient_v_y,
			      RealType* RESTRICT interface_normal_x,
			      RealType* RESTRICT interface_normal_y,
			      RealType** alphak_gradx_left,
			      RealType** alphak_gradx_right,
			      RealType** alphak_grady_bot,
			      RealType** alphak_grady_top,
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

  ComputeGradientAlpha(nx,ny,nb_mat,dx,dy,in_c_1,in_c_2,alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  
  if (BoundaryConditions == "Wall") {
    ComputeGradientAlphaBoundaryWall(nx,ny,nb_mat,dx,dy,in_c_1,in_c_2,alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  }


  for (index_t imat = 0; imat < 2; ++imat) {
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {      
	const int cell_ooo  = (nx * iy) + ix;
	std::cout <<  alphak_grady_bot  [imat][cell_ooo] << std::endl;
	std::cout <<  alphak_grady_top  [imat][cell_ooo] << std::endl;
	std::cout <<  alphak_gradx_left [imat][cell_ooo] << std::endl;
	std::cout <<  alphak_gradx_right[imat][cell_ooo] << std::endl;
      }
    }
  }
  
  std::cout << "Start Lagrange Flux ok !" << std::endl;
  exit(0);
}

  

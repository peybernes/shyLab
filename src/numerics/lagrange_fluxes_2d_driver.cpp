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
			      const int nb_cells,
			      const RealType dx,
			      const RealType dy,
			      const RealType dt,
			      const RealType halo_width,
			      RealType* gamma_mix,
			      RealType* speed_of_sound_mix,			      
			      RealType* gamma_k,
			      RealType* pi_prime_k,
			      RealType**         in_c_k,
			      RealType**         in_rho_k,
			      RealType*          in_rho,
			      const RealType* RESTRICT in_cell_volumic_fraction,
			      const RealType* RESTRICT cell_volumes,
			      // out
			      RealType* in_u_cell,
			      RealType* in_v_cell,
			      RealType* RESTRICT out_u,
			      RealType* RESTRICT out_v,
			      RealType* RESTRICT out_e,
			      RealType* RESTRICT out_e_1,
			      RealType* RESTRICT out_e_2,
			      RealType* RESTRICT out_rho,
			      RealType* RESTRICT out_rho_1,
			      RealType* RESTRICT out_rho_2,
			      RealType**         out_c_k,
			      RealType* RESTRICT out_cell_volumic_fraction,
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
			      RealType*  rho_total_energy,
			      RealType*  rho_e,
			      RealType*  rho_U,
			      RealType*  rho_V,
			      RealType*  pressure,
			      RealType*  beta,
			      RealType*  pi_prime_mix,
			      RealType*  p_plus_pi_prime,
			      RealType*  p_plus_pi_prime_gradx,
			      RealType*  p_plus_pi_prime_grady,	 
			      RealType* rho_e_gradx_left,
			      RealType* rho_e_gradx_right,
			      RealType* rho_e_grady_top,
			      RealType* rho_e_grady_bot,
			      RealType* p_gradx_left,
			      RealType* p_gradx_right,
			      RealType* p_grady_top,
			      RealType* p_grady_bot,
			      RealType* rho_gradx_left,
			      RealType* rho_gradx_right,
			      RealType* rho_grady_top,
			      RealType* rho_grady_bot,
			      RealType** alpha_beta_k,
			      RealType** alpha_beta_k_tmp,
			      RealType** masse_k,
			      RealType** masse_k_tmp,
			      RealType** alphak_gradx_left,
			      RealType** alphak_gradx_right,
			      RealType** alphak_grady_bot,
			      RealType** alphak_grady_top) {


  
  int nb_faces_x = (nx + 1) * ny;
  int nb_faces_y = (ny + 1) * nx;
  // For HLL
  RealType p_xet                [nb_faces_x];
  RealType p_yet                [nb_faces_y];
  RealType u_et                 [nb_faces_x];
  RealType v_et                 [nb_faces_y];
  //
  RealType u2                   [nb_cells];
  RealType v2                   [nb_cells];
  RealType u_gradx_left         [nb_cells];
  RealType u_gradx_right        [nb_cells];
  RealType v_grady_top          [nb_cells];
  RealType v_grady_bot          [nb_cells];
  RealType beta_tmp             [nb_cells];
  RealType rho_U_tmp            [nb_cells];
  RealType rho_V_tmp            [nb_cells];
  RealType rho_total_energy_tmp [nb_cells];
  RealType* alpha_beta_k_tmp_bis [nb_mat];
  RealType* masse_k_tmp_bis      [nb_mat];
  RealType* rhok_gradx           [nb_mat];
  RealType* rhok_grady           [nb_mat];
  for (int k = 0;k < nb_mat; k++) {
    alpha_beta_k_tmp_bis [k]  = new RealType[nb_cells];
    masse_k_tmp_bis      [k]  = new RealType[nb_cells];
    rhok_gradx           [k]  = new RealType[nb_cells];
    rhok_grady           [k]  = new RealType[nb_cells];
  }  
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      const int cell_ooo  = (nx * iy) + ix;
      u_gradx_left         [cell_ooo] = 0.;
      u_gradx_right        [cell_ooo] = 0.;
      v_grady_top          [cell_ooo] = 0.;
      v_grady_bot          [cell_ooo] = 0.;
      u2                   [cell_ooo] = in_u_cell[cell_ooo] * in_u_cell[cell_ooo];
      v2                   [cell_ooo] = in_v_cell[cell_ooo] * in_v_cell[cell_ooo];
      beta                 [cell_ooo] = 1.;
      beta_tmp             [cell_ooo] = 1.;
      rho_U_tmp            [cell_ooo] = rho_U[cell_ooo];
      rho_V_tmp            [cell_ooo] = rho_V[cell_ooo];
      rho_total_energy_tmp [cell_ooo] = rho_total_energy[cell_ooo];
      for (index_t imat = 0; imat < nb_mat; ++imat) {
	alpha_beta_k_tmp    [imat][cell_ooo] = in_c_k          [imat][cell_ooo];
	alpha_beta_k_tmp_bis[imat][cell_ooo] = alpha_beta_k_tmp[imat][cell_ooo];
	masse_k_tmp_bis     [imat][cell_ooo] =  masse_k        [imat][cell_ooo];
      }
    }
  }
  
  ComputeGradientAlpha(nx,ny,nb_mat,dx,dy,in_c_k,alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  
  if (BoundaryConditions == "Wall") {
    ComputeGradientAlphaBoundaryWall(nx,ny,nb_mat,dx,dy,in_c_k,alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  }

  ComputeGradientPplusPiPrime(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,p_plus_pi_prime_gradx,p_plus_pi_prime_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPplusPiPrimeBoundaryWall(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,p_plus_pi_prime_gradx,p_plus_pi_prime_grady);
  }

  ComputeGradientP(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPBoundaryWall(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);
  }

  ComputeGradientRho(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		     alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		     rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		     rhok_gradx,rhok_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientRhoBoundaryWall(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		       alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		       rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		       rhok_gradx,rhok_grady);
  }

  ComputeHLL(nx,ny,dx,dy,
	      rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
	      in_rho,pressure,in_u_cell,in_v_cell,speed_of_sound_mix,
	      p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top,
	      u_gradx_left,u_gradx_right,v_grady_bot,v_grady_top,
	      p_xet,u_et,p_yet,v_et);

  ComputeHLLBoundaryWall(nx,ny,dx,dy,
	      rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
	      in_rho,pressure,in_u_cell,in_v_cell,speed_of_sound_mix,
	      p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top,
	      u_gradx_left,u_gradx_right,v_grady_bot,v_grady_top,
	      p_xet,u_et,p_yet,v_et);

  //  for (index_t imat = 0; imat < 2; ++imat) {
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {      
	const int cell_ooo  = (nx * iy) + ix;
	//std::cout <<   p_xet[cell_ooo] << std::endl;
	//std::cout <<   p_yet[cell_ooo] << std::endl;
	// std::cout <<  rhok_gradx  [imat][cell_ooo] << std::endl;
	// std::cout <<  rho_grady_top  [cell_ooo] << std::endl;
	// std::cout << p_gradx_left[cell_ooo] << " " << p_grady_top[cell_ooo] << std::endl;
	// std::cout <<  alphak_grady_bot  [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_grady_top  [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_gradx_left [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_gradx_right[imat][cell_ooo] << std::endl;
      }
    }
    //  }
  
}

  

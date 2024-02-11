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
			      const RealType* RESTRICT in_cell_volumic_fraction,
			      const RealType* RESTRICT cell_volumes,
			      // out
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
			      RealType** alphak_gradx_left,
			      RealType** alphak_gradx_right,
			      RealType** alphak_grady_bot,
			      RealType** alphak_grady_top) {

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


  for (index_t imat = 0; imat < 2; ++imat) {
    for (index_t iy = 1; iy < ny-1; ++iy) {
      for (index_t ix = 1; ix < nx-1; ++ix) {      
	const int cell_ooo  = (nx * iy) + ix;
	std::cout << p_gradx_left[cell_ooo] << " " << p_grady_top[cell_ooo] << std::endl;
	// std::cout <<  alphak_grady_bot  [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_grady_top  [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_gradx_left [imat][cell_ooo] << std::endl;
	// std::cout <<  alphak_gradx_right[imat][cell_ooo] << std::endl;
      }
    }
  }
  
}

  

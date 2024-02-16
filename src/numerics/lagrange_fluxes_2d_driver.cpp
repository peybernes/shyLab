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
			      RealType dt,
			      const RealType CFL,
			      const RealType halo_width,
			      RealType* gamma_mix,
			      RealType* speed_of_sound_mix,			      
			      RealType* gamma_k,
			      RealType* pi_prime_k,
			      RealType**         in_c_k,
			      RealType**         in_rho_k,
			      RealType*          in_rho,
			      RealType* in_total_energy,
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
			      RealType** alphak_grady_top,
 			      RealType*  rho_total_energy_fluxes_x,
			      RealType*  rho_total_energy_fluxes_y,
			      RealType*  rho_U_fluxes_x,
			      RealType*  rho_U_fluxes_y,
			      RealType*  rho_V_fluxes_x,
			      RealType*  rho_V_fluxes_y,
			      RealType*  beta_fluxes_x,
			      RealType*  beta_fluxes_y,
			      RealType** masse_fluxes_k_x,
			      RealType** masse_fluxes_k_y,
			      RealType** alpha_beta_fluxes_k_x,
			      RealType** alpha_beta_fluxes_k_y
			    ) {


  
  RealType dt_tmp = dt;
  RealType dt_local;
  const RealType epsilon = 10E-14;
  int nb_faces_x = (nx + 1) * ny;
  int nb_faces_y = (ny + 1) * nx;
  // For HLL
  RealType p_xet                [nb_faces_x];
  RealType p_yet                [nb_faces_y];
  RealType u_et                 [nb_faces_x];
  RealType v_et                 [nb_faces_y];
  //
  RealType beta_gradx           [nb_cells];
  RealType beta_grady           [nb_cells];
  RealType u2                   [nb_cells];
  RealType v2                   [nb_cells];
  RealType u2_gradx             [nb_cells];
  RealType v2_gradx             [nb_cells];
  RealType u2_grady             [nb_cells];
  RealType v2_grady             [nb_cells];
  RealType u_gradx_left         [nb_cells];
  RealType u_gradx_right        [nb_cells];
  RealType v_gradx_left         [nb_cells];
  RealType v_gradx_right        [nb_cells];
  RealType u_grady_top          [nb_cells];
  RealType u_grady_bot          [nb_cells];
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
      beta_gradx           [cell_ooo] = 0.;
      beta_grady           [cell_ooo] = 0.;
      u2_gradx             [cell_ooo] = 0.;
      v2_gradx             [cell_ooo] = 0.;
      u2_grady             [cell_ooo] = 0.;
      v2_grady             [cell_ooo] = 0.;
      u_gradx_left         [cell_ooo] = 0.;
      u_gradx_right        [cell_ooo] = 0.;
      v_grady_top          [cell_ooo] = 0.;
      v_grady_bot          [cell_ooo] = 0.;
      v_gradx_left         [cell_ooo] = 0.;
      v_gradx_right        [cell_ooo] = 0.;
      u_grady_top          [cell_ooo] = 0.;
      u_grady_bot          [cell_ooo] = 0.;
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
  
  // RK2: 1st step
  ComputeGradientAlpha(nx,ny,nb_mat,dx,dy,in_c_k,
		  // out		 
		       alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);


  if (BoundaryConditions == "Wall") {
    ComputeGradientAlphaBoundaryWall(nx,ny,nb_mat,dx,dy,in_c_k,
		  // out		 
				     alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  }

  // for (int k = 0;k < nb_mat; k++) {
  //   for (index_t iy = 0; iy < ny; ++iy) {
  //     for (index_t ix = 0; ix < nx; ++ix) {
  // 	const int cell_ooo  = (nx * iy) + ix;
  // 	std::cout << "imat = " << k   << "cell = " << cell_ooo  << " alphak_gradx_left = " << alphak_gradx_left[k][cell_ooo] << std::endl;
  // 	std::cout << "imat = " << k   << "cell = " << cell_ooo  << " alphak_gradx_right = " << alphak_gradx_right[k][cell_ooo] << std::endl;
  // 	std::cout << "imat = " << k   << "cell = " << cell_ooo  << " alphak_grady_bot = " << alphak_grady_bot[k][cell_ooo] << std::endl;
  // 	std::cout << "imat = " << k   << "cell = " << cell_ooo  << " alphak_grady_top = " << alphak_grady_top[k][cell_ooo] << std::endl;
  //     }
  //   }
  // }
  ComputeGradientPplusPiPrime(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,
		  // out		 
			      p_plus_pi_prime_gradx,p_plus_pi_prime_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPplusPiPrimeBoundaryWall(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,
		  // out		 
					    p_plus_pi_prime_gradx,p_plus_pi_prime_grady);
  }

  ComputeGradientP(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		  // out		 
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPBoundaryWall(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,
		   // out		 
		   rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);
  }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "RK1 cell = " << cell_ooo  << " p_gradx_left = " << p_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_e_gradx_left = " << rho_e_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " p_gradx_right = " << p_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_e_gradx_right = " << rho_e_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " p_grady_bot = " << p_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_e_grady_bot = " << rho_e_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " p_grady_top = " << p_grady_top[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_e_grady_top = " << rho_e_grady_top[cell_ooo] << std::endl;
  //   }
  // }
  
  ComputeGradientRho(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		     alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		     // out		 
		     rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		     rhok_gradx,rhok_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientRhoBoundaryWall(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		       alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
  		       // out		 
		       rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		       rhok_gradx,rhok_grady);
  }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_gradx_left = " << rho_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_gradx_right = " << rho_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_grady_bot = " << rho_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK1 cell = " << cell_ooo  << " rho_grady_top = " << rho_grady_top[cell_ooo] << std::endl;
  //   }
  // }
  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "cell RK1 = " << cell_ooo  << " rho_gradx_left = " << rho_gradx_left[cell_ooo] << std::endl;
  //     //std::cout << "cell = " << cell_ooo  << " rho_gradx_right= " << rho_gradx_right[cell_ooo] << std::endl;
  //   }
  // }
 
  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK1 in_rho = " << in_rho[cell_ooo] << std::endl;
  //   }
  // }

      ComputeHLL(nx,ny,dx,dy,
	      rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
	      in_rho,pressure,in_u_cell,in_v_cell,speed_of_sound_mix,
	      p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top,
	      u_gradx_left,u_gradx_right,v_grady_bot,v_grady_top,
	      // out		 
	      p_xet,p_yet,u_et,v_et);
  
  if (BoundaryConditions == "Wall") {
    ComputeHLLBoundaryWall(nx,ny,dx,dy,
			   rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
			   in_rho,pressure,in_u_cell,in_v_cell,speed_of_sound_mix,
			   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top,
			   u_gradx_left,u_gradx_right,v_grady_bot,v_grady_top,
			   // out		 
			   p_xet,p_yet,u_et,v_et);
  }

  // for(int i=0;i<nb_faces_x;i++){
  //   std::cout<< "face = " << i << "p_xet = " << p_xet[i] << std::endl;
  //   std::cout<< "face = " << i << "u_et = " << u_et[i] << std::endl;
  // }
  // for(int i=0;i<nb_faces_y;i++){
  //   std::cout<< "face = " << i << "p_yet = " << p_yet[i] << std::endl;
  //   std::cout<< "face = " << i << "v_et = " << v_et[i] << std::endl;
  // }

  ComputeHLLFluxesZX(nx,ny,dx,dy,nb_mat,
		  in_rho,rho_e,beta,u_et,p_xet,in_rho_k,rhok_gradx,in_c_k,alphak_gradx_left,rho_gradx_left,rho_e_gradx_left,u_gradx_left,v_gradx_left,
		  in_u_cell,in_v_cell,u2,v2,u2_gradx,v2_gradx,alphak_gradx_right,rho_gradx_right,rho_e_gradx_right,u_gradx_right,v_gradx_right,beta_gradx,
		  // out		 
		  masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x);
  

  if (BoundaryConditions == "Wall") {
    ComputeHLLFluxesZXBoundaryWall(nx,ny,dx,dy,nb_mat,
				   in_rho,rho_e,beta,u_et,p_xet,in_rho_k,rhok_gradx,in_c_k,alphak_gradx_left,rho_gradx_left,rho_e_gradx_left,u_gradx_left,v_gradx_left,
				   in_u_cell,in_v_cell,u2,v2,u2_gradx,v2_gradx,alphak_gradx_right,rho_gradx_right,rho_e_gradx_right,u_gradx_right,v_gradx_right,beta_gradx,
				   // out		 
				   masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x);
  }
  
  // //  for (index_t imat = 0; imat < nb_mat; ++imat) {
  //   for (index_t ix = 0; ix < nb_faces_x; ++ix) {      
  //     std::cout << "face = " << ix << " masse = " <<  masse_fluxes_k_x[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  alpha_beta_fluxes_k_x[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  rho_total_energy_fluxes_x[ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_U_fluxes_x [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_V_fluxes_x [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << beta_fluxes_x [ix] << std::endl;
  //   }
  //   //  }


  ComputeHLLFluxesZY(nx,ny,dx,dy,nb_mat,
		  in_rho,rho_e,beta,v_et,p_yet,in_rho_k,rhok_grady,in_c_k,alphak_grady_bot,rho_grady_bot,rho_e_grady_bot,u_grady_bot,v_grady_bot,
		  in_u_cell,in_v_cell,u2,v2,u2_grady,v2_grady,alphak_grady_top,rho_grady_top,rho_e_grady_top,u_grady_top,v_grady_top,beta_grady,
		  // out		 
		  masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y);
  

  if (BoundaryConditions == "Wall") {
    ComputeHLLFluxesZYBoundaryWall(nx,ny,dx,dy,nb_mat,
				   masse_k,in_rho,rho_e,beta,v_et,p_yet,in_rho_k,rhok_grady,in_c_k,alphak_grady_bot,rho_grady_bot,rho_e_grady_bot,u_grady_bot,v_grady_bot,
				   in_u_cell,in_v_cell,u2,v2,u2_grady,v2_grady,alphak_grady_top,rho_grady_top,rho_e_grady_top,u_grady_top,v_grady_top,beta_grady,
				   // out		 
				   masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y);
  }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "masse = " << masse_k_tmp[0][cell_ooo] << std::endl;
  //     std::cout <<  "masse = " << masse_k_tmp[1][cell_ooo] << std::endl;
  //     std::cout <<  "masse = " << masse_k_tmp[2][cell_ooo] << std::endl;
  //   }
  // }
  dt_local = 0.5 * dt;

  ComputeVariablesFromFluxes(nx,ny,dx,dy,dt_local,nb_mat,epsilon,gamma_k,pi_prime_k,
			     rho_total_energy,rho_U,rho_V,beta,masse_k_tmp,alpha_beta_k_tmp,
			     masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x,
			     masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y,
		  // out		 
			     in_rho,beta,in_rho_k,in_c_k,in_u_cell,in_v_cell,rho_total_energy,rho_U,rho_V,masse_k,alpha_beta_k,alpha_beta_k_tmp,masse_k_tmp,in_total_energy,rho_e,pressure);
  

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK1 in_u_cell = " << in_u_cell[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 in_v_cell = " << in_v_cell[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 rho_U = " << rho_U[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 rho_V = " << rho_V[cell_ooo] << std::endl;
  //   }
  // }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK1 in_rho = " << in_rho[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 rho_total_energy = " << rho_total_energy[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 in_total_energy = " << in_total_energy[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 rho_e = " <<  rho_e[cell_ooo] << std::endl;
  //     std::cout <<  "RK1 pressure = " << pressure[cell_ooo] << std::endl;
  //   }
  // }

  // RK2: 2nd step

  dt = TimeStepLFMix(nx, ny, dx, dy, nb_mat, CFL, gamma_mix, pi_prime_mix, speed_of_sound_mix, gamma_k,pi_prime_k, in_rho, pressure, in_c_k, in_u_cell, in_v_cell);  

  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      const int cell_ooo  = (nx * iy) + ix;
      u2                   [cell_ooo] = in_u_cell[cell_ooo] * in_u_cell[cell_ooo];
      v2                   [cell_ooo] = in_v_cell[cell_ooo] * in_v_cell[cell_ooo];
    }
  }

  ComputeGradientAlpha(nx,ny,nb_mat,dx,dy,in_c_k,
		  // out		 
		       alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  
  if (BoundaryConditions == "Wall") {
    ComputeGradientAlphaBoundaryWall(nx,ny,nb_mat,dx,dy,in_c_k,
		  // out		 
				     alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top);
  }

  ComputeGradientPplusPiPrime(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,
		  // out		 
			      p_plus_pi_prime_gradx,p_plus_pi_prime_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPplusPiPrimeBoundaryWall(nx,ny,dx,dy,pressure,pi_prime_mix,p_plus_pi_prime,
		  // out		 
					    p_plus_pi_prime_gradx,p_plus_pi_prime_grady);
  }

  ComputeGradientP(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		  // out		 
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);

  if (BoundaryConditions == "Wall") {
    ComputeGradientPBoundaryWall(nx,ny,dx,dy,nb_mat,pi_prime_mix,pressure,in_c_k,
		   gamma_k,pi_prime_k,
		   alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		   p_plus_pi_prime_gradx,p_plus_pi_prime_grady,
		   rho_e,
		   // out		 
		   rho_e_gradx_left,rho_e_gradx_right,rho_e_grady_bot,rho_e_grady_top,
		   p_gradx_left,p_gradx_right,p_grady_bot,p_grady_top);
  }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "RK2 cell = " << cell_ooo  << " p_gradx_left = " << p_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_e_gradx_left = " << rho_e_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " p_gradx_right = " << p_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_e_gradx_right = " << rho_e_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " p_grady_bot = " << p_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_e_grady_bot = " << rho_e_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " p_grady_top = " << p_grady_top[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_e_grady_top = " << rho_e_grady_top[cell_ooo] << std::endl;
  //   }
  // }

  ComputeGradientRho(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		     alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
		     // out		 
		     rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		     rhok_gradx,rhok_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientRhoBoundaryWall(nx,ny,dx,dy,nb_mat,in_rho_k,in_c_k,
		       alphak_gradx_left,alphak_gradx_right,alphak_grady_bot,alphak_grady_top,
  		       // out		 
		       rho_gradx_left,rho_gradx_right,rho_grady_top,rho_grady_bot,
		       rhok_gradx,rhok_grady);
  }
  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_gradx_left = " << rho_gradx_left[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_gradx_right = " << rho_gradx_right[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_grady_bot = " << rho_grady_bot[cell_ooo] << std::endl;
  //     std::cout << "RK2 cell = " << cell_ooo  << " rho_grady_top = " << rho_grady_top[cell_ooo] << std::endl;
  //   }
  // }
  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout << "cell RK2 = " << cell_ooo  << " rho_gradx_left = " << rho_gradx_left[cell_ooo] << std::endl;
  //     //std::cout << "cell = " << cell_ooo  << " rho_gradx_right= " << rho_gradx_right[cell_ooo] << std::endl;
  //   }
  // }
  ComputeGradientBeta(nx,ny,dx,dy,beta,
		  // out		 
		       beta_gradx,beta_grady);

  if (BoundaryConditions == "Wall") {
    ComputeGradientBetaBoundaryWall(nx,ny,dx,dy,beta,
		  // out		 
		       beta_gradx,beta_grady);
  }

  ComputeHLLFluxesZX(nx,ny,dx,dy,nb_mat,
		  in_rho,rho_e,beta,u_et,p_xet,in_rho_k,rhok_gradx,in_c_k,alphak_gradx_left,rho_gradx_left,rho_e_gradx_left,u_gradx_left,v_gradx_left,
		  in_u_cell,in_v_cell,u2,v2,u2_gradx,v2_gradx,alphak_gradx_right,rho_gradx_right,rho_e_gradx_right,u_gradx_right,v_gradx_right,beta_gradx,
		  // out		 
		  masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x);
  

  if (BoundaryConditions == "Wall") {
    ComputeHLLFluxesZXBoundaryWall(nx,ny,dx,dy,nb_mat,
				   in_rho,rho_e,beta,u_et,p_xet,in_rho_k,rhok_gradx,in_c_k,alphak_gradx_left,rho_gradx_left,rho_e_gradx_left,u_gradx_left,v_gradx_left,
				   in_u_cell,in_v_cell,u2,v2,u2_gradx,v2_gradx,alphak_gradx_right,rho_gradx_right,rho_e_gradx_right,u_gradx_right,v_gradx_right,beta_gradx,
				   // out		 
				   masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x);
  }
  
  // //  for (index_t imat = 0; imat < nb_mat; ++imat) {
  //   for (index_t ix = 0; ix < nb_faces_x; ++ix) {      
  //     std::cout << "face = " << ix << " masse = " <<  masse_fluxes_k_x[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  alpha_beta_fluxes_k_x[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  rho_total_energy_fluxes_x[ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_U_fluxes_x [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_V_fluxes_x [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << beta_fluxes_x [ix] << std::endl;
  //   }
  //   //  }


  ComputeHLLFluxesZY(nx,ny,dx,dy,nb_mat,
		  in_rho,rho_e,beta,v_et,p_yet,in_rho_k,rhok_grady,in_c_k,alphak_grady_bot,rho_grady_bot,rho_e_grady_bot,u_grady_bot,v_grady_bot,
		  in_u_cell,in_v_cell,u2,v2,u2_grady,v2_grady,alphak_grady_top,rho_grady_top,rho_e_grady_top,u_grady_top,v_grady_top,beta_grady,
		  // out		 
		  masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y);
  
  if (BoundaryConditions == "Wall") {
    
    ComputeHLLFluxesZYBoundaryWall(nx,ny,dx,dy,nb_mat,
				   masse_k,in_rho,rho_e,beta,v_et,p_yet,in_rho_k,rhok_grady,in_c_k,alphak_grady_bot,rho_grady_bot,rho_e_grady_bot,u_grady_bot,v_grady_bot,
				   in_u_cell,in_v_cell,u2,v2,u2_grady,v2_grady,alphak_grady_top,rho_grady_top,rho_e_grady_top,u_grady_top,v_grady_top,beta_grady,
				   // out		 
				   masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y);
  }

  dt_local = dt;
  
  ComputeVariablesFromFluxes(nx,ny,dx,dy,dt_local,nb_mat,epsilon,gamma_k,pi_prime_k,
			     rho_total_energy_tmp,rho_U_tmp,rho_V_tmp,beta_tmp,masse_k_tmp_bis,alpha_beta_k_tmp_bis,
			     masse_fluxes_k_x,alpha_beta_fluxes_k_x,rho_U_fluxes_x,rho_V_fluxes_x,beta_fluxes_x,rho_total_energy_fluxes_x,
			     masse_fluxes_k_y,alpha_beta_fluxes_k_y,rho_U_fluxes_y,rho_V_fluxes_y,beta_fluxes_y,rho_total_energy_fluxes_y,
		  // out		 
			     in_rho,beta,in_rho_k,in_c_k,in_u_cell,in_v_cell,rho_total_energy,rho_U,rho_V,masse_k,alpha_beta_k,alpha_beta_k_tmp,masse_k_tmp,in_total_energy,rho_e,pressure);

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK2 in_rho = " << in_rho[cell_ooo] << std::endl;
  //   }
  // }

  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK2 in_u_cell = " << in_u_cell[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 in_v_cell = " << in_v_cell[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 rho_U = " << rho_U[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 rho_V = " << rho_V[cell_ooo] << std::endl;
  //   }
  // }
  // for (index_t iy = 0; iy < ny; ++iy) {
  //   for (index_t ix = 0; ix < nx; ++ix) {
  //     const int cell_ooo  = (nx * iy) + ix;
  //     std::cout <<  "RK2 in_rho = " << in_rho[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 rho_total_energy = " << rho_total_energy[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 in_total_energy = " << in_total_energy[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 rho_e = " <<  rho_e[cell_ooo] << std::endl;
  //     std::cout <<  "RK2 pressure = " << pressure[cell_ooo] << std::endl;
  //   }
  // }

  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      const int cell_ooo  = (nx * iy) + ix;
      //std::cout << " cell = " << cell_ooo <<  " rho = " << in_rho[cell_ooo] << std::endl;
      //std::cout << " cell = "  << cell_ooo <<  " E = " << in_total_energy[cell_ooo] << std::endl;
      //std::cout << " cell = "  << cell_ooo <<  " ux = " << in_u_cell[cell_ooo] << std::endl;
      //std::cout << " cell = "  << cell_ooo <<  " vy = " << in_v_cell[cell_ooo] << std::endl;
      //std::cout << " cell = "  << cell_ooo <<  " pressure = " << pressure[cell_ooo] << std::endl;
      //std::cout << " cell = "  << cell_ooo <<  " rho_e = " << rho_e[cell_ooo] << std::endl;
    }
  }



  // //  for (index_t imat = 0; imat < nb_mat; ++imat) {
  //   for (index_t ix = 0; ix < nb_faces_y; ++ix) {      
  //     std::cout << "face = " << ix << " masse = " <<  masse_fluxes_k_y[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  alpha_beta_fluxes_k_y[0][ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " <<  rho_total_energy_fluxes_y[ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_U_fluxes_y [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << rho_V_fluxes_y [ix] << std::endl;
  //     std::cout << "face = " << ix << " masse = " << beta_fluxes_y [ix] << std::endl;
  //   }
  //   //  }


  //  for (index_t imat = 0; imat < nb_mat; ++imat) {
    for (index_t ix = 0; ix < nb_faces_x; ++ix) {      
      //std::cout << "face = " << ix << " masse = " << rho_U_fluxes_x [ix] << std::endl;
      //std::cout << "face = " << ix << " masse = " << rho_V_fluxes_x [ix] << std::endl;
      //std::cout << "face = " << ix << " masse = " << beta_fluxes_x [ix] << std::endl;
    }
    //  }
	
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

  

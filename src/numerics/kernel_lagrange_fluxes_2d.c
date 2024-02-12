// -*- c++ -*- (for emacs users)
#include "kernel_lagrange_fluxes_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


// Compute time-step
RealType TimeStepLFMix(int nx,
			      int ny,
			      const RealType dx,
			      const RealType dy,
			      const int nb_mat,
			      const RealType CFL,
		              RealType* gamma_mix,
		              RealType* pi_prime_mix,
		              RealType* speed_of_sound_mix,
			      const RealType* gamma_k,
			      const RealType* pi_prime_k,
			      RealType* density,
			      RealType* pressure,
		              RealType** in_c_k,
			      const RealType* RESTRICT in_velocity_x,
			      const RealType* RESTRICT in_velocity_y) {
 
  RealType h = std::min(dx,dy);
  RealType dt = CFL * h;
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      const int cell_ooo  = (nx * iy) + ix;
      
      RealType tmp_gamma=0;
      RealType tmp_pi=0;
      speed_of_sound_mix[cell_ooo]=0;
      gamma_mix[cell_ooo]=0;
      pi_prime_mix[cell_ooo]=0;
      
      for(int k = 0;k < nb_mat; k++){
	tmp_gamma = tmp_gamma + in_c_k[k][cell_ooo]/(gamma_k[k]-1);
	tmp_pi    = tmp_pi + in_c_k[k][cell_ooo] * pi_prime_k[k] * gamma_k[k]/(gamma_k[k] - 1);
      }
      gamma_mix[cell_ooo] = 1. + 1. / tmp_gamma;
      pi_prime_mix[cell_ooo] = tmp_pi / (tmp_gamma * gamma_mix[cell_ooo]);//test2
      //m_pi[i]=tmp_pi/(tmp_gamma);
      speed_of_sound_mix[cell_ooo] = std::sqrt(gamma_mix[cell_ooo] * (pressure[cell_ooo] + pi_prime_mix[cell_ooo]) / density[cell_ooo]);
      //m_C[i]=std::sqrt((std::abs(m_p[i])*m_gamma[i]+m_pi[i])/m_rho[i]);
      //m_C[i]=350;
      dt = std::min(dt, CFL * h / (std::abs(in_velocity_x[cell_ooo]) + std::abs(in_velocity_y[cell_ooo]) + speed_of_sound_mix[cell_ooo]));
       
    }
  }
  return dt;
}


void ComputeGradientAlpha(index_t nx, 
			  index_t ny, 
			  index_t nb_mat, 
			  RealType dx,
			  RealType dy,
			  RealType** in_c_k,
			  RealType** alphak_gradx_left,
			  RealType** alphak_gradx_right,
			  RealType** alphak_grady_bot,
			  RealType** alphak_grady_top
			  ) {

  RealType gradx, grady;
  RealType phi,tmp;
  RealType zpp,zpm,zmm,zmp;
  RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
  RealType grad_tmp;
  RealType tmp_rhoe;
  RealType sum_alphax, sum_alphay;
  
  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_alpha.h"
      
    }
  }
}

void ComputeGradientPplusPiPrime(index_t nx, 
			  index_t ny, 
			  RealType dx,
			  RealType dy,
			  RealType*  pressure,
			  RealType*  pi_prime_mix, 	 
			  RealType*  p_plus_pi_prime,	 
			  RealType*  p_plus_pi_prime_gradx,
			  RealType*  p_plus_pi_prime_grady	 
			  ) {

  RealType gradx, grady;
  RealType phi,tmp;
  RealType zpp,zpm,zmm,zmp;
  RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
  RealType grad_tmp;
  
  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // all cells
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      const int cell_ooo  = (nx * iy) + ix;
       p_plus_pi_prime[cell_ooo] = pressure[cell_ooo] + pi_prime_mix[cell_ooo];
    }
  }
  // inner cells
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
    }
  }


}





void ComputeGradientP(index_t nx, 
			  index_t    ny, 
			  RealType   dx,
			  RealType   dy,
			  int        nb_mat,
		          RealType*  pi_prime_mix,
			  RealType*  pressure,
		          RealType**  in_c_k,
		          RealType*  gamma_k,
			  RealType*  pi_prime_k,
			  RealType**  alphak_gradx_left,
			  RealType**  alphak_gradx_right,
			  RealType**  alphak_grady_bot,
			  RealType**  alphak_grady_top,
			  RealType*  p_plus_pi_prime_gradx,
			  RealType*  p_plus_pi_prime_grady,
			  RealType*  rho_e,
			  RealType*  rho_e_gradx_left,
			  RealType*  rho_e_gradx_right,
			  RealType*  rho_e_grady_bot,
			  RealType*  rho_e_grady_top,
			  RealType*  p_gradx_left,
			  RealType*  p_gradx_right,
			  RealType*  p_grady_bot,
			  RealType*  p_grady_top) {

  RealType gradx, grady;
  RealType phi,tmp;
  RealType zpp,zpm,zmm,zmp;
  RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
  RealType grad_tmp;
  RealType tmp_gamma_r;
  RealType tmp_pi_r;
  RealType tmp_gamma_l;
  RealType tmp_pi_l;
  RealType tmp_gamma_t;
  RealType tmp_pi_t;
  RealType tmp_gamma_b;
  RealType tmp_pi_b;
  
  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p.h"
    }
  }


}


void ComputeGradientRho(index_t nx, 
			index_t    ny, 
			RealType   dx,
			RealType   dy,
			int        nb_mat,
			RealType** in_rho_k,
			RealType**  in_c_k,
			RealType**  alphak_gradx_left,
			RealType**  alphak_gradx_right,
			RealType**  alphak_grady_bot,
			RealType**  alphak_grady_top,
			RealType*  rho_gradx_left,
			RealType*  rho_gradx_right,
			RealType*  rho_grady_top,
			RealType*  rho_grady_bot,
			RealType**  rhok_gradx,
			RealType**  rhok_grady
			) {

  RealType gradx, grady;
  RealType phi,tmp;
  RealType zpp,zpm,zmm,zmp;
  RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
  RealType grad_tmp;
  RealType sum_alphay;
  RealType sum_alphax;
  
  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_rho.h"
    }
  }


}


void ComputeHLL(index_t nx, 
		index_t    ny, 
		RealType   dx,
		RealType   dy,
		RealType*  rho_gradx_left,
		RealType*  rho_gradx_right,
		RealType*  rho_grady_top,
		RealType*  rho_grady_bot,
		RealType*  in_rho,
		RealType*  pressure,
		RealType*  in_u_cell,
		RealType*  in_v_cell,
		RealType*  speed_of_sound_mix,
		RealType*  p_gradx_left,
		RealType*  p_gradx_right,
		RealType*  p_grady_bot,
		RealType*  p_grady_top,
		RealType*  u_gradx_left,
		RealType*  u_gradx_right,
		RealType*  v_grady_bot,
		RealType*  v_grady_top,
		RealType*  p_xet,
		RealType*  p_yet,
		RealType*  u_et,
		RealType*  v_et
			) {

  RealType gradx, grady;
  RealType phi,tmp;
  RealType zpp,zpm,zmm,zmp;
  RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
  RealType grad_tmp;
  RealType sum_alphay;
  RealType sum_alphax;
  
  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
    for (index_t iy = 0; iy < ny - 1; ++iy) {
      for (index_t ix = 0; ix < nx - 1; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const int cell_op1  = CellCellOP1 (cell_ooo, nx);
	const int cell_p1o  = CellCellP1O (cell_ooo, nx);
	
	const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
	const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"
    // if  (face_top == 140) {
    //   std::cout << in_rho [cell_ooo] << " " << rho_grady_top[cell_ooo] << " " <<  p_yet[140] << " " << speed_of_sound_mix[cell_ooo] << " " << v_grady_bot[cell_op1] << std::endl;
    //   std::cout << pressure[cell_op1] << " " << p_grady_bot   [cell_op1] << " " <<  p_yet[140] << " " << speed_of_sound_mix[cell_ooo] << " " << v_grady_bot[cell_op1] << std::endl;
    // }
	//std::cout << "k = " << cell_ooo << "facetop = " << face_top  << " and pyet = " <<  p_yet[face_top] << std::endl;
	// std::cout << "k = " << cell_ooo << "faceright = " << face_right  << " and pxet = " <<  p_xet[face_right] << std::endl;
	// std::cout << "k = " << cell_ooo << "faceright = " << face_right  << " and u_et = " <<  u_et[face_right] << std::endl;

    }
  }


}





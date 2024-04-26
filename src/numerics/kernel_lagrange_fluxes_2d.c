// -*- c++ -*- (for emacs users)
#include "kernel_lagrange_fluxes_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"
#include "offload_acc_omp.h"

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
#pragma omp parallel for collapse(2) reduction(min:dt) //present(in_c_k[:nb_mat][:nx*ny],speed_of_sound_mix[:nx*ny],gamma_mix[:nx*ny],pi_prime_mix[:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],pressure[:nx*ny],in_velocity_x[:nx*ny],in_velocity_y[:nx*ny],density[:nx*ny])
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

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
  //GPU_PARALLEL_LOOP_ALL_LEVELS_COLLAPSE2 
#pragma omp parallel for collapse(2) //present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      RealType gradx, grady;
      RealType phi,tmp;
      RealType zpp,zpm,zmm,zmp;
      RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
      RealType grad_tmp;
      RealType tmp_rhoe;
      RealType sum_alphax, sum_alphay;
  
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

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // all cells
#pragma omp parallel for collapse(2) //present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny])
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      RealType gradx, grady;
      RealType phi,tmp;
      RealType zpp,zpm,zmm,zmp;
      RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
      RealType grad_tmp;
  
      const int cell_ooo  = (nx * iy) + ix;
       p_plus_pi_prime[cell_ooo] = pressure[cell_ooo] + pi_prime_mix[cell_ooo];
    }
  }
  // inner cells
#pragma omp parallel for collapse(2) //present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      RealType gradx, grady;
      RealType phi,tmp;
      RealType zpp,zpm,zmm,zmp;
      RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
      RealType grad_tmp;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_y;
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

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
#pragma omp parallel for collapse(2) //present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
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

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
#pragma omp parallel for collapse(2) //present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      RealType gradx, grady;
      RealType phi,tmp;
      RealType zpp,zpm,zmm,zmp;
      RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
      RealType grad_tmp;
      RealType sum_alphay;
      RealType sum_alphax;
  
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
		index_t    nb_faces_x, 
		index_t    nb_faces_y, 
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

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
#pragma omp parallel for collapse(2) //present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  for (index_t iy = 0; iy < ny - 1; ++iy) {
      for (index_t ix = 0; ix < nx - 1; ++ix) {

	RealType gradx, grady;
	RealType phi,tmp;
	RealType zpp,zpm,zmm,zmp;
	RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
	RealType grad_tmp;
	RealType sum_alphay;
	RealType sum_alphax;
	
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


void ComputeHLLFluxesZX(index_t nx, 
			index_t    ny, 
			index_t    nb_faces_x, 
			RealType   dx,
			RealType   dy,
			int        nb_mat,
			RealType*  in_rho,
			RealType*  rho_e,
			RealType*  beta,
			RealType*  u_et,
			RealType*  p_xet,
			RealType** in_rho_k,
			RealType** rhok_gradx,
			RealType** in_c_k,
			RealType** alphak_gradx_left,
			RealType*  rho_gradx_left,
			RealType*  rho_e_gradx_left,
			RealType*  u_gradx_left,
			RealType*  v_gradx_left,
			RealType*  in_u_cell,
			RealType*  in_v_cell,
			RealType*  u2,
			RealType*  v2,
			RealType*  u2_gradx,
			RealType*  v2_gradx,
			RealType** alphak_gradx_right,
			RealType*  rho_gradx_right,
			RealType*  rho_e_gradx_right,
			RealType*  u_gradx_right,
			RealType*  v_gradx_right,
			RealType*  beta_gradx,
			RealType** masse_fluxes_k_x,
			RealType** alpha_beta_fluxes_k_x,
			RealType* rho_U_fluxes_x,
			RealType* rho_V_fluxes_x,
			RealType* beta_fluxes_x,
			RealType* rho_total_energy_fluxes_x
			) {

#pragma omp parallel for collapse(2) //present(in_rho[:nx*ny],rho_e[:nx*ny],beta[:nx*ny],u_et[:nb_faces_x],p_xet[:nb_faces_x],in_rho_k[:nb_mat][:nx*ny],rhok_gradx[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_e_gradx_left[:nx*ny],u_gradx_left[:nx*ny],v_gradx_left[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],u2[:nx*ny],v2[:nx*ny],u2_gradx[:nx*ny],v2_gradx[:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],rho_gradx_right[:nx*ny],rho_e_gradx_right[:nx*ny],u_gradx_right[:nx*ny],v_gradx_right[:nx*ny],beta_gradx[:nx*ny],masse_fluxes_k_x[:nb_mat][:nb_faces_x],alpha_beta_fluxes_k_x[:nb_mat][:nb_faces_x],rho_U_fluxes_x[:nb_faces_x],rho_V_fluxes_x[:nb_faces_x],beta_fluxes_x[:nb_faces_x],rho_total_energy_fluxes_x[:nb_faces_x])
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx - 1; ++ix) {
      
      const index_t cell_ooo = (nx * iy) + ix;
      
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
      
      if (u_et[face_right] >= 0.) {
#include "lf_compute_flux_z_moins_x.h"	
      }
      else if (u_et[face_right] < 0.) {
#include "lf_compute_flux_z_plus_x.h"	
      }
      //printf("face_right = %d masse = %f \n",face_right,masse_fluxes_k_x[0][face_right]);
      //printf("face_right = %d  rho = %f u_et = %f\n",face_right,in_rho_k[0][cell_p1o],u_et[face_right]);
      
    }
  }

}


void ComputeHLLFluxesZY(index_t nx, 
			index_t    ny, 
			index_t    nb_faces_y, 
			RealType   dx,
			RealType   dy,
			int        nb_mat,
			RealType*  in_rho,
			RealType*  rho_e,
			RealType*  beta,
			RealType*  v_et,
			RealType*  p_yet,
			RealType** in_rho_k,
			RealType** rhok_grady,
			RealType** in_c_k,
			RealType** alphak_grady_bot,
			RealType*  rho_grady_bot,
			RealType*  rho_e_grady_bot,
			RealType*  u_grady_bot,
			RealType*  v_grady_bot,
			RealType*  in_u_cell,
			RealType*  in_v_cell,
			RealType*  u2,
			RealType*  v2,
			RealType*  u2_grady,
			RealType*  v2_grady,
			RealType** alphak_grady_top,
			RealType*  rho_grady_top,
			RealType*  rho_e_grady_top,
			RealType*  u_grady_top,
			RealType*  v_grady_top,
			RealType*  beta_grady,
			RealType** masse_fluxes_k_y,
			RealType** alpha_beta_fluxes_k_y,
			RealType* rho_U_fluxes_y,
			RealType* rho_V_fluxes_y,
			RealType* beta_fluxes_y,
			RealType* rho_total_energy_fluxes_y
			) {

#pragma omp parallel for collapse(2) //present(in_rho[:nx*ny],rho_e[:nx*ny],beta[:nx*ny],v_et[:nb_faces_y],p_yet[:nb_faces_y],in_rho_k[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],rho_grady_bot[:nx*ny],rho_e_grady_bot[:nx*ny],u_grady_bot[:nx*ny],v_grady_bot[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],u2[:nx*ny],v2[:nx*ny],u2_grady[:nx*ny],v2_grady[:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_grady_top[:nx*ny],rho_e_grady_top[:nx*ny],u_grady_top[:nx*ny],v_grady_top[:nx*ny],beta_grady[:nx*ny],masse_fluxes_k_y[:nb_mat][:nb_faces_y],alpha_beta_fluxes_k_y[:nb_mat][:nb_faces_y],rho_U_fluxes_y[:nb_faces_y],rho_V_fluxes_y[:nb_faces_y],beta_fluxes_y[:nb_faces_y],rho_total_energy_fluxes_y[:nb_faces_y])
  for (index_t iy = 0; iy < ny - 1; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      
      const index_t cell_ooo = (nx * iy) + ix;
      
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      
      const index_t face_top = CellFaceOP1(cell_ooo, iy, nx);
      
      if (v_et[face_top] >= 0.) {
#include "lf_compute_flux_z_moins_y.h"	
      }
      else if (v_et[face_top] < 0.) {
#include "lf_compute_flux_z_plus_y.h"	
      }
    }
  }

}


void ComputeVariablesFromFluxes(index_t nx, 
			index_t    ny, 
			index_t    nb_faces_x, 
			index_t    nb_faces_y, 
			RealType   dx,
			RealType   dy,
			RealType   dt_local,
			int        nb_mat,
			RealType   epsilon,	
			RealType*  gamma_k,
			RealType*  pi_prime_k,
			RealType*  in_rho_total_energy,
			RealType*  in_rho_U,
			RealType*  in_rho_V,
			RealType*  in_beta,
			RealType** in_masse_k_tmp,
			RealType** in_alpha_beta_k_tmp,
			RealType** masse_fluxes_k_x,
			RealType** alpha_beta_fluxes_k_x,
			RealType*  rho_U_fluxes_x,
			RealType*  rho_V_fluxes_x,
			RealType*  beta_fluxes_x,
			RealType*  rho_total_energy_fluxes_x,
			RealType** masse_fluxes_k_y,
			RealType** alpha_beta_fluxes_k_y,
			RealType*  rho_U_fluxes_y,
			RealType*  rho_V_fluxes_y,
			RealType*  beta_fluxes_y,
			RealType*  rho_total_energy_fluxes_y,
			RealType*  in_rho,
			RealType*  out_beta,
			RealType** in_rho_k,
			RealType** in_c_k,
			RealType*  in_u_cell,
			RealType*  in_v_cell,
			RealType*  out_rho_total_energy,
			RealType*  out_rho_U,
			RealType*  out_rho_V,
			RealType** masse_k,
			RealType** alpha_beta_k,
			RealType** out_alpha_beta_k_tmp,
			RealType** out_masse_k_tmp,
		        RealType*  in_total_energy,
			RealType*  rho_e,
			RealType*  pressure
			) {


  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;

  // inner cells
#pragma omp parallel for collapse(2) //present(gamma_k[:nb_mat],pi_prime_k[:nb_mat],in_rho_total_energy[:nx*ny],in_rho_U[:nx*ny],in_rho_V[:nx*ny],in_beta[:nx*ny],in_masse_k_tmp[:nb_mat][:nx*ny],in_alpha_beta_k_tmp[:nb_mat][:nx*ny],masse_fluxes_k_x[:nb_mat][:nb_faces_x],alpha_beta_fluxes_k_x[:nb_mat][:nb_faces_x],rho_U_fluxes_x[:nb_faces_x],rho_V_fluxes_x[:nb_faces_x],beta_fluxes_x[:nb_faces_x],rho_total_energy_fluxes_x[:nb_faces_x],masse_fluxes_k_y[:nb_mat][:nb_faces_y],alpha_beta_fluxes_k_y[:nb_mat][:nb_faces_y],rho_U_fluxes_y[:nb_faces_y],rho_V_fluxes_y[:nb_faces_y],beta_fluxes_y[:nb_faces_y],rho_total_energy_fluxes_y[:nb_faces_y],in_rho[:nx*ny],out_beta[:nx*ny],in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],out_rho_total_energy[:nx*ny],out_rho_U[:nx*ny],out_rho_V[:nx*ny],masse_k[:nb_mat][:nx*ny],alpha_beta_k[:nb_mat][:nx*ny],out_alpha_beta_k_tmp[:nb_mat][:nx*ny],out_masse_k_tmp[:nb_mat][:nx*ny],in_total_energy[:nx*ny],rho_e[:nx*ny],pressure[:nx*ny])
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      
      const int cell_ooo  = (nx * iy) + ix;
      
      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
      const index_t face_bot   = CellFaceOM1(cell_ooo, iy, nx);
      const index_t face_left  = CellFaceM1O(cell_ooo, iy, nx);
      
      out_rho_total_energy[cell_ooo] = in_rho_total_energy[cell_ooo] - dt_local * (
				    h_x * rho_total_energy_fluxes_x[face_right] - h_x * rho_total_energy_fluxes_x[face_left]
		                  + h_y * rho_total_energy_fluxes_y[face_top]   - h_y * rho_total_energy_fluxes_y[face_bot]
				  );

      out_rho_U[cell_ooo] = in_rho_U[cell_ooo] - dt_local * (
				    h_x * rho_U_fluxes_x[face_right] - h_x * rho_U_fluxes_x[face_left]
		                  + h_y * rho_U_fluxes_y[face_top]   - h_y * rho_U_fluxes_y[face_bot]
				  );

      out_rho_V[cell_ooo] = in_rho_V[cell_ooo]	- dt_local * (
				    h_x * rho_V_fluxes_x[face_right] - h_x * rho_V_fluxes_x[face_left]
		                  + h_y * rho_V_fluxes_y[face_top]   - h_y * rho_V_fluxes_y[face_bot]
				  );

      out_beta[cell_ooo]  = in_beta[cell_ooo] - dt_local * (
			            h_x * beta_fluxes_x[face_right] - h_x * beta_fluxes_x[face_left]
		                  + h_y * beta_fluxes_y[face_top]   - h_y * beta_fluxes_y[face_bot]
				  );
      
      in_rho[cell_ooo] = 0;
      
      for(int k = 0; k < nb_mat; k++){

	out_masse_k_tmp     [k][cell_ooo] = in_masse_k_tmp[k][cell_ooo] - dt_local * (
				    h_x * masse_fluxes_k_x[k][face_right]       - h_x * masse_fluxes_k_x[k][face_left]
		                  + h_y * masse_fluxes_k_y[k][face_top]         - h_y * masse_fluxes_k_y[k][face_bot]
				  );

	out_alpha_beta_k_tmp[k][cell_ooo] = in_alpha_beta_k_tmp[k][cell_ooo] - dt_local * (
				    h_x * alpha_beta_fluxes_k_x[k][face_right] - h_x * alpha_beta_fluxes_k_x[k][face_left]
		                  + h_y * alpha_beta_fluxes_k_y[k][face_top]   - h_y * alpha_beta_fluxes_k_y[k][face_bot]
				  );

	
	in_rho    [cell_ooo] = in_rho[cell_ooo] + out_masse_k_tmp[k][cell_ooo];
	masse_k[k][cell_ooo] = out_masse_k_tmp[k][cell_ooo];
      }

      for(int k = 0; k < nb_mat; k++){
	
	alpha_beta_k[k][cell_ooo] = out_alpha_beta_k_tmp[k][cell_ooo];
	in_c_k      [k][cell_ooo] = out_alpha_beta_k_tmp[k][cell_ooo] / out_beta[cell_ooo];
	
	if(in_c_k[k][cell_ooo] > epsilon) {
	  in_rho_k[k][cell_ooo] = masse_k[k][cell_ooo] / in_c_k[k][cell_ooo];
	}
	else{
	  in_rho_k[k][cell_ooo] = 0.;
	}
      }
      
      in_total_energy[cell_ooo] = out_rho_total_energy[cell_ooo] / in_rho[cell_ooo];
      in_u_cell[cell_ooo]    = out_rho_U[cell_ooo] / in_rho[cell_ooo];
      in_v_cell[cell_ooo]    = out_rho_V[cell_ooo] / in_rho[cell_ooo];
      
    }
  }

  // inner cells
#pragma omp parallel for collapse(2) //present(gamma_k[:nb_mat],pi_prime_k[:nb_mat],in_rho_total_energy[:nx*ny],in_rho_U[:nx*ny],in_rho_V[:nx*ny],in_beta[:nx*ny],in_masse_k_tmp[:nb_mat][:nx*ny],in_alpha_beta_k_tmp[:nb_mat][:nx*ny],masse_fluxes_k_x[:nb_mat][:nb_faces_x],alpha_beta_fluxes_k_x[:nb_mat][:nb_faces_x],rho_U_fluxes_x[:nb_faces_x],rho_V_fluxes_x[:nb_faces_x],beta_fluxes_x[:nb_faces_x],rho_total_energy_fluxes_x[:nb_faces_x],masse_fluxes_k_y[:nb_mat][:nb_faces_y],alpha_beta_fluxes_k_y[:nb_mat][:nb_faces_y],rho_U_fluxes_y[:nb_faces_y],rho_V_fluxes_y[:nb_faces_y],beta_fluxes_y[:nb_faces_y],rho_total_energy_fluxes_y[:nb_faces_y],in_rho[:nx*ny],out_beta[:nx*ny],in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],out_rho_total_energy[:nx*ny],out_rho_U[:nx*ny],out_rho_V[:nx*ny],masse_k[:nb_mat][:nx*ny],alpha_beta_k[:nb_mat][:nx*ny],out_alpha_beta_k_tmp[:nb_mat][:nx*ny],out_masse_k_tmp[:nb_mat][:nx*ny],in_total_energy[:nx*ny],rho_e[:nx*ny],pressure[:nx*ny])
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      
      RealType gamma_tmp = 0;

      const int cell_ooo  = (nx * iy) + ix;
      
      RealType tmp = out_rho_total_energy[cell_ooo]
	    - 0.5 * in_rho[cell_ooo] * in_u_cell[cell_ooo] * in_u_cell[cell_ooo]
	    - 0.5 * in_rho[cell_ooo] * in_v_cell[cell_ooo] * in_v_cell[cell_ooo];

      rho_e[cell_ooo] = tmp;

      gamma_tmp = 0;

      for(int k = 0; k < nb_mat; k++){
	tmp = tmp - in_c_k[k][cell_ooo] * gamma_k[k] * pi_prime_k[k]
	      / (gamma_k[k] - 1);
	gamma_tmp = gamma_tmp + in_c_k[k][cell_ooo] / (gamma_k[k] - 1);
	// std::cout<<m_M->getyk(j*m_NbCell+i)<<std::endl;
      }
      
      pressure[cell_ooo] = tmp / gamma_tmp;

      
    }
  }

  // std::cout << "c_3 257 = " << in_c_k[0][257] <<  " c_3 258 = "<< in_c_k[0][258] << std::endl; 
  // std::cout << "pressure 257 = " << pressure[257] <<  " pressure 258 = "<< pressure[258] << std::endl; 

}




void ComputeGradientBeta(index_t nx, 
			 index_t ny, 
			 RealType dx,
			 RealType dy,
			 RealType*  beta, 	 
			 RealType*  beta_gradx,
			 RealType*  beta_grady	 
			 ) {

  const RealType twelth  = 1./12.;
  const RealType third   = 1./3. ;
  const RealType h_x     = 1./dx ;
  const RealType h_y     = 1./dy ;
  
  // inner cells
#pragma omp parallel for collapse(2) //present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
      
      RealType gradx, grady;
      RealType phi,tmp;
      RealType zpp,zpm,zmm,zmp;
      RealType zpp_tmp,zpm_tmp,zmp_tmp,zmm_tmp;
      RealType grad_tmp;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_y;
    }
  }


}

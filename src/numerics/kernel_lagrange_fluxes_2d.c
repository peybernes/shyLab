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




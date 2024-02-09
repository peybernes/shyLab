// -*- c++ -*- (for emacs users)
#include "kernel_lagrange_fluxes_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


// Compute time-step
RealType TimeStepLFMix(int nx,
			      int ny,
			      const RealType dx,
			      const RealType dy,
			      const RealType CFL,
			      const RealType gamma_1,
			      const RealType gamma_2,
			      const RealType pi_1,
			      const RealType pi_2,
			      const RealType* RESTRICT density_1,
			      const RealType* RESTRICT density_2,
			      const RealType* RESTRICT pressure_1,
			      const RealType* RESTRICT pressure_2,
			      const RealType* RESTRICT in_c_1,
			      const RealType* RESTRICT in_c_2,
			      const RealType* RESTRICT in_velocity_x,
			      const RealType* RESTRICT in_velocity_y) {
 
  RealType max_velocity = 0.0;
  RealType speed_of_sound_1 = 0.0;
  RealType speed_of_sound_2 = 0.0;
  
  //    #pragma omp for reduction(max:max_velocity) nowait
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      const RealType   p_1_ooo = pressure_1[iy * nx + ix];
      const RealType   p_2_ooo = pressure_2[iy * nx + ix];
      const RealType rho_1_ooo =  density_1[iy * nx + ix];
      const RealType rho_2_ooo =  density_2[iy * nx + ix];
      if (in_c_1[iy * nx + ix] > 0.01) { 
	speed_of_sound_1 = SpeedOfSound(gamma_1, rho_1_ooo, p_1_ooo, pi_1);
      }
      else { 
	speed_of_sound_2 = SpeedOfSound(gamma_2, rho_2_ooo, p_2_ooo, pi_2);
      }
      max_velocity = std::max(  max_velocity , fabs(in_velocity_x[iy * nx + ix]) + fabs(in_velocity_y[iy * nx + ix])  + std::max(speed_of_sound_1,speed_of_sound_2));
    }
  }

  RealType dt = CFL * std::min(dx,dy)/max_velocity;
  
  return dt;
}


void ComputeGradientAlpha(index_t nx, 
			  index_t ny, 
			  index_t nb_mat, 
			  RealType dx,
			  RealType dy,
			  const RealType* RESTRICT in_c_1,
			  const RealType* RESTRICT in_c_2,
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




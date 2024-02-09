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
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {
      
      int cell_ooo  = (nx * iy) + ix;
      int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      int cell_om1  = CellCellOM1 (cell_ooo, nx);
      int cell_op1  = CellCellOP1 (cell_ooo, nx);
      int cell_m1o  = CellCellM1O (cell_ooo, nx);
      int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
      sum_alphax = 0;
      sum_alphay = 0;

      for (int k = 0;k<nb_mat;k++) {

	RealType alpha_ooo;  
	RealType alpha_m1m1;
	RealType alpha_p1p1;
	RealType alpha_m1p1;
	RealType alpha_p1m1;
	RealType alpha_om1;
	RealType alpha_op1;
	RealType alpha_m1o;
	RealType alpha_p1o;
	
	if (k==0) {
	  alpha_ooo  = in_c_1[cell_ooo] ;
	  alpha_m1m1 = in_c_1[cell_m1m1];
	  alpha_p1p1 = in_c_1[cell_p1p1];
	  alpha_m1p1 = in_c_1[cell_m1p1];
	  alpha_p1m1 = in_c_1[cell_p1m1];
	  alpha_om1  = in_c_1[cell_om1] ;
	  alpha_op1  = in_c_1[cell_op1] ;
	  alpha_m1o  = in_c_1[cell_m1o] ;
	  alpha_p1o  = in_c_1[cell_p1o] ;
	}
	else {
	  alpha_ooo  = in_c_2[cell_ooo] ;
	  alpha_m1m1 = in_c_2[cell_m1m1];
	  alpha_p1p1 = in_c_2[cell_p1p1];
	  alpha_m1p1 = in_c_2[cell_m1p1];
	  alpha_p1m1 = in_c_2[cell_p1m1];
	  alpha_om1  = in_c_2[cell_om1] ;
	  alpha_op1  = in_c_2[cell_op1] ;
	  alpha_m1o  = in_c_2[cell_m1o] ;
	  alpha_p1o  = in_c_2[cell_p1o] ;
	}
	
	zpp=2.;
	zpm=2.;
	zmp=2.;
	zmm=2.;
	
	gradx = twelth * (alpha_p1p1 - alpha_m1p1 ) +
	  third  * (alpha_p1o  - alpha_m1o  ) +
	  twelth * (alpha_p1m1 - alpha_m1m1 );
	
	grady = twelth * (alpha_p1p1 - alpha_p1m1 ) +
	  third  * (alpha_op1  - alpha_om1  ) +
	  twelth * (alpha_m1p1 - alpha_m1m1 );
	
	grad_tmp = 0.5 * (grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zpp     = std::max(alpha_ooo,alpha_p1o);
	    zpp_tmp = std::max(alpha_op1,alpha_p1p1);
	    zpp     = std::max(zpp,zpp_tmp);
	    zpp     = (zpp - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zpp     = std::min(alpha_ooo,alpha_p1o);
	    zpp_tmp = std::min(alpha_op1,alpha_p1p1);
	    zpp     = std::min(zpp,zpp_tmp);
	    zpp     = (zpp - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = 0.5 * (-grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zpm     = std::max(alpha_ooo,alpha_p1o);
	    zpm_tmp = std::max(alpha_om1,alpha_p1m1);
	    zpm     = std::max(zpm,zpm_tmp);
	    zpm     = (zpm - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zpm     = std::min(alpha_ooo,alpha_p1o);
	    zpm_tmp = std::min(alpha_om1,alpha_p1m1);
	    zpm     = std::min(zpm,zpm_tmp);
	    zpm     = (zpm - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = 0.5 * (grady - gradx);
	
	if(grad_tmp > 0)
	  {
	    zmp     = std::max(alpha_ooo,alpha_m1o);
	    zmp_tmp = std::max(alpha_op1,alpha_m1p1);
	    zmp     = std::max(zmp,zmp_tmp);
	    zmp     = (zmp - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zmp     = std::min(alpha_ooo,alpha_m1o);
	    zmp_tmp = std::min(alpha_op1,alpha_m1p1);
	    zmp     = std::min(zmp,zmp_tmp);
	    zmp     = (zmp - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = -0.5 * (grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zmm     = std::max(alpha_ooo,alpha_m1o);
	    zmm_tmp = std::max(alpha_om1,alpha_m1m1);
	    zmm     = std::max(zmm,zmm_tmp);
	    zmm     = (zmm - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zmm     = std::min(alpha_ooo,alpha_m1o);
	    zmm_tmp = std::min(alpha_om1,alpha_m1m1);
	    zmm     = std::min(zmm,zmm_tmp);
	    zmm     = (zmm - alpha_ooo) / (grad_tmp);
	  }
	
	tmp = std::min(zpp,zmm);
	phi = std::min(zmp,zpm);
	phi = std::min(phi,tmp);
	phi = std::min(phi,2.);
	
	sum_alphax = sum_alphax + 0.5 * phi * gradx;
	sum_alphay = sum_alphay + 0.5 * phi * grady;
	
	alphak_grady_top  [k][cell_ooo] = phi * grady * h_y;
	alphak_gradx_right[k][cell_ooo] = phi * gradx * h_x;
      
      }

      for (int k = 0;k<nb_mat;k++) {
      
	RealType alpha_ooo;  

	if (k==0) {
	  alpha_ooo  = in_c_1[cell_ooo] ;
	}
	else {
	  alpha_ooo  = in_c_2[cell_ooo] ;
	}

	alphak_grady_bot  [k][cell_ooo] = -2 * h_y * (alpha_ooo * (sum_alphay)  - 0.5 * dx * alphak_grady_top  [k][cell_ooo]) / (1 - sum_alphay);
      
	alphak_gradx_left [k][cell_ooo] = -2 * h_x * (alpha_ooo * (sum_alphax)  - 0.5 * dx * alphak_gradx_right[k][cell_ooo]) / (1 - sum_alphax);
      
	alphak_grady_top  [k][cell_ooo] =  2 * h_y * (alpha_ooo * (-sum_alphay) + 0.5 * dx * alphak_grady_top  [k][cell_ooo]) / (1 + sum_alphay);
      
	alphak_gradx_right[k][cell_ooo] =  2 * h_x * (alpha_ooo * (-sum_alphax) + 0.5 * dx * alphak_gradx_right[k][cell_ooo]) / (1 + sum_alphax);

      }
    }
  }


}




// -*- c++ -*- (for emacs users)
#include "kernel_lagrange_fluxes_2d_boundary_conditions.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


// Boundary Wall

void ComputeGradientAlphaBoundaryWall(index_t nx, 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_alpha.h"

    }
    
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_alpha.h"

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const int cell_om1  = cell_ooo;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_alpha.h"
      
    }
    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = cell_ooo;
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_alpha.h"
      
    }
  }
  //Corners
    
  //Xmin Ymin
#pragma acc parallel present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  { const index_t ix = 0;
    const index_t iy = 0;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = cell_ooo;
    const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
    const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
    
#include "lf_compute_gradient_alpha.h"
      
  }

  //Xmax Ymin
#pragma acc parallel present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const int cell_p1m1 = cell_ooo;
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      

#include "lf_compute_gradient_alpha.h"
      
  }

  //Xmin Ymax
#pragma acc parallel present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
    const int cell_m1p1 = cell_ooo;
    const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      

#include "lf_compute_gradient_alpha.h"
      
  }

  //Xmax Ymax
#pragma acc parallel present(in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const int cell_p1p1 = cell_ooo;
    const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
    

#include "lf_compute_gradient_alpha.h"
      
  }

}


void ComputeGradientPplusPiPrimeBoundaryWall(index_t nx, 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;

    }
    
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const int cell_om1  = cell_ooo;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
    }
    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = cell_ooo;
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
    }
  }
  //Corners
    
  //Xmin Ymin
#pragma acc parallel present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = 0;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = cell_ooo;
    const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
    const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
    
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymin
#pragma acc parallel present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const int cell_p1m1 = cell_ooo;
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmin Ymax
#pragma acc parallel present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
    const int cell_m1p1 = cell_ooo;
    const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymax
#pragma acc parallel present(p_plus_pi_prime[:nx*ny],pressure[:nx*ny],pi_prime_mix[:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const int cell_p1p1 = cell_ooo;
    const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
    
#include "lf_compute_gradient_p_plus_pi_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

}

void ComputeGradientPBoundaryWall(index_t nx, 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p.h"

    }
    
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_p.h"

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const int cell_om1  = cell_ooo;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p.h"
      
    }
    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = cell_ooo;
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p.h"
      
    }
  }
  //Corners
    
  //Xmin Ymin
#pragma acc parallel present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = 0;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = cell_ooo;
    const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
    const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
    
#include "lf_compute_gradient_p.h"
      
  }

  //Xmax Ymin
#pragma acc parallel present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const int cell_p1m1 = cell_ooo;
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_p.h"
      
  }

  //Xmin Ymax
#pragma acc parallel present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
    const int cell_m1p1 = cell_ooo;
    const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_p.h"
      
  }

  //Xmax Ymax
#pragma acc parallel present(pi_prime_mix[:nx*ny],pressure[:nx*ny],in_c_k[:nb_mat][:nx*ny],gamma_k[:nb_mat],pi_prime_k[:nb_mat],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],p_plus_pi_prime_gradx[:nx*ny],p_plus_pi_prime_grady[:nx*ny],rho_e[:nx*ny],rho_e_gradx_left[:nx*ny],rho_e_gradx_right[:nx*ny],rho_e_grady_bot[:nx*ny],rho_e_grady_top[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const int cell_p1p1 = cell_ooo;
    const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
    
#include "lf_compute_gradient_p.h"
      
  }

}


void ComputeGradientRhoBoundaryWall(index_t nx, 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_rho.h"

    }
    
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_rho.h"

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const int cell_om1  = cell_ooo;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_rho.h"
      
    }
    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = cell_ooo;
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_rho.h"
      
    }
  }
  //Corners
    
  //Xmin Ymin
#pragma acc parallel present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  { const index_t ix = 0;
    const index_t iy = 0;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = cell_ooo;
    const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
    const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
    
#include "lf_compute_gradient_rho.h"
      
  }

  //Xmax Ymin
#pragma acc parallel present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const int cell_p1m1 = cell_ooo;
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_rho.h"
      
  }

  //Xmin Ymax
#pragma acc parallel present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
    const int cell_m1p1 = cell_ooo;
    const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_rho.h"
      
  }

  //Xmax Ymax
#pragma acc parallel present(in_rho_k[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],rhok_gradx[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const int cell_p1p1 = cell_ooo;
    const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
    
#include "lf_compute_gradient_rho.h"
      
  }

}

void ComputeHLLBoundaryWall(index_t nx, 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_op1  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);

      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"
      
    }
  }
  //Corners
    
  //Xmax Ymin
#pragma acc parallel present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"
      
  }

  //Xmin Ymax
#pragma acc parallel present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
      
#include "lf_compute_hll.h"
      
  }

  //Xmax Ymax
#pragma acc parallel present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = cell_ooo;
    const int cell_p1o  = cell_ooo;
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
    
#include "lf_compute_hll.h"
      
  }

  //Xmin and Xmax 
#pragma acc parallel loop present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  for (index_t iy = 0; iy < ny; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1o  = cell_ooo;
      const int cell_oppo = cell_ooo + nx - 1;
      
      const index_t face_left  = CellFaceM1O(cell_ooo, iy, nx);
      const index_t face_oppo  = CellFaceP1O(cell_oppo, iy, nx);
      
      p_xet[face_left]   = ((in_rho[cell_ooo] * pressure[cell_m1o] + in_rho[cell_m1o] * pressure[cell_ooo]) -
			   ( in_rho[cell_ooo] * in_rho[cell_m1o])
			   * std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_m1o])
			   * (in_u_cell[cell_ooo] - in_u_cell[cell_m1o])) / (in_rho[cell_ooo] + in_rho[cell_m1o]);
    
      u_et[face_left]    = -(in_rho[cell_m1o] * in_u_cell[cell_m1o] + in_rho[cell_ooo] * in_u_cell[cell_ooo] - (pressure[cell_ooo] - pressure[cell_m1o])
			   / std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_m1o])) / (in_rho[cell_ooo] + in_rho[cell_m1o]);

      u_et[face_oppo]    = u_et[face_oppo];
      
    }
  }

#pragma acc parallel loop present(rho_gradx_left[:nx*ny],rho_gradx_right[:nx*ny],rho_grady_bot[:nx*ny],rho_grady_top[:nx*ny],in_rho[:nx*ny],pressure[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],speed_of_sound_mix[:nx*ny],p_gradx_left[:nx*ny],p_gradx_right[:nx*ny],p_grady_bot[:nx*ny],p_grady_top[:nx*ny],u_gradx_left[:nx*ny],u_gradx_right[:nx*ny],v_grady_top[:nx*ny],v_grady_bot[:nx*ny],p_xet[:nb_faces_x],p_yet[:nb_faces_y],u_et[:nb_faces_x],v_et[:nb_faces_y])
  for (index_t ix = 0; ix < nx; ++ix) {

    const index_t iy = 0;
    
    const int cell_ooo   = (nx * iy) + ix;
    const int cell_om1   = cell_ooo;
    const int cell_oppo  = cell_ooo + (ny - 1) * nx;
    
    const index_t face_bot    = CellFaceOM1(cell_ooo, iy, nx);
    const index_t face_oppo   = CellFaceOP1(cell_oppo, iy, nx);
      

    p_yet [face_bot]   = ((in_rho[cell_ooo] * pressure[cell_om1] + in_rho[cell_om1] * pressure[cell_ooo]) -
			  (in_rho[cell_ooo] * in_rho  [cell_om1])
			  *std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_om1])
			  *(in_v_cell[cell_ooo] - in_v_cell[cell_om1])) / (in_rho[cell_ooo] + in_rho[cell_om1]);

    v_et [face_bot]    = - (in_rho[cell_om1] * in_v_cell[cell_om1] + in_rho[cell_ooo] * in_v_cell[cell_ooo] - (pressure[cell_ooo] - pressure[cell_om1])
			   / std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_om1])) / (in_rho[cell_ooo] + in_rho[cell_om1]);

    v_et [face_oppo]   = - v_et [face_oppo];

  }
  
}

void ComputeHLLFluxesZXBoundaryWall(index_t nx, 
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

  //Xmin and Xmax 
#pragma acc parallel loop present(in_rho[:nx*ny],rho_e[:nx*ny],beta[:nx*ny],u_et[:nb_faces_x],p_xet[:nb_faces_x],in_rho_k[:nb_mat][:nx*ny],rhok_gradx[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_gradx_left[:nb_mat][:nx*ny],rho_gradx_left[:nx*ny],rho_e_gradx_left[:nx*ny],u_gradx_left[:nx*ny],v_gradx_left[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],u2[:nx*ny],v2[:nx*ny],u2_gradx[:nx*ny],v2_gradx[:nx*ny],alphak_gradx_right[:nb_mat][:nx*ny],rho_gradx_right[:nx*ny],rho_e_gradx_right[:nx*ny],u_gradx_right[:nx*ny],v_gradx_right[:nx*ny],beta_gradx[:nx*ny],masse_fluxes_k_x[:nb_mat][:nb_faces_x],alpha_beta_fluxes_k_x[:nb_mat][:nb_faces_x],rho_U_fluxes_x[:nb_faces_x],rho_V_fluxes_x[:nb_faces_x],beta_fluxes_x[:nb_faces_x],rho_total_energy_fluxes_x[:nb_faces_x])
  for (index_t iy = 0; iy < ny; ++iy) {
    {
      {
	const index_t ix = nx - 1;
	
	const int cell_ooo  = (nx * iy) + ix;
	const int cell_p1o  = cell_ooo;
	
	const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
	
	if (u_et[face_right] >= 0.) {
#include "lf_compute_flux_z_moins_x.h"	
	}
	else if (u_et[face_right] < 0.) {
#include "lf_compute_flux_z_plus_x.h"	
	}
	
      }

      { const index_t ix = 0;
	
	const int cell_ooo   = (nx * iy) + ix;
	const int face_left  = CellFaceM1O(cell_ooo, iy, nx);

	for(int k = 0; k < nb_mat; k++){
	  
	  masse_fluxes_k_x[k][face_left] = u_et[face_left] * (in_c_k[k][cell_ooo]
	                                  - 0.5 * dx * alphak_gradx_left[k][cell_ooo])
	                                  * (rhok_gradx[k][cell_ooo]
					  - 0.5 * dx * rhok_gradx[k][cell_ooo]);
	  
	  alpha_beta_fluxes_k_x[k][face_left] = u_et[face_left] *
	                                       (in_c_k[k][cell_ooo]
                                 	       - 0.5 * dx * alphak_gradx_left[k][cell_ooo])
                                               * (beta[cell_ooo] - 0.5 * dx * beta_gradx[cell_ooo]);

	}
	rho_U_fluxes_x[face_left] = u_et[face_left]
	                          * (in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo])
	                          * (in_u_cell[cell_ooo] - 0.5 * dx *  u_gradx_left[cell_ooo]) + p_xet[face_left];

	rho_V_fluxes_x[face_left] = u_et[face_left] *
	                           (in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo])
	                          * (in_v_cell[cell_ooo] - 0.5 * dx * v_gradx_left[cell_ooo]);//+p_xet[face_left];

        rho_total_energy_fluxes_x[face_left] = u_et[face_left] *
	                                     ((in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo])
					      * 0.5 * (u2[cell_ooo] - 0.5 * dx * u2_gradx[cell_ooo]
					      + v2[cell_ooo] - 0.5 * dx * v2_gradx[cell_ooo])
					      + (rho_e[cell_ooo] - 0.5 * dx * rho_e_gradx_left[cell_ooo]))
                                              + p_xet[face_left] * u_et[face_left];
				     
	beta_fluxes_x[face_left] = u_et[face_left] * (beta[cell_ooo] - 0.5 * dx * beta_gradx[cell_ooo]);
      }				     
	  
      
    }
  }
  
}


void ComputeHLLFluxesZYBoundaryWall(index_t nx, 
			index_t    ny, 
			index_t    nb_faces_y, 
			RealType   dx,
			RealType   dy,
			int        nb_mat,
			RealType**  masse_k,
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

  //Ymin and Ymax 
#pragma acc parallel loop present(in_rho[:nx*ny],rho_e[:nx*ny],beta[:nx*ny],v_et[:nb_faces_y],p_yet[:nb_faces_y],in_rho_k[:nb_mat][:nx*ny],rhok_grady[:nb_mat][:nx*ny],in_c_k[:nb_mat][:nx*ny],alphak_grady_bot[:nb_mat][:nx*ny],rho_grady_bot[:nx*ny],rho_e_grady_bot[:nx*ny],u_grady_bot[:nx*ny],v_grady_bot[:nx*ny],in_u_cell[:nx*ny],in_v_cell[:nx*ny],u2[:nx*ny],v2[:nx*ny],u2_grady[:nx*ny],v2_grady[:nx*ny],alphak_grady_top[:nb_mat][:nx*ny],rho_grady_top[:nx*ny],rho_e_grady_top[:nx*ny],u_grady_top[:nx*ny],v_grady_top[:nx*ny],beta_grady[:nx*ny],masse_fluxes_k_y[:nb_mat][:nb_faces_y],alpha_beta_fluxes_k_y[:nb_mat][:nb_faces_y],rho_U_fluxes_y[:nb_faces_y],rho_V_fluxes_y[:nb_faces_y],beta_fluxes_y[:nb_faces_y],rho_total_energy_fluxes_y[:nb_faces_y])
  for (index_t ix = 0; ix < nx; ++ix) {

    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_op1  = cell_ooo;

      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      
      if (v_et[face_top] >= 0.) {
#include "lf_compute_flux_z_moins_y.h"	
      }
      else if (v_et[face_top] < 0.) {
#include "lf_compute_flux_z_plus_y.h"	
      }        
    }


    { const index_t iy = 0;
      
      const int cell_ooo   = (nx * iy) + ix;
      const int face_bot   = CellFaceOM1(cell_ooo, iy, nx);
      
      for(int k = 0; k < nb_mat; k++){
	// modif compared to phd Mahy mk ? rho k in following line ?
	masse_fluxes_k_y[k][face_bot] = v_et[face_bot] * masse_k [k][cell_ooo];
	  //(in_rho_k[k][cell_ooo] - 0.5 * dy * rhok_grady[k][cell_ooo]);

	alpha_beta_fluxes_k_y[k][face_bot] = v_et[face_bot] * (in_c_k[k][cell_ooo]
					    - 0.5 * dy * alphak_grady_bot[k][cell_ooo])
 	                                    * (beta[cell_ooo] - 0.5 * dy * beta_grady[cell_ooo]);
	
      }
      
      rho_U_fluxes_y[face_bot] = v_et[face_bot] *
	                            (in_rho[cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo] )*
				 (in_u_cell[cell_ooo] - 0.5 * dy * u_grady_bot[cell_ooo]);//+p_yet[f_posb[l1]];

      rho_V_fluxes_y[face_bot] = v_et[face_bot] *
	                        (in_rho[cell_ooo]    - 0.5 * dy * rho_grady_bot[cell_ooo]) *
	                        (in_v_cell[cell_ooo] - 0.5 * dy * v_grady_bot [cell_ooo]) + p_yet[face_bot];

      rho_total_energy_fluxes_y[face_bot] = v_et[face_bot] *
                     	       ((in_rho[cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo]) *
				0.5 * (u2[cell_ooo] - 0.5 * dy * u2_grady[cell_ooo]
			       + v2[cell_ooo] - 0.5 * dy * v2_grady[cell_ooo])
				+ (rho_e[cell_ooo] - 0.5 * dy * rho_e_grady_bot[cell_ooo]))
                         	+ p_yet[face_bot] * v_et[face_bot];

      beta_fluxes_y[face_bot] = v_et[face_bot] * (beta[cell_ooo] -0.5 * dy * beta_grady[cell_ooo]);

    }

    
  }

}

void ComputeGradientBetaBoundaryWall(index_t nx, 
			 index_t ny, 
			 RealType dx,
			 RealType dy,
			 RealType*  beta, 	 
			 RealType*  beta_gradx,
			 RealType*  beta_grady	 
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
  
  //Xmin and Xmax 
#pragma acc parallel loop present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;

    }
    
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;

    }
  }

  //Ymin and Ymax 
#pragma acc parallel loop present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const int cell_om1  = cell_ooo;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
    }
    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const int cell_om1  = CellCellOM1 (cell_ooo, nx);
      const int cell_op1  = cell_ooo;
      const int cell_m1o  = CellCellM1O (cell_ooo, nx);
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
    }
  }
  //Corners
    
  //Xmin Ymin
#pragma acc parallel present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = 0;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = cell_ooo;
    const int cell_p1p1 = CellCellP1P1(cell_ooo, nx);
    const int cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
    
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymin
#pragma acc parallel present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1p1 = CellCellOP1(cell_ooo, nx);
    const int cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const int cell_p1m1 = cell_ooo;
    const int cell_om1  = cell_ooo;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmin Ymax
#pragma acc parallel present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_p1p1 = CellCellP1O(cell_ooo, nx);
    const int cell_m1p1 = cell_ooo;
    const int cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymax
#pragma acc parallel present(beta_gradx[:nx*ny],beta_grady[:nx*ny],beta[:nx*ny])
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const int cell_p1p1 = cell_ooo;
    const int cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const int cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const int cell_om1  = CellCellOM1 (cell_ooo, nx);
    const int cell_op1  = cell_ooo;
    const int cell_m1o  = CellCellM1O (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
    
#include "lf_compute_gradient_beta.h"
      beta_gradx[cell_ooo] = phi * gradx * h_x;
      beta_grady[cell_ooo] = phi * grady * h_x;
      
  }

}














//       PERIODIC ////////////////////////////////////////////

 void ComputeHLLBoundaryPeriodic(index_t nx, 
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
  
  //Xmin and Xmax 
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    { const index_t ix = nx - 1;
      
      const int cell_ooo  = (nx * iy) + ix;
      const int cell_op1  = CellCellOP1 (cell_ooo, nx);
      const int cell_p1o  = cell_ooo;
      
      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"

    }
  }

  //Ymin and Ymax 
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = ny - 1;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_op1  = cell_ooo;
      const int cell_p1o  = CellCellP1O (cell_ooo, nx);

      const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
      const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"
      
    }
  }
  //Corners
    
  //Xmax Ymin
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = CellCellOP1 (cell_ooo, nx);
    const int cell_p1o  = cell_ooo;
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);

#include "lf_compute_hll.h"
      
  }

  //Xmin Ymax
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = cell_ooo;
    const int cell_p1o  = CellCellP1O (cell_ooo, nx);
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
      
#include "lf_compute_hll.h"
      
  }

  //Xmax Ymax
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
      
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_op1  = cell_ooo;
    const int cell_p1o  = cell_ooo;
      
    const index_t face_top   = CellFaceOP1(cell_ooo, iy, nx);
    const index_t face_right = CellFaceP1O(cell_ooo, iy, nx);
    
#include "lf_compute_hll.h"
      
  }

  //Xmin and Xmax 
  for (index_t iy = 0; iy < ny - 1; ++iy) {
    { const index_t ix = 0;

      const int cell_ooo  = (nx * iy) + ix;
      const int cell_m1o  = cell_ooo;
      
      const index_t face_left  = CellFaceM1O(cell_ooo, iy, nx);
      
      p_xet[face_left]   = (((in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left [cell_ooo]) * (pressure[cell_m1o] - 0.5 * dx * p_gradx_left  [cell_m1o])
			    + (in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left [cell_ooo]) * (pressure[cell_m1o] + 0.5 * dx * p_gradx_right [cell_m1o]))-
			   ((  in_rho[cell_m1o] + 0.5 * dx * rho_gradx_right[cell_m1o]) * (in_rho  [cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo]))
			   * std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_m1o])*
			   (in_u_cell[cell_ooo] - 0.5 * dx * u_gradx_left[cell_ooo] - in_u_cell[cell_m1o] - 0.5 * dx * u_gradx_right[cell_m1o]))
	                   /(in_rho[cell_m1o] + 0.5 * dx * rho_gradx_right[cell_m1o] + in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo]);
      
      u_et[face_left]    = ((in_rho[cell_ooo]  - 0.5 * dx * rho_gradx_left [cell_ooo]) * (in_u_cell[cell_ooo] - 0.5 * dx * u_gradx_left [cell_ooo])
			  + (in_rho[cell_m1o]  + 0.5 * dx * rho_gradx_right[cell_m1o]) * (in_u_cell[cell_m1o] + 0.5 * dx * u_gradx_right[cell_m1o])
			  - (pressure[cell_ooo]- 0.5 * dx * p_gradx_left   [cell_ooo]  -  pressure[cell_m1o]  - 0.5 * dx * p_gradx_right[cell_m1o])
			  / std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_m1o])) /
	                  ( in_rho[cell_m1o] + 0.5 * dx * rho_gradx_right[cell_m1o] + in_rho[cell_ooo] - 0.5 * dx * rho_gradx_left[cell_ooo]);
      
    }
  }

  for (index_t ix = 0; ix < nx - 1; ++ix) {

    const index_t iy = 0;
    
    const int cell_ooo  = (nx * iy) + ix;
    const int cell_om1  = cell_ooo;
    
    const index_t face_bot   = CellFaceOM1(cell_ooo, iy, nx);
      
    p_yet [face_bot] = (((in_rho[cell_om1] + 0.5 * dy * rho_grady_top[cell_om1]) * (pressure[cell_ooo] - 0.5 * dy * p_grady_bot  [cell_ooo])
		       + (in_rho[cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo]) * (pressure[cell_om1] + 0.5 * dy * p_grady_top  [cell_om1])) -
			((in_rho[cell_om1] + 0.5 * dy * rho_grady_top[cell_om1]) * (in_rho  [cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo]))
			 * std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_om1]) *
			 (in_v_cell[cell_ooo] - 0.5 * dy * v_grady_bot  [cell_ooo] -  in_v_cell[cell_om1] - 0.5 * dy * v_grady_top  [cell_om1]))
                        / (in_rho   [cell_om1] + 0.5 * dy * rho_grady_top[cell_om1] +  in_rho   [cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo]);
    v_et [face_bot] = ((   in_rho[cell_ooo]   - 0.5 * dy * rho_grady_bot[cell_ooo]) * (in_v_cell[cell_ooo] - 0.5 * dy * v_grady_bot[cell_ooo])
			+ (in_rho[cell_om1]   + 0.5 * dy * rho_grady_top[cell_om1]) * (in_v_cell[cell_om1] + 0.5 * dy * v_grady_top[cell_om1])
			- (pressure[cell_ooo] - 0.5 * dy * p_grady_bot  [cell_ooo] -   pressure [cell_om1] - 0.5 * dy * p_grady_top[cell_om1])
		        / std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_om1]))
                        / (in_rho[cell_om1] + 0.5 * dy * rho_grady_top[cell_om1] + in_rho[cell_ooo] - 0.5 * dy * rho_grady_bot[cell_ooo]);

  }

  }





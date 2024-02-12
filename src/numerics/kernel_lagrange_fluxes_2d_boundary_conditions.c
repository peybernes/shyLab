// -*- c++ -*- (for emacs users)
#include "kernel_lagrange_fluxes_2d_boundary_conditions.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;

    }
  }

  //Ymin and Ymax 
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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
    }
  }
  //Corners
    
  //Xmin Ymin
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
    
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymin
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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmin Ymax
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
      
#include "lf_compute_gradient_p_plus_p_prime.h"
      p_plus_pi_prime_gradx[cell_ooo] = phi * gradx * h_x;
      p_plus_pi_prime_grady[cell_ooo] = phi * grady * h_x;
      
  }

  //Xmax Ymax
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
    
#include "lf_compute_gradient_p_plus_p_prime.h"
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





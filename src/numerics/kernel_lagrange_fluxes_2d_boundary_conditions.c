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




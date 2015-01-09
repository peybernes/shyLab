// -*- c++ -*- (for emacs users)
#include "kernel_multimat_interface_2d_boundary_conditions.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"




  // Periodic conditions

void ReconstructGradientMultimatXPeriodicBoundary(index_t nx,
						  index_t ny,
						  RealType dx,
						  RealType dy,
						  const RealType* RESTRICT in_cell_volumic_fraction,
						  const RealType* RESTRICT volume_fluxes_1,
						  const RealType* RESTRICT volume_fluxes_2,
						  const RealType* RESTRICT lag_variable_1,
						  const RealType* RESTRICT lag_variable_2,
						  RealType* RESTRICT gradient_variable_1,
						  RealType* RESTRICT gradient_variable_2) {

  //x min
  for (index_t iy = 0; iy < ny; ++iy) {
    index_t ix = 0;
    
    const int cell_ooo = iy * nx + ix;
    const int cell_m1o = cell_ooo + nx - 1;
    const int cell_p1o = CellCellP1O(cell_ooo, nx);
      
    const int face_m2o = CellFaceM1O(cell_m1o, iy, nx);
    const int face_m1o = CellFaceM1O(cell_ooo, iy, nx);
    const int face_p1o = CellFaceP1O(cell_ooo, iy, nx);
    const int face_p2o = CellFaceP1O(cell_p1o, iy, nx);

#include "reconstruct_variable_gradient_multimat_o2_X.h"  
   
    gradient_variable_1[cell_ooo] = limited_grad_variable_1;
    gradient_variable_2[cell_ooo] = limited_grad_variable_2;
      
  }
 
  // x max
  for (index_t iy = 0; iy < ny; ++iy) {
    index_t ix = nx - 1;
    
    const int cell_ooo = iy * nx + ix;
    const int cell_m1o = CellCellM1O(cell_ooo, nx);
    const int cell_p1o = iy * nx;
      
    const int face_m2o = CellFaceM1O(cell_m1o, iy, nx);
    const int face_m1o = CellFaceM1O(cell_ooo, iy, nx);
    const int face_p1o = CellFaceP1O(cell_ooo, iy, nx);
    const int face_p2o = CellFaceP1O(cell_p1o, iy, nx);

#include "reconstruct_variable_gradient_multimat_o2_X.h"

    gradient_variable_1[cell_ooo] = limited_grad_variable_1;
    gradient_variable_2[cell_ooo] = limited_grad_variable_2;
      
  }

} // end ReconstructGradientMultimatXPeriodicBoundary


void ReconstructGradientMultimatYPeriodicBoundary(index_t nx,
						  index_t ny,
						  RealType dx,
						  RealType dy,
						  const RealType* RESTRICT in_cell_volumic_fraction,
						  const RealType* RESTRICT volume_fluxes_1,
						  const RealType* RESTRICT volume_fluxes_2,
						  const RealType* RESTRICT lag_variable_1,
						  const RealType* RESTRICT lag_variable_2,
						  RealType* RESTRICT gradient_variable_1,
						  RealType* RESTRICT gradient_variable_2) {
  
  //y min
  index_t iy = 0;
  index_t iy_sym = ny - 1;
  for (index_t ix = 0; ix < nx; ++ix) {

    const int cell_ooo = iy * nx + ix;
    const int cell_om1 = iy_sym * nx + iy;
    const int cell_op1 = CellCellOP1(cell_ooo, nx);
      
    const int face_om2 = CellFaceOM1(cell_om1, iy_sym, nx);
    const int face_om1 = CellFaceOM1(cell_ooo, iy, nx);
    const int face_op1 = CellFaceOP1(cell_ooo, iy, nx);
    const int face_op2 = CellFaceOP1(cell_op1, iy + 1, nx);

#include "reconstruct_variable_gradient_multimat_o2_Y.h"  

    gradient_variable_1[cell_ooo] = limited_grad_variable_1;
    gradient_variable_2[cell_ooo] = limited_grad_variable_2;
  }

  //y max
  iy = ny - 1;
  iy_sym = 0;
  for (index_t ix = 0; ix < nx; ++ix) {
    
    const int cell_ooo = iy * nx + ix;
    const int cell_om1 = CellCellOM1(cell_ooo, nx);
    const int cell_op1 = iy_sym * nx + ix;
    
    const int face_om2 = CellFaceOM1(cell_om1, iy - 1, nx);
    const int face_om1 = CellFaceOM1(cell_ooo, iy, nx);
    const int face_op1 = CellFaceOP1(cell_ooo, iy, nx);
    const int face_op2 = CellFaceOP1(cell_op1, iy_sym, nx);

#include "reconstruct_variable_gradient_multimat_o2_Y.h"  

    gradient_variable_1[cell_ooo] = limited_grad_variable_1;
    gradient_variable_2[cell_ooo] = limited_grad_variable_2;
  }

} // end ReconstructGradientMultimatYPeriodicBoundary


void ReconstructGradientMultimatDiagPeriodicBoundary(index_t nx,
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* RESTRICT predicted_u,
						     const RealType* RESTRICT predicted_v,
						     const RealType* RESTRICT in_cell_volumic_fraction,
						     const RealType* RESTRICT volume_fluxes_1_corner,
						     const RealType* RESTRICT volume_fluxes_2_corner,
						     const RealType* RESTRICT lag_variable_1,
						     const RealType* RESTRICT lag_variable_2,
						     RealType* RESTRICT gradient_variable_1_diag,
						     RealType* RESTRICT gradient_variable_2_diag,
						     RealType* RESTRICT gradient_variable_1_antidiag,
						     RealType* RESTRICT gradient_variable_2_antidiag) {


  //Xmin and Xmax 
  for (index_t iy = 1; iy < ny - 1; ++iy) {

    { const index_t ix = 0;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_sym = iy * nx + nx - 1;
      const index_t cell_m1m1 = CellCellOM1(cell_sym, nx);
      const index_t cell_m1p1 = CellCellOP1(cell_sym, nx);
      const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_multimat_diag.h"

      gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
      gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
      gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
      gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
    }

    { const index_t ix = nx - 1;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_sym = iy * nx;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellOM1(cell_sym, nx);
      const index_t cell_p1p1 = CellCellOP1(cell_sym, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);     

#include "reconstruct_variable_gradient_multimat_diag.h"

      gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
      gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
      gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
      gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
    }
  }

  //Ymin and Ymax 
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;
      const index_t iy_sym = ny - 1;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_sym = iy_sym * nx + ix;
      const index_t cell_m1m1 = CellCellM1O(cell_sym, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellP1O(cell_sym, nx);
      const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy_sym, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy_sym, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_multimat_diag.h"

      gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
      gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
      gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
      gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
    }

    { const index_t iy = ny - 1;
      const index_t iy_sym = 0;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_sym = iy_sym * nx + ix;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1O(cell_sym, nx);
      const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1O(cell_sym, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy_sym, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy_sym, nx);

#include "reconstruct_variable_gradient_multimat_diag.h"

      gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
      gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
      gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
      gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
    }
  }

  //Corners

  //Xmin Ymin
  { const index_t ix = 0;
    const index_t iy = 0;
    const index_t ix_sym = nx - 1;
    const index_t iy_sym = ny - 1;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_sym_x = iy * nx + ix_sym;
    const index_t cell_sym_y = iy_sym * nx + ix;
    const index_t cell_sym_xy = iy_sym * nx + ix_sym;

    const index_t cell_m1m1 = cell_sym_xy;
    const index_t cell_m1p1 = CellCellOP1(cell_sym_x, nx);
    const index_t cell_p1m1 = CellCellP1O(cell_sym_y, nx);
    const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy_sym, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy_sym, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);
 
#include "reconstruct_variable_gradient_multimat_diag.h"

    gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
    gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
    gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
    gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
  }

  //Xmin Ymax
  { const index_t ix = nx - 1;
    const index_t iy = 0;
    const index_t ix_sym = 0;
    const index_t iy_sym = ny - 1;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_sym_x = iy * nx + ix_sym;
    const index_t cell_sym_y = iy_sym * nx + ix;
    const index_t cell_sym_xy = iy_sym * nx + ix_sym;

    const index_t cell_m1m1 = CellCellM1O(cell_sym_y, nx);
    const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const index_t cell_p1m1 = cell_sym_xy;
    const index_t cell_p1p1 = CellCellOP1(cell_sym_x, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy_sym, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy_sym, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);
 
#include "reconstruct_variable_gradient_multimat_diag.h"

    gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
    gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
    gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
    gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
  }

  //Xmax Ymin
  { const index_t ix = 0;
    const index_t iy = ny - 1;
    const index_t ix_sym = nx - 1;
    const index_t iy_sym = 0;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_sym_x = iy * nx + ix_sym;
    const index_t cell_sym_y = iy_sym * nx + ix;
    const index_t cell_sym_xy = iy_sym * nx + ix_sym;

    const index_t cell_m1m1 = CellCellOM1(cell_sym_x, nx);
    const index_t cell_m1p1 = cell_sym_xy;
    const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const index_t cell_p1p1 = CellCellP1O(cell_sym_y, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy_sym, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy_sym, nx);
 
#include "reconstruct_variable_gradient_multimat_diag.h"

    gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
    gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
    gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
    gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
  }

  //Xmax Ymax
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;
    const index_t ix_sym = 0;
    const index_t iy_sym = 0;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_sym_x = iy * nx + ix_sym;
    const index_t cell_sym_y = iy_sym * nx + ix;
    const index_t cell_sym_xy = iy_sym * nx + ix_sym;

    const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const index_t cell_m1p1 = CellCellM1O(cell_sym_y, nx);
    const index_t cell_p1m1 = CellCellOM1(cell_sym_x, nx);
    const index_t cell_p1p1 = cell_sym_xy;

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy_sym, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy_sym, nx);
 
#include "reconstruct_variable_gradient_multimat_diag.h"

    gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
    gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
    gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
    gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
  }

}//end ReconstructGradientMultimatDiagPeriodicBoundary


void ReconstructNormalInterfacePeriodicBoundary(// in
						index_t nx,
						index_t ny,
						RealType dx,
						RealType dy,
						const RealType* in_variable,
						// out
						RealType* interface_normal_x,
						RealType* interface_normal_y) {


  // x min
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {
      
      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellOM1(cell_m1oo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_m1oo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	} 
      
    }
  }

  // x max
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = nx - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_p1oo, nx);
      const int cell_p1p1 = CellCellOP1(cell_p1oo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	} 

    }
  }

  // y min
  for (index_t ix = 1; ix < nx - 1; ++ix) {
    
    index_t iy = 0;
    index_t iy_sym = ny - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oom1, nx);
      const int cell_p1m1 = CellCellP1O(cell_oom1, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	} 

    }
  }

  // y max
  for (index_t ix = 1; ix < nx - 1; ++ix) {
    
    index_t iy = ny - 1;
    index_t iy_sym = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oop1, nx);
      const int cell_m1p1 = CellCellM1O(cell_oop1, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	} 
      
    }
  }

  // corner x min, y min
  { index_t ix = 0;
    index_t iy = 0;
    index_t iy_sym = ny - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = nx * ny - 1;
      const int cell_p1m1 = CellCellP1O(cell_oom1, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_m1oo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // corner x max, y min
  { index_t ix = nx - 1;
    index_t iy = 0;
    index_t iy_sym = ny - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oom1, nx);
      const int cell_p1m1 = nx * (ny - 1);
      const int cell_p1p1 = CellCellOP1(cell_p1oo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    } 
  }

  // corner x min, y max
  { index_t ix = 0;
    index_t iy = ny - 1;
    index_t iy_sym = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellOM1(cell_m1oo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oop1, nx);
      const int cell_m1p1 = nx - 1; 
      
#include "reconstruct_variable_gradient_Youngs.h"
   
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // corner x max, y max
  { index_t ix = nx - 1;
    index_t iy = ny - 1;
    index_t iy_sym = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_p1oo, nx);
      const int cell_p1p1 = 0;
      const int cell_m1p1 = CellCellM1O(cell_oop1, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

} // end ReconstructNormalInterfacePeriodicBoundary



void ComputeDirectionalLagrangianFractionalVolumeFluxesXPeriodicBoundary(// in
									 index_t nx, 
									 index_t ny, 
									 RealType dx,
									 RealType dy,
									 const RealType* RESTRICT volume_fluxes,
									 const RealType* RESTRICT directional_lagrangian_volume,
									 const RealType* RESTRICT cell_volumic_fraction,
									 const RealType* RESTRICT interface_normal_x,
									 const RealType* RESTRICT interface_normal_y,
									 // out
									 RealType* RESTRICT volume_fluxes_1,
									 RealType* RESTRICT volume_fluxes_2) {


  for (index_t iy = 0; iy < ny; ++iy) {

    // x min and x max (symetry)
    
      const int ix = 0;
      const index_t face_ooo = ((nx + 1) * iy) + ix;
      const index_t face_sym = ((nx + 1) * iy) + nx;

      if ( volume_fluxes[face_ooo] == 0.0 ) {

	volume_fluxes_1[face_ooo] = 0.0;
	volume_fluxes_2[face_ooo] = 0.0;
	volume_fluxes_1[face_sym] = 0.0;
	volume_fluxes_2[face_sym] = 0.0;

      } else {

        index_t cell_ooo = FaceXCellM1(face_sym, iy, nx); 

	if ( volume_fluxes[face_ooo] < 0.0 ) {
	  cell_ooo = FaceXCellP1(face_ooo, iy, nx);
	}

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	RealType vol_flux = volume_fluxes[face_ooo];
	  
	// If only one material in cell_ooo
	if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	  volume_fluxes_1[face_ooo] = vol_fraction_1 * volume_fluxes[face_ooo];
	  volume_fluxes_2[face_ooo] = vol_fraction_2 * volume_fluxes[face_ooo];

	} else { 

	  index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	  index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
	    
	  // Algebraic values of left and right displacement
	  RealType delta_x_prev = volume_fluxes[prev_face] / dy;
	  RealType delta_x_next = volume_fluxes[next_face] / dy;
	    
	  // Vertices of the rectangular approximation of lagrangian cell_ooo
	  RealType xmin = (cell_ooo - cell_ooo / nx * nx ) * dx + delta_x_prev;
	  RealType xmax = (cell_ooo - cell_ooo / nx * nx + 1) * dx + delta_x_next;
	  RealType ymin = cell_ooo / nx * dy;
	  RealType ymax = ymin + dy;

	  // Init of variables for placement of the interface on rectangular approximation of lagrangian cell_ooo
	  RealType normal_x = interface_normal_x[cell_ooo];
	  RealType normal_y = interface_normal_y[cell_ooo];
	  // A and B are the two points that define the interface
	  RealType pointA[2];
	  RealType pointB[2];

	  PlaceInterface(xmin, xmax, ymin, ymax, normal_x, normal_y, vol_fraction_1, pointA, pointB);
	  
	  // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	  RealType pointA_d[2];
	  RealType pointB_d[2];
	  RealType x = ix * dx;
	  RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * (x + nx * dx) + (1 - sign(vol_flux)) * (x + delta_x_prev) );
	  RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * (x + nx * dx + delta_x_next) + (1 - sign(vol_flux)) * x );

	  IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin, ymax, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    if (vol_flux > 0.0) {
	      if (normal_x < 0.0) {
		volume_fluxes_1[face_ooo] = vol_flux;
		volume_fluxes_2[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1[face_ooo] = 0.0;
		volume_fluxes_2[face_ooo] = vol_flux;
	      }
	    } else {
	      if (normal_x > 0.0) {
		volume_fluxes_1[face_ooo] = vol_flux;
		volume_fluxes_2[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1[face_ooo] = 0.0;
		volume_fluxes_2[face_ooo] = vol_flux;
	      }
	    }
	  	    
	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

	      // If the whole quantity of one material initially in cell_ooo is in dVol

	      if ( (vol_flux > 0.0) && (pointA_d[0] > nx * dx) && (pointB_d[0] > nx * dx) ) {
		if (normal_x < 0.0) {
		  volume_fluxes_1[face_ooo] = vol_1_lag;
		  volume_fluxes_2[face_ooo] = vol_flux - vol_1_lag;
		} else { 
		  volume_fluxes_1[face_ooo] = vol_flux - vol_2_lag;
		  volume_fluxes_2[face_ooo] = vol_2_lag;
		}
	   
	      } else if ( (vol_flux < 0.0) && (pointA_d[0] < x) && (pointB_d[0] < x) ) {		  
		if (normal_x > 0.0) {
		  volume_fluxes_1[face_ooo] = - vol_1_lag;
		  volume_fluxes_2[face_ooo] = vol_flux + vol_1_lag;
		} else { 
		  volume_fluxes_1[face_ooo] = vol_flux + vol_2_lag;
		  volume_fluxes_2[face_ooo] = - vol_2_lag;
		}
	      
	      } else {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin, ymax, normal_x, normal_y);

		if (vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		}

		volume_fluxes_1[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2[face_ooo] = vol_fraction_2_d * vol_flux;
	    	       
	      }
	      
	    }
	}

	  volume_fluxes_1[face_sym] = volume_fluxes_1[face_ooo];
	  volume_fluxes_2[face_sym] = volume_fluxes_2[face_ooo];
      }
  }
} // end ComputeDirectionalLagrangianFractionalVolumeFluxesXPeriodicBoundary



void ComputeDirectionalLagrangianFractionalVolumeFluxesYPeriodicBoundary(// in
									 index_t nx, 
									 index_t ny, 
									 RealType dx,
									 RealType dy,
									 const RealType* RESTRICT volume_fluxes_y,
									 const RealType* RESTRICT directional_lagrangian_volume,
									 const RealType* RESTRICT cell_volumic_fraction,
									 const RealType* RESTRICT interface_normal_x,
									 const RealType* RESTRICT interface_normal_y,
									 // out
									 RealType* RESTRICT volume_fluxes_1_y,
									 RealType* RESTRICT volume_fluxes_2_y) {
  
  for (index_t ix = 0; ix < nx; ++ix) {

    // y min and y max (symetry)

      const int iy = 0;
      const int iy_sym = ny;
      const index_t face_ooo = nx * iy + ix;
      const index_t face_sym = nx * iy_sym + ix;

      if ( volume_fluxes_y[face_ooo] == 0.0 ) {

	volume_fluxes_1_y[face_ooo] = 0.0;
	volume_fluxes_2_y[face_ooo] = 0.0;
	volume_fluxes_1_y[face_sym] = 0.0;
	volume_fluxes_2_y[face_sym] = 0.0;

      } else {

	index_t cell_ooo = FaceYCellM1(face_sym, iy, nx);

	if ( volume_fluxes_y[face_ooo] < 0.0 ) {
	  cell_ooo = FaceYCellP1(face_ooo, iy, nx);
	}

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	RealType vol_flux = volume_fluxes_y[face_ooo];

	// If only one material in cell_ooo
	if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	  volume_fluxes_1_y[face_ooo] = vol_fraction_1 * volume_fluxes_y[face_ooo];
	  volume_fluxes_2_y[face_ooo] = vol_fraction_2 * volume_fluxes_y[face_ooo];
	  volume_fluxes_1_y[face_sym] = volume_fluxes_1_y[face_ooo];
	  volume_fluxes_2_y[face_sym] = volume_fluxes_2_y[face_ooo];

	} else { 

	  index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	  index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	  // Algebraic values of up and down displacement
	  RealType delta_y_prev = volume_fluxes_y[prev_face] / dx;
	  RealType delta_y_next = volume_fluxes_y[next_face] / dx;

	  // Vertices of the rectangular approximation of lagrangian cell_ooo
	  RealType xmin = (cell_ooo - cell_ooo / nx * nx ) * dx;
	  RealType xmax = xmin + dx;
	  RealType ymin = cell_ooo / nx * dy + delta_y_prev;
	  RealType ymax = (cell_ooo / nx + 1) * dy + delta_y_next;

	  // Init of variables for placement of the interface on rectangular approximation of lagrangian cell_ooo
	  RealType normal_x = interface_normal_x[cell_ooo];
	  RealType normal_y = interface_normal_y[cell_ooo];
	  RealType pointA[2];
	  RealType pointB[2];
	    
	  PlaceInterface(xmin, xmax, ymin, ymax, normal_x, normal_y, vol_fraction_1, pointA, pointB);

	  // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	  RealType pointA_d[2];
	  RealType pointB_d[2];
	  RealType y = iy * dy;
	  RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * (y + ny * dy) + (1 - sign(vol_flux)) * (y + delta_y_prev) );
	  RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * (y + ny * dy + delta_y_next) + (1 - sign(vol_flux)) * y );

	  IntersectLineRectangle(pointA, pointB, xmin, xmax, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    if (vol_flux > 0.0) {
	      if (normal_y < 0.0) {
		volume_fluxes_1_y[face_ooo] = vol_flux;
		volume_fluxes_2_y[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_y[face_ooo] = 0.0;
		volume_fluxes_2_y[face_ooo] = vol_flux;
	      }
	    } else {
	      if (normal_y > 0.0) {
		volume_fluxes_1_y[face_ooo] = vol_flux;
		volume_fluxes_2_y[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_y[face_ooo] = 0.0;
		volume_fluxes_2_y[face_ooo] = vol_flux;
	      }
	    }
	  	    
	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

	      // If the whole quantity of one material initially in cell_ooo is in dVol

	      if ( (vol_flux > 0.0) && (pointA_d[1] > ny * dy) && (pointB_d[1] > ny * dy) ) {
		if (normal_y < 0.0) {
		  volume_fluxes_1_y[face_ooo] = vol_1_lag;
		  volume_fluxes_2_y[face_ooo] = vol_flux - vol_1_lag;
		} else { 
		  volume_fluxes_1_y[face_ooo] = vol_flux - vol_2_lag;
		  volume_fluxes_2_y[face_ooo] = vol_2_lag;
		}
	   
	      } else if ( (vol_flux < 0.0) && (pointA_d[1] < y) && (pointB_d[1] < y) ) {		  
		if (normal_y > 0.0) {
		  volume_fluxes_1_y[face_ooo] = - vol_1_lag;
		  volume_fluxes_2_y[face_ooo] = vol_flux + vol_1_lag;
		} else { 
		  volume_fluxes_1_y[face_ooo] = vol_flux + vol_2_lag;
		  volume_fluxes_2_y[face_ooo] = - vol_2_lag;
		}
	      
	      } else {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin, xmax, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		}

		volume_fluxes_1_y[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2_y[face_ooo] = vol_fraction_2_d * vol_flux;
	    	       
	      }
	      
	    }
	}

	volume_fluxes_1_y[face_sym] = volume_fluxes_1_y[face_ooo];
	volume_fluxes_2_y[face_sym] = volume_fluxes_2_y[face_ooo];
      
      } 
  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesYPeriodicBoundary



void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerPeriodicBoundary(// in
										index_t nx, 
										index_t ny, 
										RealType dx,
										RealType dy,
										RealType dt,
										const RealType* RESTRICT volume_fluxes,
										const RealType* RESTRICT volume_fluxes_y,
										const RealType* RESTRICT directional_lagrangian_volume,
										const RealType* RESTRICT cell_volumic_fraction,
										const RealType* RESTRICT interface_normal_x,
										const RealType* RESTRICT interface_normal_y,
										const RealType* RESTRICT in_vx,
										const RealType* RESTRICT in_vy,
										// out
										RealType* RESTRICT volume_fluxes_1,
										RealType* RESTRICT volume_fluxes_1_y,
										RealType* RESTRICT volume_fluxes_1_corner,
										RealType* RESTRICT volume_fluxes_2,
										RealType* RESTRICT volume_fluxes_2_y,
										RealType* RESTRICT volume_fluxes_2_corner) {
  RealType pointA[nx * ny][2];
  RealType pointB[nx * ny][2];

  // Place Interface on each reconstructed lagragian cell 
  //X edges
  for (index_t iy = 0; iy < ny; ++iy) { 
    // Xmin
    { const index_t ix = 0;
      const index_t cell_ooo = nx * iy + ix;
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo,iy,nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo,iy,nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo,iy,nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo,iy,nx);

      const RealType normal_x = interface_normal_x[cell_ooo];
      const RealType normal_y = interface_normal_y[cell_ooo];
     
      RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
      
      if ((vol_fraction_1 != 0.0) && (vol_fraction_2 != 0.0)) { 

	const RealType xmin_lag = ix * dx + 0.5 * (in_vx[node_m1m1] + in_vx[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (in_vx[node_p1m1] + in_vx[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (in_vy[node_m1m1] + in_vy[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (in_vy[node_p1m1] + in_vy[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA[cell_ooo], pointB[cell_ooo]);

      } else {
	
	pointA[cell_ooo][0] = 1.0e20;
	pointA[cell_ooo][1] = 1.0e20;
	pointB[cell_ooo][0] = 1.0e20;
	pointB[cell_ooo][1] = 1.0e20;
      }	
    }
    // Xmax
    { const index_t ix = nx - 1;
      const index_t cell_ooo = nx * iy + ix;
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo,iy,nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo,iy,nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo,iy,nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo,iy,nx);

      const RealType normal_x = interface_normal_x[cell_ooo];
      const RealType normal_y = interface_normal_y[cell_ooo];
     
      RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      if ((vol_fraction_1 != 0.0) && (vol_fraction_2 != 0.0)) { 

	const RealType xmin_lag = ix * dx + 0.5 * (in_vx[node_m1m1] + in_vx[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (in_vx[node_p1m1] + in_vx[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (in_vy[node_m1m1] + in_vy[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (in_vy[node_p1m1] + in_vy[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA[cell_ooo], pointB[cell_ooo]);

      } else {
	
	pointA[cell_ooo][0] = 1.0e20;
	pointA[cell_ooo][1] = 1.0e20;
	pointB[cell_ooo][0] = 1.0e20;
	pointB[cell_ooo][1] = 1.0e20;
      }
    }
  }

  // Y edges
  for (index_t ix = 0; ix < nx; ++ix) { 
    // Ymin
    { const index_t iy = 0;
      const index_t cell_ooo = nx * iy + ix;
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo,iy,nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo,iy,nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo,iy,nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo,iy,nx);

      const RealType normal_x = interface_normal_x[cell_ooo];
      const RealType normal_y = interface_normal_y[cell_ooo];
     
      RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
      
      if ((vol_fraction_1 != 0.0) && (vol_fraction_2 != 0.0)) { 

	const RealType xmin_lag = ix * dx + 0.5 * (in_vx[node_m1m1] + in_vx[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (in_vx[node_p1m1] + in_vx[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (in_vy[node_m1m1] + in_vy[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (in_vy[node_p1m1] + in_vy[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA[cell_ooo], pointB[cell_ooo]);

      } else {
	
	pointA[cell_ooo][0] = 1.0e20;
	pointA[cell_ooo][1] = 1.0e20;
	pointB[cell_ooo][0] = 1.0e20;
	pointB[cell_ooo][1] = 1.0e20;
      }  
    }
    // Ymax
    { const index_t iy = ny - 1;
      const index_t cell_ooo = nx * iy + ix;
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo,iy,nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo,iy,nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo,iy,nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo,iy,nx);

      const RealType normal_x = interface_normal_x[cell_ooo];
      const RealType normal_y = interface_normal_y[cell_ooo];
     
      RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
      
      if ((vol_fraction_1 != 0.0) && (vol_fraction_2 != 0.0)) { 

	const RealType xmin_lag = ix * dx + 0.5 * (in_vx[node_m1m1] + in_vx[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (in_vx[node_p1m1] + in_vx[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (in_vy[node_m1m1] + in_vy[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (in_vy[node_p1m1] + in_vy[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA[cell_ooo], pointB[cell_ooo]);

      } else {
	
	pointA[cell_ooo][0] = 1.0e20;
	pointA[cell_ooo][1] = 1.0e20;
	pointB[cell_ooo][0] = 1.0e20;
	pointB[cell_ooo][1] = 1.0e20;
      }
    }
  }

  // Faces X
  for (index_t iy = 0; iy < ny; ++iy) {
  
        const index_t ix = 0;
	const index_t face_ooo = ((nx + 1) * iy) + ix;
	const index_t face_sym = ((nx + 1) * iy) + nx;

	if ( volume_fluxes[face_ooo] == 0.0 ) {

	  volume_fluxes_1[face_ooo] = 0.0;
	  volume_fluxes_2[face_ooo] = 0.0;
	  volume_fluxes_1[face_sym] = 0.0;
	  volume_fluxes_2[face_sym] = 0.0;

	} else {

	 index_t cell_ooo = FaceXCellM1(face_sym, iy, nx); 

	  if ( volume_fluxes[face_ooo] < 0.0 ) {
	     cell_ooo = FaceXCellP1(face_ooo, iy, nx);
	  }

	  RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	  RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	  RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	  RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
	 
	  RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_flux = volume_fluxes[face_ooo];

	  // If only one material in cell_ooo
	  if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	    volume_fluxes_1[face_ooo] = vol_fraction_1 * vol_flux;
	    volume_fluxes_2[face_ooo] = vol_fraction_2 * vol_flux;

	  } else { 

	    const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	    const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
	    
	    const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	    const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	    const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	    const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);
	    
	    const RealType normal_x = interface_normal_x[cell_ooo];
	    const RealType normal_y = interface_normal_y[cell_ooo];

	    // Algebraic values of left and right displacement
	    RealType ymin = (cell_ooo / nx) * dy;
	    RealType ymax = (cell_ooo / nx + 1) * dy;
	    RealType ymin_lag_prev = ymin + in_vy[node_m1m1] * dt;
	    RealType ymax_lag_prev = ymax + in_vy[node_m1p1] * dt;
	    RealType ymin_lag_next = ymin + in_vy[node_p1m1] * dt;
	    RealType ymax_lag_next = ymax + in_vy[node_p1p1] * dt;
	    RealType ymin_prev = std::max(ymin_lag_prev,ymin);
	    RealType ymin_next = std::max(ymin_lag_next,ymin);
	    RealType ymax_prev = std::min(ymax_lag_prev,ymax);
	    RealType ymax_next = std::min(ymax_lag_next,ymax); 
	    RealType dy_lag_prev = ymax_prev - ymin_prev;
	    RealType dy_lag_next = ymax_next - ymin_next;
	    
	    RealType delta_x_prev = volume_fluxes[prev_face] / dy_lag_prev;
	    RealType delta_x_next = volume_fluxes[next_face] / dy_lag_next;
	  
	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = ix * dx;
	    RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * (x + nx * dx) + (1 - sign(vol_flux)) * (x + delta_x_prev) );
	    RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * (x + nx * dx + delta_x_next) + (1 - sign(vol_flux)) * x );
	    RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * ymin_next + (1 - sign(vol_flux)) * ymin_prev );
	    RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * ymax_next + (1 - sign(vol_flux)) * ymax_prev );

	    IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[cell_ooo][0] + pointB[cell_ooo][0]) / 2.;
	    RealType yI = (pointA[cell_ooo][1] + pointB[cell_ooo][1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      volume_fluxes_1[face_ooo] = vol_flux;
	      volume_fluxes_2[face_ooo] = 0.0;
	    } else {
	      volume_fluxes_1[face_ooo] = 0.0;
	      volume_fluxes_2[face_ooo] = vol_flux;
	    }	      
	    	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		}

		volume_fluxes_1[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2[face_ooo] = vol_fraction_2_d * vol_flux;  
	    }
	  }

	  volume_fluxes_1[face_sym] = volume_fluxes_1[face_ooo];
	  volume_fluxes_2[face_sym] = volume_fluxes_2[face_ooo]; 
	}

  } //end Faces X


  // Faces Y
  for (index_t ix = 0; ix < nx; ++ix) {

	const index_t iy = 0;
	const index_t iy_sym = ny;
	const index_t face_ooo = nx * iy + ix;
	const index_t face_sym = nx * iy_sym + ix;

	if ( volume_fluxes_y[face_ooo] == 0.0 ) {

	  volume_fluxes_1_y[face_ooo] = 0.0;
	  volume_fluxes_2_y[face_ooo] = 0.0;
	  volume_fluxes_1_y[face_sym] = 0.0;
	  volume_fluxes_2_y[face_sym] = 0.0;

	} else {

	  index_t cell_ooo = FaceYCellM1(face_sym, iy_sym, nx); 

	  if ( volume_fluxes_y[face_ooo] < 0.0 ) {
	    cell_ooo = FaceYCellP1(face_ooo, iy, nx);
	  }

	  const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	  const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	  const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	  const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	  RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_flux = volume_fluxes_y[face_ooo];

	  // If only one material in cell_ooo
	  if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	    volume_fluxes_1_y[face_ooo] = vol_fraction_1 * volume_fluxes_y[face_ooo];
	    volume_fluxes_2_y[face_ooo] = vol_fraction_2 * volume_fluxes_y[face_ooo];

	  } else { 

	    const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	    const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	    const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	    const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	    const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	    const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

	    const RealType normal_x = interface_normal_x[cell_ooo];
	    const RealType normal_y = interface_normal_y[cell_ooo];

	    // Algebraic values of up and down displacement
	    RealType xmin = (cell_ooo - cell_ooo / nx * nx) * dx;
	    RealType xmax = (cell_ooo - cell_ooo / nx * nx + 1) * dx;
	    RealType xmin_lag_prev = xmin + in_vx[node_m1m1] * dt;
	    RealType xmax_lag_prev = xmax + in_vx[node_p1m1] * dt;
	    RealType xmin_lag_next = xmin + in_vx[node_m1p1] * dt;
	    RealType xmax_lag_next = xmax + in_vx[node_p1p1] * dt;
	    RealType xmin_prev = std::max(xmin_lag_prev,xmin);
	    RealType xmin_next = std::max(xmin_lag_next,xmin);
	    RealType xmax_prev = std::min(xmax_lag_prev,xmax);
	    RealType xmax_next = std::min(xmax_lag_next,xmax); 
	    RealType dx_lag_prev = xmax_prev - xmin_prev;
	    RealType dx_lag_next = xmax_next - xmin_next;

	    RealType delta_y_prev = volume_fluxes_y[prev_face] / dx_lag_prev;
	    RealType delta_y_next = volume_fluxes_y[next_face] / dx_lag_next;

	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType y = iy * dy;
	    RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * (y + ny * dy) + (1 - sign(vol_flux)) * (y + delta_y_prev) );
	    RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * (y + ny * dy + delta_y_next) + (1 - sign(vol_flux)) * y );
	    RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * xmin_next + (1 - sign(vol_flux)) * xmin_prev );
	    RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * xmax_next + (1 - sign(vol_flux)) * xmax_prev );

	    IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[cell_ooo][0] + pointB[cell_ooo][0]) / 2.;
	    RealType yI = (pointA[cell_ooo][1] + pointB[cell_ooo][1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      volume_fluxes_1_y[face_ooo] = vol_flux;
	      volume_fluxes_2_y[face_ooo] = 0.0;
	    } else {
	      volume_fluxes_1_y[face_ooo] = 0.0;
	      volume_fluxes_2_y[face_ooo] = vol_flux;
	    }	      
	  	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
		}

		volume_fluxes_1_y[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2_y[face_ooo] = vol_fraction_2_d * vol_flux;	      
	    }

	  }

	  volume_fluxes_1_y[face_sym] = volume_fluxes_1_y[face_ooo];
	  volume_fluxes_2_y[face_sym] = volume_fluxes_2_y[face_ooo];
	}  

  } //end Faces Y 


  // Corners
  //Xmin and Xmax
  for (index_t iy = 1; iy < ny; ++iy) {

      const index_t ix = 0;
      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * iy + nx;

      const RealType delta_x = in_vx[node_ooo] * dt;
      const RealType delta_y = in_vy[node_ooo] * dt;
      const RealType vol_flux = fabs(delta_x) * fabs(delta_y);
      
      if (vol_flux / (dx * dy) < epsilon_corner) {

	volume_fluxes_1_corner[node_ooo] = 0.0;
	volume_fluxes_2_corner[node_ooo] = 0.0;
	volume_fluxes_1_corner[node_sym] = 0.0;
	volume_fluxes_2_corner[node_sym] = 0.0;
	
      } else {
	
	index_t cell_ooo = NodeCellM1M1(node_sym,iy,nx);
	
	if (delta_x  < 0.0) {
	  if (delta_y > 0.0) {
	    cell_ooo = NodeCellP1M1(node_ooo, iy, nx);
	  } else {
	    cell_ooo = NodeCellP1P1(node_ooo, iy, nx);
	  }
	} else if (delta_y < 0.0) {
	  cell_ooo = NodeCellM1P1(node_sym, iy, nx);
	}

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

	// If only one material in cell_ooo
	if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	  volume_fluxes_1_corner[node_ooo] = vol_fraction_1 * vol_flux;
	  volume_fluxes_2_corner[node_ooo] = vol_fraction_2 * vol_flux;

	} else { 

	  const RealType normal_x = interface_normal_x[cell_ooo];
	  const RealType normal_y = interface_normal_y[cell_ooo];

	  // Compute the vertices of lagrangian dVol at node_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	  RealType pointA_d[2];
	  RealType pointB_d[2];
	  RealType x = 0.5 * ( (1 - sign(delta_x)) * ix * dx + (1 + sign(delta_x)) * ((nx + ix) * dx) );
	  RealType x_lag = x + delta_x;
	  RealType y = iy * dy;
	  RealType y_lag = y + delta_y;
	  RealType xmin_d = std::min(x,x_lag);
	  RealType xmax_d = std::max(x,x_lag);
	  RealType ymin_d = std::min(y,y_lag);
	  RealType ymax_d = std::max(y,y_lag);

	  IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	  // Default value : only one material in dVol
	  RealType x_inter = IntersectLineHorizontalLine(pointA[cell_ooo], pointB[cell_ooo], y_lag);

	  // If the interface is not horinzontal
	  if (x_inter < 1.0e19) {
	    // If the corner is on the right side of the lagrangian interface
	    if (x_inter < x_lag) {
	      if (normal_x < 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    // Else it is on the left side of the lagrangian interface
	    } else {
	      if (normal_x > 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    }
	  // If the interface is horizontal
	  } else {
	    if (normal_y * delta_y < 0.0) {
	      volume_fluxes_1_corner[node_ooo] = vol_flux;
	      volume_fluxes_2_corner[node_ooo] = 0.0;
	    } else { 
	      volume_fluxes_1_corner[node_ooo] = 0.0;
	      volume_fluxes_2_corner[node_ooo] = vol_flux;
	    }
	  }	  	    
	    
	  // If the interface intersects dVol, we compute the volumic fraction in dVol
	  if ( pointA_d[0] < 1.0e19 ) {

	      // Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at node_ooo
	      RealType vol_fraction_1_d, vol_fraction_2_d;

	      vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

	      if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 0.0;
		vol_fraction_2_d = 1.0;
		
	      } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 1.0;
		vol_fraction_2_d = 0.0;
		
	      } else {
	      
		RealType vol_fraction_11_d = vol_fraction_1_d;
		RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
	      }

	      volume_fluxes_1_corner[node_ooo] = vol_fraction_1_d * vol_flux;
	      volume_fluxes_2_corner[node_ooo] = vol_fraction_2_d * vol_flux;	    	       
	  }
	}
	  volume_fluxes_1_corner[node_sym] = volume_fluxes_1_corner[node_ooo];
	  volume_fluxes_2_corner[node_sym] = volume_fluxes_2_corner[node_ooo];
      }  

  } //end Xmin and Xmax

  //Ymin and Ymax
  for (index_t ix = 1; ix < nx; ++ix) {

      const index_t iy = 0;
      const index_t iy_sym = ny;
      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * iy_sym + ix;

      const RealType delta_x = in_vx[node_ooo] * dt;
      const RealType delta_y = in_vy[node_ooo] * dt;
      const RealType vol_flux = fabs(delta_x) * fabs(delta_y);
      
      if (vol_flux / (dx * dy) < epsilon_corner) {

	volume_fluxes_1_corner[node_ooo] = 0.0;
	volume_fluxes_2_corner[node_ooo] = 0.0;
	volume_fluxes_1_corner[node_sym] = 0.0;
	volume_fluxes_2_corner[node_sym] = 0.0;
	
      } else {
	
	index_t cell_ooo = NodeCellM1M1(node_sym,iy_sym,nx);
	
	if (delta_x  < 0.0) {
	  if (delta_y > 0.0) {
	    cell_ooo = NodeCellP1M1(node_sym, iy_sym, nx);
	  } else {
	    cell_ooo = NodeCellP1P1(node_ooo, iy, nx);
	  }
	} else if (delta_y < 0.0) {
	  cell_ooo = NodeCellM1P1(node_ooo, iy, nx);
	}

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

	// If only one material in cell_ooo
	if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	  volume_fluxes_1_corner[node_ooo] = vol_fraction_1 * vol_flux;
	  volume_fluxes_2_corner[node_ooo] = vol_fraction_2 * vol_flux;

	} else { 

	  const RealType normal_x = interface_normal_x[cell_ooo];
	  const RealType normal_y = interface_normal_y[cell_ooo];

	  // Compute the vertices of lagrangian dVol at node_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	  RealType pointA_d[2];
	  RealType pointB_d[2];
	  RealType y = 0.5 * ( (1 - sign(delta_y)) * iy * dy + (1 + sign(delta_y)) * ((ny + iy) * dy) );
	  RealType y_lag = y + delta_y;
	  RealType x = ix * dx;
	  RealType x_lag = x + delta_x;
	  RealType xmin_d = std::min(x,x_lag);
	  RealType xmax_d = std::max(x,x_lag);
	  RealType ymin_d = std::min(y,y_lag);
	  RealType ymax_d = std::max(y,y_lag);

	  IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	  // Default value : only one material in dVol
	  RealType x_inter = IntersectLineHorizontalLine(pointA[cell_ooo], pointB[cell_ooo], y_lag);

	  // If the interface is not horinzontal
	  if (x_inter < 1.0e19) {
	    // If the corner is on the right side of the lagrangian interface
	    if (x_inter < x_lag) {
	      if (normal_x < 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    // Else it is on the left side of the lagrangian interface
	    } else {
	      if (normal_x > 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    }
	  // If the interface is horizontal
	  } else {
	    if (normal_y * delta_y < 0.0) {
	      volume_fluxes_1_corner[node_ooo] = vol_flux;
	      volume_fluxes_2_corner[node_ooo] = 0.0;
	    } else { 
	      volume_fluxes_1_corner[node_ooo] = 0.0;
	      volume_fluxes_2_corner[node_ooo] = vol_flux;
	    }
	  }	  	    
	    
	  // If the interface intersects dVol, we compute the volumic fraction in dVol
	  if ( pointA_d[0] < 1.0e19 ) {

	      // Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at node_ooo
	      RealType vol_fraction_1_d, vol_fraction_2_d;

	      vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

	      if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 0.0;
		vol_fraction_2_d = 1.0;
		
	      } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 1.0;
		vol_fraction_2_d = 0.0;
		
	      } else {
	      
		RealType vol_fraction_11_d = vol_fraction_1_d;
		RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
	      }

	      volume_fluxes_1_corner[node_ooo] = vol_fraction_1_d * vol_flux;
	      volume_fluxes_2_corner[node_ooo] = vol_fraction_2_d * vol_flux;	    	       
	  }
	}
	  volume_fluxes_1_corner[node_sym] = volume_fluxes_1_corner[node_ooo];
	  volume_fluxes_2_corner[node_sym] = volume_fluxes_2_corner[node_ooo];
      }  

  } //end Ymin and Ymax

 // Corner Xmin, Ymin
  { const index_t ix = 0;
    const index_t iy = 0;
    const index_t iy_sym = ny;
    const index_t node_ooo = (nx + 1) * iy + ix;
    const index_t node_sym_x = (nx + 1) * iy + nx;
    const index_t node_sym_y = (nx + 1) * iy_sym + ix;
    const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;      

    const RealType delta_x = in_vx[node_ooo] * dt;
    const RealType delta_y = in_vy[node_ooo] * dt;
    const RealType vol_flux = fabs(delta_x) * fabs(delta_y);
	
    if (vol_flux / (dx * dy) < epsilon_corner) {

      volume_fluxes_1_corner[node_ooo] = 0.0;
      volume_fluxes_2_corner[node_ooo] = 0.0;
      volume_fluxes_1_corner[node_sym_x] = 0.0;
      volume_fluxes_2_corner[node_sym_x] = 0.0;
      volume_fluxes_1_corner[node_sym_y] = 0.0;
      volume_fluxes_2_corner[node_sym_y] = 0.0;	
      volume_fluxes_1_corner[node_sym_xy] = 0.0;
      volume_fluxes_2_corner[node_sym_xy] = 0.0;

    } else {
	
      index_t cell_ooo = NodeCellM1M1(node_sym_xy,iy_sym,nx);
	
      if (delta_x  < 0.0) {
	if (delta_y > 0.0) {
	  cell_ooo = NodeCellP1M1(node_sym_y, iy_sym, nx);
	} else {
	  cell_ooo = NodeCellP1P1(node_ooo, iy, nx);
	}
      } else if (delta_y < 0.0) {
	cell_ooo = NodeCellM1P1(node_sym_x, iy, nx);
      }

      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      // If only one material in cell_ooo
      if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	volume_fluxes_1_corner[node_ooo] = vol_fraction_1 * vol_flux;
	volume_fluxes_2_corner[node_ooo] = vol_fraction_2 * vol_flux;

      } else { 

	const RealType normal_x = interface_normal_x[cell_ooo];
	const RealType normal_y = interface_normal_y[cell_ooo];

	// Compute the vertices of lagrangian dVol at node_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	RealType pointA_d[2];
	RealType pointB_d[2];
	RealType y = 0.5 * ( (1 - sign(delta_y)) * iy * dy + (1 + sign(delta_y)) * ((ny + iy) * dy) );
	RealType y_lag = y + delta_y;
	RealType x =  0.5 * ( (1 - sign(delta_x)) * ix * dx + (1 + sign(delta_x)) * ((nx + ix) * dx) );
	RealType x_lag = x + delta_x;
	RealType xmin_d = std::min(x,x_lag);
	RealType xmax_d = std::max(x,x_lag);
	RealType ymin_d = std::min(y,y_lag);
	RealType ymax_d = std::max(y,y_lag);

	IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	// Default value : only one material in dVol
	RealType x_inter = IntersectLineHorizontalLine(pointA[cell_ooo], pointB[cell_ooo], y_lag);

	// If the interface is not horinzontal
	if (x_inter < 1.0e19) {
	  // If the corner is on the right side of the lagrangian interface
	  if (x_inter < x_lag) {
	    if (normal_x < 0.0) {
	      volume_fluxes_1_corner[node_ooo] = vol_flux;
	      volume_fluxes_2_corner[node_ooo] = 0.0;
	    } else { 
	      volume_fluxes_1_corner[node_ooo] = 0.0;
	      volume_fluxes_2_corner[node_ooo] = vol_flux;
	    }
	    // Else it is on the left side of the lagrangian interface
	  } else {
	    if (normal_x > 0.0) {
	      volume_fluxes_1_corner[node_ooo] = vol_flux;
	      volume_fluxes_2_corner[node_ooo] = 0.0;
	    } else { 
	      volume_fluxes_1_corner[node_ooo] = 0.0;
	      volume_fluxes_2_corner[node_ooo] = vol_flux;
	    }
	  }
	  // If the interface is horizontal
	} else {
	  if (normal_y * delta_y < 0.0) {
	    volume_fluxes_1_corner[node_ooo] = vol_flux;
	    volume_fluxes_2_corner[node_ooo] = 0.0;
	  } else { 
	    volume_fluxes_1_corner[node_ooo] = 0.0;
	    volume_fluxes_2_corner[node_ooo] = vol_flux;
	  }
	}	  	    
	    
	// If the interface intersects dVol, we compute the volumic fraction in dVol
	if ( pointA_d[0] < 1.0e19 ) {

	  // Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at node_ooo
	  RealType vol_fraction_1_d, vol_fraction_2_d;

	  vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

	  if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
	    vol_fraction_1_d = 0.0;
	    vol_fraction_2_d = 1.0;
		
	  } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
	    vol_fraction_1_d = 1.0;
	    vol_fraction_2_d = 0.0;
		
	  } else {
	      
	    RealType vol_fraction_11_d = vol_fraction_1_d;
	    RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
	    vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
	    vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

	    assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
	  }

	  volume_fluxes_1_corner[node_ooo] = vol_fraction_1_d * vol_flux;
	  volume_fluxes_2_corner[node_ooo] = vol_fraction_2_d * vol_flux;	    	       
	}
      }

      volume_fluxes_1_corner[node_sym_x] = volume_fluxes_1_corner[node_ooo];
      volume_fluxes_2_corner[node_sym_x] = volume_fluxes_2_corner[node_ooo];
      volume_fluxes_1_corner[node_sym_y] = volume_fluxes_1_corner[node_ooo];
      volume_fluxes_2_corner[node_sym_y] = volume_fluxes_2_corner[node_ooo];        
      volume_fluxes_1_corner[node_sym_xy] = volume_fluxes_1_corner[node_ooo];
      volume_fluxes_2_corner[node_sym_xy] = volume_fluxes_2_corner[node_ooo];
    }
  } //end Corner Xmin, Ymin

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerPeriodicBoundary




    // Rayleigh Taylor 

void ReconstructNormalInterfaceBoundaryRt(// in
					  index_t nx,
					  index_t ny,
					  RealType dx,
					  RealType dy,
					  const RealType* in_variable,
					  // out
					  RealType* interface_normal_x,
					  RealType* interface_normal_y) {


  // x min
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellOM1(cell_m1oo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_m1oo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // x max
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = nx - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_p1oo, nx);
      const int cell_p1p1 = CellCellOP1(cell_p1oo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
    
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // y min
  for (index_t ix = 1; ix < nx - 1; ++ix) {
 
    index_t iy = 0;
    index_t iy_sym = 1;

    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oom1, nx);
      const int cell_p1m1 = CellCellP1O(cell_oom1, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // y max
  for (index_t ix = 1; ix < nx - 1; ++ix) {
    
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oop1, nx);
      const int cell_m1p1 = CellCellM1O(cell_oop1, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}

    }
  }

  // corner x min, y min
  { index_t ix = 0;
    index_t iy = 0;
    index_t iy_sym = 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_oom1, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_m1oo, nx);
      const int cell_m1m1 = cell_m1p1;
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x max, y min
  { index_t ix = nx - 1;
    index_t iy = 0;
    index_t iy_sym = 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = iy_sym * nx + ix;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oom1, nx);
      const int cell_p1p1 = CellCellOP1(cell_p1oo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      const int cell_p1m1 = cell_p1p1;
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x min, y max
  { index_t ix = 0;
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo + nx - 1;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellOM1(cell_m1oo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oop1, nx);
      const int cell_m1p1 = cell_m1m1; 
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x max, y max
  { index_t ix = nx - 1;
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = iy * nx;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = iy_sym * nx + ix;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_p1oo, nx);
      const int cell_p1p1 = cell_p1m1;
      const int cell_m1p1 = CellCellM1O(cell_oop1, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

} // end ReconstructNormalInterfaceBoundaryRt



    // Wall boundary conditions

void ReconstructNormalInterfaceBoundaryWall(// in
					    index_t nx,
					    index_t ny,
					    RealType dx,
					    RealType dy,
					    const RealType* in_variable,
					    // out
					    RealType* interface_normal_x,
					    RealType* interface_normal_y) {


  // x min
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = 0;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellOM1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
      
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}      
    }
  }

  // x max
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    
    index_t ix = nx - 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = cell_oooo;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_oooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // y min
  for (index_t ix = 1; ix < nx - 1; ++ix) {
 
    index_t iy = 0;
    index_t iy_sym = 1;

    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = cell_oooo;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"
  
	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}      
    }
  }

  // y max
  for (index_t ix = 1; ix < nx - 1; ++ix) {
    
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = cell_oooo;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1O(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}      
    }
  }

  // corner x min, y min
  { index_t ix = 0;
    index_t iy = 0;
    index_t iy_sym = 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = cell_oooo;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1O(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
      const int cell_m1p1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = cell_oooo;
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x max, y min
  { index_t ix = nx - 1;
    index_t iy = 0;
    index_t iy_sym = 1;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = cell_oooo;
      const int cell_oom1 = cell_oooo;
      const int cell_oop1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1m1 = CellCellM1O(cell_oooo, nx);
      const int cell_p1p1 = CellCellOP1(cell_oooo, nx);
      const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
      const int cell_p1m1 = cell_oooo;
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x min, y max
  { index_t ix = 0;
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = cell_oooo;
      const int cell_p1oo = CellCellP1O(cell_oooo, nx);
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = cell_oooo;
      const int cell_m1m1 = CellCellOM1(cell_oooo, nx);
      const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
      const int cell_p1p1 = CellCellP1O(cell_oooo, nx);
      const int cell_m1p1 = cell_oooo; 
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

  // corner x max, y max
  { index_t ix = nx - 1;
    index_t iy = ny - 1;
    index_t iy_sym = ny - 2;
    
    const int cell_oooo = iy * nx + ix;

    if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

      interface_normal_x[cell_oooo] = 0.0;
      interface_normal_y[cell_oooo] = 0.0;
      
    } else {

      const int cell_m1oo = CellCellM1O(cell_oooo, nx);
      const int cell_p1oo = cell_oooo;
      const int cell_oom1 = CellCellOM1(cell_oooo, nx);
      const int cell_oop1 = cell_oooo;
      const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
      const int cell_p1m1 = CellCellOM1(cell_oooo, nx);
      const int cell_p1p1 = cell_oooo;
      const int cell_m1p1 = CellCellM1O(cell_oooo, nx);
      
#include "reconstruct_variable_gradient_Youngs.h"

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
    }
  }

} // end ReconstructNormalInterfaceBoundaryWall



void ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundaryWall(// in
								     index_t nx, 
								     index_t ny, 
								     RealType dx,
								     RealType dy,
								     const RealType* RESTRICT volume_fluxes,
								     const RealType* RESTRICT directional_lagrangian_volume,
								     const RealType* RESTRICT cell_volumic_fraction,
								     const RealType* RESTRICT interface_normal_x,
								     const RealType* RESTRICT interface_normal_y,
								     // out
								     RealType* RESTRICT volume_fluxes_1,
								     RealType* RESTRICT volume_fluxes_2) {


  for (index_t iy = 0; iy < ny; ++iy) {

    // x min
    { index_t ix = 0;

      index_t face_ooo = ((nx + 1) * iy) + ix; 

      volume_fluxes_1[face_ooo] = 0;
      volume_fluxes_2[face_ooo] = 0;
    }

    // x max
    { index_t ix = nx;

      index_t face_ooo = ((nx + 1) * iy) + ix; 

      volume_fluxes_1[face_ooo] = 0;
      volume_fluxes_2[face_ooo] = 0;
    }   

  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundaryWall



void ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundaryWall(// in
								     index_t nx, 
								     index_t ny, 
								     RealType dx,
								     RealType dy,
								     const RealType* RESTRICT volume_fluxes_y,
								     const RealType* RESTRICT directional_lagrangian_volume,
								     const RealType* RESTRICT cell_volumic_fraction,
								     const RealType* RESTRICT interface_normal_x,
								     const RealType* RESTRICT interface_normal_y,
								     // out
								     RealType* RESTRICT volume_fluxes_1_y,
								     RealType* RESTRICT volume_fluxes_2_y) {


  for (index_t ix = 0; ix < nx; ++ix) {

    // y min
    { index_t iy = 0;

      index_t face_ooo = nx * iy + ix; 

      volume_fluxes_1_y[face_ooo] = 0;
      volume_fluxes_2_y[face_ooo] = 0;
    }

    // y max
    { index_t iy = ny;

      index_t face_ooo = nx * iy + ix; 

      volume_fluxes_1_y[face_ooo] = 0;
      volume_fluxes_2_y[face_ooo] = 0;
    }   

  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundaryWall


void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerBoundaryWall(// in
									    index_t nx, 
									    index_t ny, 
									    RealType dx,
									    RealType dy,
									    RealType dt,
									    const RealType* RESTRICT volume_fluxes,
									    const RealType* RESTRICT volume_fluxes_y,
									    const RealType* RESTRICT directional_lagrangian_volume,
									    const RealType* RESTRICT cell_volumic_fraction,
									    const RealType* RESTRICT interface_normal_x,
									    const RealType* RESTRICT interface_normal_y,
									    const RealType* RESTRICT in_vx,
									    const RealType* RESTRICT in_vy,
									    // out
									    RealType* RESTRICT volume_fluxes_1,
									    RealType* RESTRICT volume_fluxes_1_y,
									    RealType* RESTRICT volume_fluxes_1_corner,
									    RealType* RESTRICT volume_fluxes_2,
									    RealType* RESTRICT volume_fluxes_2_y,
									    RealType* RESTRICT volume_fluxes_2_corner) {

  // Faces X
  for (index_t iy = 0; iy < ny; ++iy) {
 
        const index_t ix = 0;
	const index_t face_ooo = ((nx + 1) * iy) + ix;
	const index_t face_sym = ((nx + 1) * iy) + nx;

	volume_fluxes_1[face_ooo] = 0.0;
	volume_fluxes_2[face_ooo] = 0.0;
	volume_fluxes_1[face_sym] = 0.0;
	volume_fluxes_2[face_sym] = 0.0;
  }

  // Faces Y
  for (index_t ix = 0; ix < nx; ++ix) {

	const index_t iy = 0;
	const index_t iy_sym = ny;
	const index_t face_ooo = nx * iy + ix;
	const index_t face_sym = nx * iy_sym + ix;

	volume_fluxes_1_y[face_ooo] = 0.0;
	volume_fluxes_2_y[face_ooo] = 0.0;
	volume_fluxes_1_y[face_sym] = 0.0;
	volume_fluxes_2_y[face_sym] = 0.0;
  }


  // Corners
  //Xmin and Xmax
  for (index_t iy = 1; iy < ny; ++iy) {

      const index_t ix = 0;
      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * iy + nx;

      volume_fluxes_1_corner[node_ooo] = 0.0;
      volume_fluxes_2_corner[node_ooo] = 0.0;
      volume_fluxes_1_corner[node_sym] = 0.0;
      volume_fluxes_2_corner[node_sym] = 0.0;
  }

  //Ymin and Ymax
  for (index_t ix = 1; ix < nx; ++ix) {

      const index_t iy = 0;
      const index_t iy_sym = ny;
      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * iy_sym + ix;

      volume_fluxes_1_corner[node_ooo] = 0.0;
      volume_fluxes_2_corner[node_ooo] = 0.0;
      volume_fluxes_1_corner[node_sym] = 0.0;
      volume_fluxes_2_corner[node_sym] = 0.0;
  }

 // Corner Xmin, Ymin
  {   const index_t ix = 0;
      const index_t iy = 0;
      const index_t iy_sym = ny;
      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym_x = (nx + 1) * iy + nx;
      const index_t node_sym_y = (nx + 1) * iy_sym + ix;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;      

      volume_fluxes_1_corner[node_ooo] = 0.0;
      volume_fluxes_2_corner[node_ooo] = 0.0;
      volume_fluxes_1_corner[node_sym_x] = 0.0;
      volume_fluxes_2_corner[node_sym_x] = 0.0;
      volume_fluxes_1_corner[node_sym_y] = 0.0;
      volume_fluxes_2_corner[node_sym_y] = 0.0;	
      volume_fluxes_1_corner[node_sym_xy] = 0.0;
      volume_fluxes_2_corner[node_sym_xy] = 0.0;
  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerBoundaryWall

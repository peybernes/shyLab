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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
   
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
    }
  }

} // end ReconstructNormalInterfacePeriodicBoundary




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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
    
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;

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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
    
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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
      
      interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
      interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
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

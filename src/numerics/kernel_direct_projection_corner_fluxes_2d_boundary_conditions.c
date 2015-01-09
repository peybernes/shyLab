// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_corner_fluxes_2d_boundary_conditions.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


// Periodic Boundary


void ReconstructGradientDiagPeriodicBoundary(index_t nx,
					     index_t ny,
					     RealType dx,
					     RealType dy,
					     RealType dt,
					     const RealType* RESTRICT predicted_u,
					     const RealType* RESTRICT predicted_v,
					     const RealType* RESTRICT lag_variable,
					     RealType* RESTRICT gradient_variable_diag,
					     RealType* RESTRICT gradient_variable_antidiag) {

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

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store  
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
  
#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store  
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

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
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

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
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

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
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

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
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

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
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

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
  }

} //end ReconstructGradientDiagPeriodicBoundary


void ReconstructMassFluxOrder2XCornerFluxesPeriodicBoundary(index_t nx, 
							    index_t ny, 
							    index_t halo_width,
							    const RealType dx,
							    const RealType dy,
							    const RealType dt,
							    const RealType* RESTRICT predicted_u,
							    const RealType* RESTRICT predicted_v,
							    const RealType* RESTRICT volume_fluxes_tot,
							    const RealType* RESTRICT volume_fluxes,
							    const RealType* RESTRICT cell_density,
							    const RealType* RESTRICT cell_density_gradient,
							    RealType* RESTRICT reconstructed_density_faces_x,
							    RealType* RESTRICT mass_flux) {
  
  //x min only //
  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
    const index_t node_om1 = FaceXNodeM1(face_ooo, iy, nx);
    const index_t node_op1 = FaceXNodeP1(face_ooo, iy, nx);
    const index_t face_sym = face_ooo + nx;

    const index_t prev_cell = FaceXCellM1(face_sym, iy, nx);
    const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
    const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
    const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);

    // --.h file
    const RealType prev_cell_variable = cell_density[prev_cell];
    const RealType next_cell_variable = cell_density[next_cell];
    const RealType prev_cell_gradient = cell_density_gradient[prev_cell];
    const RealType next_cell_gradient = cell_density_gradient[next_cell];

    const RealType vol_flux = volume_fluxes[face_ooo];
    const RealType vol_flux_tot = volume_fluxes_tot[face_ooo];
    const RealType vol_flux_m1o = volume_fluxes_tot[face_m1o];
    const RealType vol_flux_p1o = volume_fluxes_tot[face_p1o];

    const RealType half = 0.5;
      
    const RealType y_om1 = (node_om1 / (nx + 1)) * dy;
    const RealType y_op1 = y_om1 + dy;
    const RealType y_om1_lag = y_om1 + predicted_v[node_om1] * dt;
    const RealType y_op1_lag = y_op1 + predicted_v[node_op1] * dt;
    const RealType dy_lag = std::min(y_op1,y_op1_lag) - std::max(y_om1,y_om1_lag);          
    const RealType one_over_dy_lag = 1.0 / dy_lag;                // 1 div
    
    const RealType dx_lag_prev_corrected = dx - vol_flux_m1o *  one_over_dy_lag;   // 1 fma
    const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o * one_over_dy_lag;	
    const RealType prev_cell_variable_o2 = prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
    const RealType mass_flux_ooo =
      (0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
      (0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
    mass_flux[face_ooo] = mass_flux_ooo;
    mass_flux[face_sym] = mass_flux_ooo;
    reconstructed_density_faces_x[face_ooo] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);    
    reconstructed_density_faces_x[face_sym] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);

  }

} //end  ReconstructMassFluxOrder2XCornerFluxesBoundary


void ReconstructMassFluxOrder2YCornerFluxesPeriodicBoundary(index_t nx, 
							    index_t ny, 
							    index_t halo_width,
							    const RealType dx,
							    const RealType dy,
							    const RealType dt,
							    const RealType* RESTRICT predicted_u,
							    const RealType* RESTRICT predicted_v,
							    const RealType* RESTRICT volume_fluxes_tot,
							    const RealType* RESTRICT volume_fluxes,
							    const RealType* RESTRICT cell_density,
							    const RealType* RESTRICT cell_density_gradient,
							    RealType* RESTRICT reconstructed_density_faces_y,
							    RealType* RESTRICT mass_flux) {

  //y min only 
  const index_t iy = 0;
  const index_t iy_sym = ny;
  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect
       
    const index_t face_ooo = (nx * iy) + ix;
    const index_t node_m1o = FaceYNodeM1(face_ooo, iy, nx);
    const index_t node_p1o = FaceYNodeP1(face_ooo, iy, nx);
    const index_t face_sym = iy_sym * nx + ix;

    const index_t prev_cell = FaceYCellM1(face_sym, iy_sym, nx);
    const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
    const index_t face_om1 = CellFaceOM1( prev_cell, iy_sym, nx);
    const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
    const RealType prev_cell_variable = cell_density[prev_cell];
    const RealType next_cell_variable = cell_density[next_cell];
    const RealType prev_cell_gradient = cell_density_gradient[prev_cell];
    const RealType next_cell_gradient = cell_density_gradient[next_cell];

    const RealType vol_flux = volume_fluxes[face_ooo];
    const RealType vol_flux_tot = volume_fluxes_tot[face_ooo];
    const RealType vol_flux_om1 = volume_fluxes_tot[face_om1];
    const RealType vol_flux_op1 = volume_fluxes_tot[face_op1];
      
    const RealType half = 0.5;

    const RealType x_m1o = (node_m1o - (node_m1o / (nx + 1)) * (nx + 1)) * dx;
    const RealType x_p1o = x_m1o + dx;
    const RealType x_m1o_lag = x_m1o + predicted_u[node_m1o] * dt;
    const RealType x_p1o_lag = x_p1o + predicted_u[node_p1o] * dt;
    const RealType dx_lag = std::min(x_p1o,x_p1o_lag) - std::max(x_m1o,x_m1o_lag); 
    const RealType one_over_dx_lag = 1.0 / dx_lag;                         
    
    const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx_lag; // 1fma
    const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx_lag;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
    const RealType mass_flux_ooo =
      (0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
      (0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
    mass_flux[face_ooo] = mass_flux_ooo;
    mass_flux[face_sym] = mass_flux_ooo;
    reconstructed_density_faces_y[face_ooo] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);
    reconstructed_density_faces_y[face_sym] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);

  }
  
} //end  ReconstructMassFluxOrder2YCornerFluxesBoundary


void ReconstructMassFluxOrder2CornerPeriodicBoundary(//in
						     index_t nx, 
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* RESTRICT directional_lagrangian_density,
						     const RealType* RESTRICT reconstructed_density_faces_x,
						     const RealType* RESTRICT reconstructed_density_faces_y,
						     const RealType* RESTRICT volume_fluxes_corner, 
						     const RealType* RESTRICT in_vx, 
						     const RealType* RESTRICT in_vy,
						     //out
						     RealType* RESTRICT mass_corner_fluxes) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      const index_t face_m1o = NodeFaceYM1OO(node_sym,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_ooo,iy,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, receiving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	  face_x_giv = face_op1;
	}
      }

      const RealType cell_density = directional_lagrangian_density[giving_cell];
      const RealType cell_density_sym = directional_lagrangian_density[receiving_cell];
      const RealType reconstructed_density_x = reconstructed_density_faces_x[face_x_giv];
      const RealType reconstructed_density_y = reconstructed_density_faces_y[face_y_giv];
      const RealType vol_flux_corner = volume_fluxes_corner[node_ooo];

      const RealType density_corner_o2 = (reconstructed_density_x + reconstructed_density_y) - cell_density;
      //const RealType density_corner_o2 = cell_density;
      //const RealType density_corner_lin = 0.5 * (cell_density + cell_density_sym);
      //const RealType density_corner_o2 = std::min(density_corner_lin,cell_density);
      const RealType dmass = density_corner_o2 * vol_flux_corner;
 
      mass_corner_fluxes[node_ooo] = dmass;
      mass_corner_fluxes[node_sym] = dmass;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      const index_t face_m1o = NodeFaceYM1OO(node_ooo,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_sym,iy_sym,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

     index_t giving_cell, receiving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	  face_x_giv = face_op1;
	}
      }

      const RealType cell_density = directional_lagrangian_density[giving_cell];
      const RealType cell_density_sym = directional_lagrangian_density[receiving_cell];
      const RealType reconstructed_density_x = reconstructed_density_faces_x[face_x_giv];
      const RealType reconstructed_density_y = reconstructed_density_faces_y[face_y_giv];
      const RealType vol_flux_corner = volume_fluxes_corner[node_ooo];

      const RealType density_corner_o2 = (reconstructed_density_x + reconstructed_density_y) - cell_density;
      //const RealType density_corner_o2 = cell_density;
      //const RealType density_corner_lin = 0.5 * (cell_density + cell_density_sym);
      //const RealType density_corner_o2 = std::min(density_corner_lin,cell_density);
      const RealType dmass = density_corner_o2 * vol_flux_corner;

      mass_corner_fluxes[node_ooo] = dmass;
      mass_corner_fluxes[node_sym] = dmass;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      const index_t cell_m1m1 = nx * (ny - 1) + nx - 1;
      const index_t cell_m1p1 = nx - 1;
      const index_t cell_p1m1 = nx * (ny - 1);
      const index_t cell_p1p1 = 0;

      const index_t face_m1o = NodeFaceYM1OO(node_sym_x,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_sym_y,ny,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;
	
      index_t giving_cell, receiving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	  face_x_giv = face_op1;
	}
      } 

      const RealType cell_density = directional_lagrangian_density[giving_cell];
      const RealType cell_density_sym = directional_lagrangian_density[receiving_cell];
      const RealType reconstructed_density_x = reconstructed_density_faces_x[face_x_giv];
      const RealType reconstructed_density_y = reconstructed_density_faces_y[face_y_giv];
      const RealType vol_flux_corner = volume_fluxes_corner[node_ooo];

      const RealType density_corner_o2 = (reconstructed_density_x + reconstructed_density_y) - cell_density;
      //const RealType density_corner_o2 = cell_density;
      //const RealType density_corner_lin = 0.5 * (cell_density + cell_density_sym);
      //const RealType density_corner_o2 = std::min(density_corner_lin,cell_density);
      const RealType dmass = density_corner_o2 * vol_flux_corner;
 
      mass_corner_fluxes[node_ooo] = dmass;
      mass_corner_fluxes[node_sym_x] = dmass;
      mass_corner_fluxes[node_sym_y] = dmass;
      mass_corner_fluxes[node_sym_xy] = dmass; 
  }

} //end ReconstructMassFluxOrder2CornerPeriodicBoundary



void ReconstructMassFluxOrder2DiagCornerPeriodicBoundary(//in
							 index_t nx, 
							 index_t ny,
							 RealType dx,
							 RealType dy,
							 RealType dt,
							 const RealType* RESTRICT predicted_u,
							 const RealType* RESTRICT predicted_v,
							 const RealType* RESTRICT directional_lagrangian_density,
							 const RealType* RESTRICT volume_fluxes_corner, 
							 const RealType* RESTRICT density_gradient_diag, 
							 const RealType* RESTRICT density_gradient_antidiag, 
							 //out
							 RealType* RESTRICT mass_corner_fluxes) {
 
 // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = predicted_u[node_ooo] * dt;
      RealType delta_y = predicted_v[node_ooo] * dt;

      if ((delta_x == 0.) || (delta_y == 0.)) {
	mass_corner_fluxes[node_ooo] = 0.;

      } else {

	index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;
	    node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy - 1, nx);

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	    node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;
	    node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy - 1, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	    node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	  }
	}

	const RealType half = 0.5;     
	const RealType cell_density = directional_lagrangian_density[giving_cell];
	const RealType vol_flux_corner = sign(delta_x) * volume_fluxes_corner[node_ooo];
	RealType dmass;

	if ((delta_x * delta_y) > 0.) {

	  const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_diag = density_gradient_diag[giving_cell];
	  const RealType prev_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_prev_corrected;
	  const RealType next_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_next_corrected;
	
	  const RealType mass_flux_corner =
	    (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	    (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	  dmass = fabs(mass_flux_corner);
      
	} else {

	  const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_antidiag = density_gradient_antidiag[giving_cell];
	  const RealType prev_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	  const RealType next_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	  const RealType mass_flux_corner =
	    (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	    (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	  dmass = fabs(mass_flux_corner);
	}

	assert((dmass > 0.) || (dmass == 0.));
	
	mass_corner_fluxes[node_ooo] = dmass;  
	mass_corner_fluxes[node_sym] = dmass;  
      } 
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = predicted_u[node_ooo] * dt;
      RealType delta_y = predicted_v[node_ooo] * dt;

      if ((delta_x == 0.) || (delta_y == 0.)) {
	mass_corner_fluxes[node_ooo] = 0.;

      } else {

	index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;
	    node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy_sym - 1, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy_sym - 1, nx);

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	    node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;
	    node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy_sym - 1, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy_sym - 1, nx);

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	    node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	  }
	}

	const RealType half = 0.5;     
	const RealType cell_density = directional_lagrangian_density[giving_cell];
	const RealType vol_flux_corner = sign(delta_x) * volume_fluxes_corner[node_ooo];
	RealType dmass;

	if ((delta_x * delta_y) > 0.) {

	  const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_diag = density_gradient_diag[giving_cell];
	  const RealType prev_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_prev_corrected;
	  const RealType next_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_next_corrected;
	
	  const RealType mass_flux_corner =
	    (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	    (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	  dmass = fabs(mass_flux_corner);
      
	} else {

	  const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_antidiag = density_gradient_antidiag[giving_cell];
	  const RealType prev_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	  const RealType next_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	  const RealType mass_flux_corner =
	    (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	    (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	  dmass = fabs(mass_flux_corner);
	}

	assert((dmass > 0.) || (dmass == 0.));
	
	mass_corner_fluxes[node_ooo] = dmass;  
	mass_corner_fluxes[node_sym] = dmass;
      }  
  } 

  // corners
  { index_t ix = 0; 
    index_t iy = 0;
    index_t ix_sym = nx;
    index_t iy_sym = ny;

    const index_t node_ooo = 0;
    const index_t node_sym_x = nx;
    const index_t node_sym_y = (nx + 1) * ny;
    const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

    const index_t cell_m1m1 = nx * (ny - 1) + nx - 1;
    const index_t cell_m1p1 = nx - 1;
    const index_t cell_p1m1 = nx * (ny - 1);
    const index_t cell_p1p1 = 0;

    RealType delta_x = predicted_u[node_ooo] * dt;
    RealType delta_y = predicted_v[node_ooo] * dt;

    if ((delta_x == 0.) || (delta_y == 0.)) {
      mass_corner_fluxes[node_ooo] = 0.;

    } else {

      index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
      if (delta_x > 0.0) {

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;
	  node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy_sym - 1, nx);
	  node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy_sym - 1, nx);

	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	  node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	  node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	}

      } else {

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;
	  node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy_sym - 1, nx);
	  node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy_sym - 1, nx);

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	  node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	  node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	}
      }

      const RealType half = 0.5;     
      const RealType cell_density = directional_lagrangian_density[giving_cell];
      const RealType vol_flux_corner = sign(delta_x) * volume_fluxes_corner[node_ooo];
      RealType dmass;

      if ((delta_x * delta_y) > 0.) {

	const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	const RealType cell_gradient_diag = density_gradient_diag[giving_cell];
	const RealType prev_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_prev_corrected;
	const RealType next_cell_density_corner = cell_density + half * cell_gradient_diag * diag_lag_next_corrected;
	
	const RealType mass_flux_corner =
	  (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	  (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	dmass = fabs(mass_flux_corner);
      
      } else {

	const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	const RealType cell_gradient_antidiag = density_gradient_antidiag[giving_cell];
	const RealType prev_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	const RealType next_cell_density_corner = cell_density + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	const RealType mass_flux_corner =
	  (half * (vol_flux_corner + fabs(vol_flux_corner)) * prev_cell_density_corner) + 
	  (half * (vol_flux_corner - fabs(vol_flux_corner)) * next_cell_density_corner);
       
	dmass = fabs(mass_flux_corner);
      }

      assert((dmass > 0.) || (dmass == 0.));
	
      mass_corner_fluxes[node_ooo] = dmass;  
      mass_corner_fluxes[node_sym_x] = dmass;  
      mass_corner_fluxes[node_sym_y] = dmass;  
      mass_corner_fluxes[node_sym_xy] = dmass;  
    } 
  }

} //end ReconstructMassFluxOrder2DiagCornerPeriodicBoundary



void ProjectMassCornerFluxesPeriodicBoundary(index_t nx, 
					     index_t ny,
					     RealType dt,
					     const RealType* RESTRICT mass_corner_fluxes, 
					     const RealType* RESTRICT in_vx, 
					     const RealType* RESTRICT in_vy, 
					     RealType* RESTRICT out_cell_mass) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, receiving_cell;

      if (delta_x > 0.0) {

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;

	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	}

      } else {

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	}
      }

      RealType dmass = mass_corner_fluxes[node_ooo];

      out_cell_mass[giving_cell] = out_cell_mass[giving_cell] - dmass; 
      out_cell_mass[receiving_cell] = out_cell_mass[receiving_cell] + dmass;     
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, receiving_cell;

      if (delta_x > 0.0) {

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;

	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	}

      } else {

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	}
      }

      RealType dmass = mass_corner_fluxes[node_ooo]; 

      out_cell_mass[giving_cell] = out_cell_mass[giving_cell] - dmass; 
      out_cell_mass[receiving_cell] = out_cell_mass[receiving_cell] + dmass;     
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      const index_t cell_m1m1 = nx * (ny - 1) + nx - 1;
      const index_t cell_m1p1 = nx - 1;
      const index_t cell_p1m1 = nx * (ny - 1);
      const index_t cell_p1p1 = 0;

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;
 
      index_t giving_cell, receiving_cell;

      if (delta_x > 0.0) {

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;

	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	}

      } else {

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	}
      }

      RealType dmass = mass_corner_fluxes[node_ooo];

      out_cell_mass[giving_cell] = out_cell_mass[giving_cell] - dmass; 
      out_cell_mass[receiving_cell] = out_cell_mass[receiving_cell] + dmass;     
  }

} //end ProjectMassCornerFluxesPeriodicBoundary



void ReconstructIntensiveVariableFluxOrder2XCornerFluxesPeriodicBoundary(index_t nx, 
									 index_t ny, 
									 index_t halo_width,
									 const RealType dx,
									 const RealType dy,
									 const RealType dt,
									 const RealType* RESTRICT predicted_u,
									 const RealType* RESTRICT predicted_v,
									 const RealType* RESTRICT volume_fluxes,
									 const RealType* RESTRICT mass_flux,
									 const RealType* RESTRICT cell_variable,
									 const RealType* RESTRICT cell_variable_gradient,
									 RealType* RESTRICT reconstructed_variable_faces_x,
									 RealType* RESTRICT variable_flux) {
  

  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
    const index_t node_om1 = FaceXNodeM1(face_ooo, iy, nx);
    const index_t node_op1 = FaceXNodeP1(face_ooo, iy, nx);
    const index_t face_sym = face_ooo + nx;

    const index_t prev_cell = FaceXCellM1(face_sym, iy, nx);
    const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
    const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
    const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);  
      
    const RealType prev_cell_variable = cell_variable[prev_cell];
    const RealType next_cell_variable = cell_variable[next_cell];
    const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];
    const RealType next_cell_gradient = cell_variable_gradient[next_cell];

    const RealType mass_flux_face = mass_flux[face_ooo];
      
    const RealType vol_flux = volume_fluxes[face_ooo];
    const RealType vol_flux_m1o = volume_fluxes[face_m1o];
    const RealType vol_flux_p1o = volume_fluxes[face_p1o];
      
    const RealType half = 0.5;

    const RealType y_om1 = (node_om1 / (nx + 1)) * dy;
    const RealType y_op1 = y_om1 + dy;
    const RealType y_om1_lag = y_om1 + predicted_v[node_om1] * dt;
    const RealType y_op1_lag = y_op1 + predicted_v[node_op1] * dt;
    const RealType dy_lag = std::min(y_op1,y_op1_lag) - std::max(y_om1,y_om1_lag);          
    const RealType one_over_dy_lag = 1.0 / dy_lag;                // 1 div
    
    const RealType dx_lag_prev_corrected = dx - vol_flux_m1o *  one_over_dy_lag;   // 1 fma
    const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o * one_over_dy_lag;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
    const RealType variable_flux_ooo =
      (0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
      (0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
    variable_flux[face_ooo] = variable_flux_ooo;
    variable_flux[face_sym] = variable_flux_ooo;
    reconstructed_variable_faces_x[face_ooo] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2); 
    reconstructed_variable_faces_x[face_sym] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2); 
    
  }

}// end ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary



void ReconstructIntensiveVariableFluxOrder2YCornerFluxesPeriodicBoundary(index_t nx, 
									 index_t ny, 
									 index_t halo_width,
									 const RealType dx,
									 const RealType dy,
									 const RealType dt,
									 const RealType* RESTRICT predicted_u,
									 const RealType* RESTRICT predicted_v,
									 const RealType* RESTRICT volume_fluxes,
									 const RealType* RESTRICT mass_flux,
									 const RealType* RESTRICT cell_variable,
									 const RealType* RESTRICT cell_variable_gradient,
									 RealType* RESTRICT reconstructed_variable_faces_y,
									 RealType* RESTRICT variable_flux) {


  index_t iy = 0;
  index_t iy_sym = ny;

  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect but boundary
      
    const index_t face_ooo = (iy * nx) + ix;
    const index_t node_m1o = FaceYNodeM1(face_ooo, iy, nx);
    const index_t node_p1o = FaceYNodeP1(face_ooo, iy, nx);
    const index_t face_sym = (iy_sym * nx) + ix;

    const index_t prev_cell = FaceYCellM1(face_sym, iy_sym, nx);
    const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
    const index_t face_om1 = CellFaceOM1( prev_cell, iy_sym, nx);
    const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
    const RealType prev_cell_variable = cell_variable[prev_cell];
    const RealType next_cell_variable = cell_variable[next_cell];
    const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];
    const RealType next_cell_gradient = cell_variable_gradient[next_cell];
 
    const RealType mass_flux_face = mass_flux[face_ooo];
      
    const RealType vol_flux = volume_fluxes[face_ooo];
    const RealType vol_flux_om1 = volume_fluxes[face_om1];
    const RealType vol_flux_op1 = volume_fluxes[face_op1];
      
    const RealType half = 0.5;

    const RealType x_m1o = (node_m1o - (node_m1o / (nx + 1)) * (nx + 1)) * dx;
    const RealType x_p1o = x_m1o + dx;
    const RealType x_m1o_lag = x_m1o + predicted_u[node_m1o] * dt;
    const RealType x_p1o_lag = x_p1o + predicted_u[node_p1o] * dt;
    const RealType dx_lag = std::min(x_p1o,x_p1o_lag) - std::max(x_m1o,x_m1o_lag); 
    const RealType one_over_dx_lag = 1.0 / dx_lag;                         
    
    const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx_lag; // 1fma
    const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx_lag;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
    const RealType variable_flux_ooo =
      (0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
      (0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
    variable_flux[face_ooo] = variable_flux_ooo;
    variable_flux[face_sym] = variable_flux_ooo;
    reconstructed_variable_faces_y[face_ooo] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2); 
    reconstructed_variable_faces_y[face_sym] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2); 

  }
  
}// end ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary


void ReconstructIntensiveVariableFluxOrder2CornerPeriodicBoundary(//in
								  index_t nx, 
								  index_t ny,
								  RealType dx,
								  RealType dy,
								  RealType dt,
								  const RealType* RESTRICT cell_variable,
								  const RealType* RESTRICT reconstructed_variable_faces_x,
								  const RealType* RESTRICT reconstructed_variable_faces_y,
								  const RealType* RESTRICT mass_corner_fluxes, 
								  const RealType* RESTRICT in_vx, 
								  const RealType* RESTRICT in_vy,
								  //out
								  RealType* RESTRICT variable_flux_corner) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      const index_t face_m1o = NodeFaceYM1OO(node_sym,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_ooo,iy,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  face_x_giv = face_op1;
	}
      }

      const RealType cell_variable_giv = cell_variable[giving_cell];
      const RealType reconstructed_variable_x = reconstructed_variable_faces_x[face_x_giv];
      const RealType reconstructed_variable_y = reconstructed_variable_faces_y[face_y_giv];
      const RealType mass_flux_corner = mass_corner_fluxes[node_ooo];

      const RealType variable_corner_o2 = (reconstructed_variable_x + reconstructed_variable_y) - cell_variable_giv;
      const RealType variable_flux_ooo = variable_corner_o2 * mass_flux_corner;

      variable_flux_corner[node_ooo] = variable_flux_ooo;
 
      variable_flux_corner[node_ooo] = variable_flux_ooo;
      variable_flux_corner[node_sym] = variable_flux_ooo;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      const index_t face_m1o = NodeFaceYM1OO(node_ooo,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_sym,iy_sym,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  face_x_giv = face_op1;
	}
      }

      const RealType cell_variable_giv = cell_variable[giving_cell];
      const RealType reconstructed_variable_x = reconstructed_variable_faces_x[face_x_giv];
      const RealType reconstructed_variable_y = reconstructed_variable_faces_y[face_y_giv];
      const RealType mass_flux_corner = mass_corner_fluxes[node_ooo];

      const RealType variable_corner_o2 = (reconstructed_variable_x + reconstructed_variable_y) - cell_variable_giv;
      const RealType variable_flux_ooo = variable_corner_o2 * mass_flux_corner;

      variable_flux_corner[node_ooo] = variable_flux_ooo;
      variable_flux_corner[node_sym] = variable_flux_ooo;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      const index_t cell_m1m1 = nx * (ny - 1) + nx - 1;
      const index_t cell_m1p1 = nx - 1;
      const index_t cell_p1m1 = nx * (ny - 1);
      const index_t cell_p1p1 = 0;

      const index_t face_m1o = NodeFaceYM1OO(node_sym_x,iy,nx);
      const index_t face_p1o = NodeFaceYP1OO(node_ooo,iy,nx);
      const index_t face_om1 = NodeFaceXOOM1(node_sym_y,ny,nx);
      const index_t face_op1 = NodeFaceXOOP1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      index_t giving_cell, face_x_giv, face_y_giv;
	
      if (delta_x > 0.0) {
	face_y_giv = face_m1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  face_x_giv = face_om1;
	    
	} else {
	  giving_cell = cell_m1p1;
	  face_x_giv = face_op1;
	}

      } else {
	face_y_giv = face_p1o;

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  face_x_giv = face_om1;

	} else {
	  giving_cell = cell_p1p1;
	  face_x_giv = face_op1;
	}
      }

      const RealType cell_variable_giv = cell_variable[giving_cell];
      const RealType reconstructed_variable_x = reconstructed_variable_faces_x[face_x_giv];
      const RealType reconstructed_variable_y = reconstructed_variable_faces_y[face_y_giv];
      const RealType mass_flux_corner = mass_corner_fluxes[node_ooo];

      const RealType variable_corner_o2 = (reconstructed_variable_x + reconstructed_variable_y) - cell_variable_giv;
      const RealType variable_flux_ooo = variable_corner_o2 * mass_flux_corner;

      variable_flux_corner[node_ooo] = variable_flux_ooo;
      variable_flux_corner[node_sym_x] = variable_flux_ooo;
      variable_flux_corner[node_sym_y] = variable_flux_ooo;
      variable_flux_corner[node_sym_xy] = variable_flux_ooo; 
  }

} //end ReconstructIntensiveVariableFluxOrder2CornerPeriodicBoundary



void ReconstructIntensiveVariableFluxOrder2DiagCornerPeriodicBoundary(//in
								      index_t nx, 
								      index_t ny,
								      RealType dx,
								      RealType dy,
								      RealType dt,
								      const RealType* RESTRICT predicted_u,
								      const RealType* RESTRICT predicted_v,
								      const RealType* RESTRICT cell_variable,
								      const RealType* RESTRICT mass_corner_fluxes, 
								      const RealType* RESTRICT variable_gradient_diag, 
								      const RealType* RESTRICT variable_gradient_antidiag, 
								      //out
								      RealType* RESTRICT variable_flux_corner) {


 // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = predicted_u[node_ooo] * dt;
      RealType delta_y = predicted_v[node_ooo] * dt;

      if ((delta_x == 0.) || (delta_y == 0.)) {
	variable_flux_corner[node_ooo] = 0.;

      } else {

	index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;
	    node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy - 1, nx);

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	    node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;
	    node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy - 1, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	    node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	  }
	}
	
	const RealType half = 0.5;     
	const RealType cell_variable_ooo = cell_variable[giving_cell];
	const RealType mass_flux_corner = sign(delta_x) * mass_corner_fluxes[node_ooo];
	RealType var_flux_corner;

	if ((delta_x * delta_y) > 0.) {

	  const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_diag = variable_gradient_diag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
      
	} else {

	  const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_antidiag = variable_gradient_antidiag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
	}

	assert((var_flux_corner > 0.) || (var_flux_corner == 0.));
	
	variable_flux_corner[node_ooo] = var_flux_corner;  
	variable_flux_corner[node_sym] = var_flux_corner;  
      } 
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = predicted_u[node_ooo] * dt;
      RealType delta_y = predicted_v[node_ooo] * dt;

      if ((delta_x == 0.) || (delta_y == 0.)) {
	variable_flux_corner[node_ooo] = 0.;

      } else {

	index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;
	    node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy_sym - 1, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy_sym - 1, nx);

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	    node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;
	    node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy_sym - 1, nx);
	    node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy_sym - 1, nx);

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	    node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	    node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	  }
	}

	const RealType half = 0.5;     
	const RealType cell_variable_ooo = cell_variable[giving_cell];
	const RealType mass_flux_corner = sign(delta_x) * mass_corner_fluxes[node_ooo];
	RealType var_flux_corner;

	if ((delta_x * delta_y) > 0.) {

	  const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_diag = variable_gradient_diag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
      
	} else {

	  const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_antidiag = variable_gradient_antidiag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
	}

	assert((var_flux_corner > 0.) || (var_flux_corner == 0.));
	
	variable_flux_corner[node_ooo] = var_flux_corner;  	
	variable_flux_corner[node_sym] = var_flux_corner;  	
      }  
  } 

  // corners
  { index_t ix = 0; 
    index_t iy = 0;
    index_t ix_sym = nx;
    index_t iy_sym = ny;

    const index_t node_ooo = 0;
    const index_t node_sym_x = nx;
    const index_t node_sym_y = (nx + 1) * ny;
    const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

    const index_t cell_m1m1 = nx * (ny - 1) + nx - 1;
    const index_t cell_m1p1 = nx - 1;
    const index_t cell_p1m1 = nx * (ny - 1);
    const index_t cell_p1p1 = 0;

    RealType delta_x = predicted_u[node_ooo] * dt;
    RealType delta_y = predicted_v[node_ooo] * dt;

    if ((delta_x == 0.) || (delta_y == 0.)) {
	variable_flux_corner[node_ooo] = 0.;

    } else {

      index_t giving_cell, receiving_cell, node_giv_m1m1, node_giv_m1p1, node_giv_p1m1, node_giv_p1p1;
	
      if (delta_x > 0.0) {

	if (delta_y > 0.0) {
	  giving_cell = cell_m1m1;
	  receiving_cell = cell_p1p1;
	  node_giv_m1m1 = CellNodeM1M1(cell_m1m1, iy_sym - 1, nx);
	  node_giv_p1p1 = CellNodeP1P1(cell_m1m1, iy_sym - 1, nx);

	} else {
	  giving_cell = cell_m1p1;
	  receiving_cell = cell_p1m1;
	  node_giv_m1p1 = CellNodeM1P1(cell_m1p1, iy, nx);
	  node_giv_p1m1 = CellNodeP1M1(cell_m1p1, iy, nx);
	}

      } else {

	if (delta_y > 0.0) {
	  giving_cell = cell_p1m1;
	  receiving_cell = cell_m1p1;
	  node_giv_m1p1 = CellNodeM1P1(cell_p1m1, iy_sym - 1, nx);
	  node_giv_p1m1 = CellNodeP1M1(cell_p1m1, iy_sym - 1, nx);

	} else {
	  giving_cell = cell_p1p1;
	  receiving_cell = cell_m1m1;
	  node_giv_m1m1 = CellNodeM1M1(cell_p1p1, iy, nx);
	  node_giv_p1p1 = CellNodeP1P1(cell_p1p1, iy, nx);
	}
      }
	
	const RealType half = 0.5;     
	const RealType cell_variable_ooo = cell_variable[giving_cell];
	const RealType mass_flux_corner = sign(delta_x) * mass_corner_fluxes[node_ooo];
	RealType var_flux_corner;

	if ((delta_x * delta_y) > 0.) {

	  const RealType diag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1m1] * dx + predicted_v[node_giv_m1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType diag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1p1] * dx + predicted_v[node_giv_p1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_diag = variable_gradient_diag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_diag * diag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
      
	} else {

	  const RealType antidiag_lag_prev_corrected = std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_m1p1] * dx + predicted_v[node_giv_m1p1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType antidiag_lag_next_corrected = - std::sqrt(dx * dx + dy * dy) - (predicted_u[node_giv_p1m1] * dx + predicted_v[node_giv_p1m1] * dy) * dt / std::sqrt(dx * dx + dy * dy);
	  const RealType cell_gradient_antidiag = variable_gradient_antidiag[giving_cell];
	  const RealType prev_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_prev_corrected;
	  const RealType next_cell_variable_corner = cell_variable_ooo + half * cell_gradient_antidiag * antidiag_lag_next_corrected;
	
	  const RealType var_flux_corner_ooo =
	    (half * (mass_flux_corner + fabs(mass_flux_corner)) * prev_cell_variable_corner) + 
	    (half * (mass_flux_corner - fabs(mass_flux_corner)) * next_cell_variable_corner);
       
	  var_flux_corner = fabs(var_flux_corner_ooo);
	}

	assert((var_flux_corner > 0.) || (var_flux_corner == 0.));
	
	variable_flux_corner[node_ooo] = var_flux_corner;  
	variable_flux_corner[node_sym_x] = var_flux_corner;  
	variable_flux_corner[node_sym_y] = var_flux_corner;  
	variable_flux_corner[node_sym_xy] = var_flux_corner;   
    } 
  }
} //end ReconstructIntensiveVariableFluxOrder2DiagCornerPeriodicBoundary



void MassProjectIntensiveVariableCornerFluxesPeriodicBoundary(index_t nx, 
							      index_t ny,
							      RealType dt,
							      const RealType* RESTRICT in_cell_variable,
							      const RealType* RESTRICT variable_flux_corner, 
							      const RealType* RESTRICT in_vx, 
							      const RealType* RESTRICT in_vy, 
							      const RealType* RESTRICT out_cell_mass,
							      RealType* RESTRICT out_cell_variable) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      if ((delta_x != 0.0) && (delta_y != 0.0)) {

	index_t giving_cell, receiving_cell;

	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	  }
	}

	const RealType variable_flux = variable_flux_corner[node_ooo];
	const RealType out_cell_mass_giv = out_cell_mass[giving_cell];
	const RealType out_cell_mass_rec = out_cell_mass[receiving_cell];

	RealType out_cell_variable_giv, out_cell_variable_rec;

	assert((0.0 < out_cell_mass_giv) || (0.0 == out_cell_mass_giv));
	if (out_cell_mass_giv == 0.0) {
	  out_cell_variable_giv = 0.0;
	} else {	
	  out_cell_variable_giv = out_cell_variable[giving_cell] - variable_flux / out_cell_mass_giv;
	}

	assert((0.0 < out_cell_mass_rec) || (0.0 == out_cell_mass_rec));
	if (out_cell_mass_rec == 0.0) {
	  out_cell_variable_rec = 0.0;
	} else {	
	  out_cell_variable_rec = out_cell_variable[receiving_cell] + variable_flux / out_cell_mass_rec;
	}

	out_cell_variable[giving_cell] = out_cell_variable_giv;
	out_cell_variable[receiving_cell] = out_cell_variable_rec;
      }
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1) * iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      if ((delta_x != 0.0) && (delta_y != 0.0)) {

	index_t giving_cell, receiving_cell;

	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	  }
	}

	const RealType variable_flux = variable_flux_corner[node_ooo];
	const RealType out_cell_mass_giv = out_cell_mass[giving_cell];
	const RealType out_cell_mass_rec = out_cell_mass[receiving_cell];

	RealType out_cell_variable_giv, out_cell_variable_rec;

	assert((0.0 < out_cell_mass_giv) || (0.0 == out_cell_mass_giv));
	if (out_cell_mass_giv == 0.0) {
	  out_cell_variable_giv = 0.0;
	} else {	
	  out_cell_variable_giv = out_cell_variable[giving_cell] - variable_flux / out_cell_mass_giv;
	}

	assert((0.0 < out_cell_mass_rec) || (0.0 == out_cell_mass_rec));
	if (out_cell_mass_rec == 0.0) {
	  out_cell_variable_rec = 0.0;
	} else {	
	  out_cell_variable_rec = out_cell_variable[receiving_cell] + variable_flux / out_cell_mass_rec;
	}

	out_cell_variable[giving_cell] = out_cell_variable_giv;
	out_cell_variable[receiving_cell] = out_cell_variable_rec;
      }
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;

      const index_t cell_m1m1 = nx * ny - 1;
      const index_t cell_m1p1 = nx - 1;
      const index_t cell_p1m1 = nx * (ny - 1);
      const index_t cell_p1p1 = 0;

      RealType delta_x = in_vx[node_ooo] * dt;
      RealType delta_y = in_vy[node_ooo] * dt;

      if ((delta_x != 0.0) && (delta_y != 0.0)) {

	index_t giving_cell, receiving_cell;

	if (delta_x > 0.0) {

	  if (delta_y > 0.0) {
	    giving_cell = cell_m1m1;
	    receiving_cell = cell_p1p1;

	  } else {
	    giving_cell = cell_m1p1;
	    receiving_cell = cell_p1m1;
	  }

	} else {

	  if (delta_y > 0.0) {
	    giving_cell = cell_p1m1;
	    receiving_cell = cell_m1p1;

	  } else {
	    giving_cell = cell_p1p1;
	    receiving_cell = cell_m1m1;
	  }
	}

	const RealType variable_flux = variable_flux_corner[node_ooo];
	const RealType out_cell_mass_giv = out_cell_mass[giving_cell];
	const RealType out_cell_mass_rec = out_cell_mass[receiving_cell];

	RealType out_cell_variable_giv, out_cell_variable_rec;

	assert((0.0 < out_cell_mass_giv) || (0.0 == out_cell_mass_giv));
	if (out_cell_mass_giv == 0.0) {
	  out_cell_variable_giv = 0.0;
	} else {	
	  out_cell_variable_giv = out_cell_variable[giving_cell] - variable_flux / out_cell_mass_giv;
	}

	assert((0.0 < out_cell_mass_rec) || (0.0 == out_cell_mass_rec));
	if (out_cell_mass_rec == 0.0) {
	  out_cell_variable_rec = 0.0;
	} else {	
	  out_cell_variable_rec = out_cell_variable[receiving_cell] + variable_flux / out_cell_mass_rec;
	}

	out_cell_variable[giving_cell] = out_cell_variable_giv;
	out_cell_variable[receiving_cell] = out_cell_variable_rec;
      }
  }

} //end MassProjectIntensiveVariableCornerFluxesPeriodicBoundary


void ProjectNodalIntensiveVariableCornerFluxesPeriodicBoundary(index_t nx, 
							       index_t ny, 
							       RealType dt,
							       index_t halo_width,
							       const RealType* RESTRICT out_cell_mass,
							       const RealType* RESTRICT in_variable,
							       const RealType* RESTRICT in_u,
							       const RealType* RESTRICT in_v,
							       const RealType* RESTRICT mass_corner_fluxes,
							       const index_t* RESTRICT sign_x_corner_fluxes,
							       const index_t* RESTRICT sign_y_corner_fluxes,
							       RealType* RESTRICT out_variable) {

 
  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_oooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_oooo + nx;

      const index_t node_m1oo = NodeNodeM1O(node_oooo, iy, nx);
      const index_t node_p1oo = NodeNodeP1O(node_oooo, iy, nx);
      const index_t node_oom1 = NodeNodeOM1(node_sym, iy, nx);
      const index_t node_oop1 = NodeNodeOP1(node_oooo, iy, nx);
      const index_t node_m1m1 = NodeNodeM1M1(node_sym, iy, nx);
      const index_t node_p1m1 = NodeNodeP1M1(node_oooo, iy, nx);
      const index_t node_m1p1 = NodeNodeM1P1(node_sym, iy, nx);
      const index_t node_p1p1 = NodeNodeP1P1(node_oooo, iy, nx);

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_sym,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_oooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_oooo,iy,nx);

#include "direct_projection_corner_fluxes_nodal_2d_computation.h"
      
    out_variable[node_oooo] = out_variable[node_oooo] + total_variable_flux;
    out_variable[node_sym] = out_variable[node_oooo]; 
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_oooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      const index_t node_m1oo = NodeNodeM1O(node_oooo, iy, nx);
      const index_t node_p1oo = NodeNodeP1O(node_oooo, iy, nx);
      const index_t node_oom1 = NodeNodeOM1(node_sym, iy_sym, nx);
      const index_t node_oop1 = NodeNodeOP1(node_oooo, iy, nx);
      const index_t node_m1m1 = NodeNodeM1M1(node_sym, iy_sym, nx);
      const index_t node_p1m1 = NodeNodeP1M1(node_sym, iy_sym, nx);
      const index_t node_m1p1 = NodeNodeM1P1(node_oooo, iy, nx);
      const index_t node_p1p1 = NodeNodeP1P1(node_oooo, iy, nx);

      const index_t cell_m1m1 = NodeCellM1M1(node_sym,iy_sym,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_oooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_sym,iy_sym,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_oooo,iy,nx);

#include "direct_projection_corner_fluxes_nodal_2d_computation.h"
      
    out_variable[node_oooo] = out_variable[node_oooo] + total_variable_flux;
    out_variable[node_sym] = out_variable[node_oooo]; 
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_oooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      const index_t node_m1oo = NodeNodeM1O(node_sym_x, iy, nx);
      const index_t node_p1oo = NodeNodeP1O(node_oooo, iy, nx);
      const index_t node_oom1 = NodeNodeOM1(node_sym_y, iy_sym, nx);
      const index_t node_oop1 = NodeNodeOP1(node_oooo, iy, nx);
      const index_t node_m1m1 = NodeNodeM1M1(node_sym_xy, iy_sym, nx);
      const index_t node_p1m1 = NodeNodeP1M1(node_sym_y, iy_sym, nx);
      const index_t node_m1p1 = NodeNodeM1P1(node_sym_x, iy, nx);
      const index_t node_p1p1 = NodeNodeP1P1(node_oooo, iy, nx);

      const index_t cell_m1m1 = nx * ny - 1;
      const index_t cell_m1p1 = nx - 1;
      const index_t cell_p1m1 = nx * (ny - 1);
      const index_t cell_p1p1 = 0;

#include "direct_projection_corner_fluxes_nodal_2d_computation.h"
      
    out_variable[node_oooo] = out_variable[node_oooo] + total_variable_flux;
    out_variable[node_sym_x] = out_variable[node_oooo]; 
    out_variable[node_sym_y] = out_variable[node_oooo]; 
    out_variable[node_sym_xy] = out_variable[node_oooo]; 
  }

} //end ProjectNodalIntensiveVariableCornerFluxesPeriodicBoundary






     // Boundary Wall


void ReconstructGradientDiagBoundaryWall(index_t nx,
					 index_t ny,
					 RealType dx,
					 RealType dy,
					 RealType dt,
					 const RealType* RESTRICT predicted_u,
					 const RealType* RESTRICT predicted_v,
					 const RealType* RESTRICT lag_variable,
					 RealType* RESTRICT gradient_variable_diag,
					 RealType* RESTRICT gradient_variable_antidiag) {

  //Xmin and Xmax 
  for (index_t iy = 1; iy < ny - 1; ++iy) {

    { const index_t ix = 0;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellOM1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellOP1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);;

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store  
    }

    { const index_t ix = nx - 1;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellOM1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellOP1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);
  
#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store  
    }   
  }

  //Ymin and Ymax 
  for (index_t ix = 1; ix < nx - 1; ++ix) {

    { const index_t iy = 0;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellM1O(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellP1O(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
    }

    { const index_t iy = ny - 1;

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1O(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1O(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy, nx);

#include "reconstruct_variable_gradient_diag.h"

      gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
      gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
    }
  }
   
  //Corners

  //Xmin Ymin
  { const index_t ix = 0;
    const index_t iy = 0;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_m1m1 = cell_ooo;
    const index_t cell_m1p1 = CellCellOP1(cell_ooo, nx);
    const index_t cell_p1m1 = CellCellP1O(cell_ooo, nx);
    const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
  }

  //Xmin Ymax
  { const index_t ix = nx - 1;
    const index_t iy = 0;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_m1m1 = CellCellM1O(cell_ooo, nx);
    const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
    const index_t cell_p1m1 = cell_ooo;
    const index_t cell_p1p1 = CellCellOP1(cell_ooo, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
  }

  //Xmax Ymin
  { const index_t ix = 0;
    const index_t iy = ny - 1;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_m1m1 = CellCellOM1(cell_ooo, nx);
    const index_t cell_m1p1 = cell_ooo;
    const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
    const index_t cell_p1p1 = CellCellP1O(cell_ooo, nx);

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy, nx);

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
  }

  //Xmax Ymax
  { const index_t ix = nx - 1;
    const index_t iy = ny - 1;

    const index_t cell_ooo = iy * nx + ix;
    const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
    const index_t cell_m1p1 = CellCellM1O(cell_ooo, nx);
    const index_t cell_p1m1 = CellCellOM1(cell_ooo, nx);
    const index_t cell_p1p1 = cell_ooo;

    const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
    const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
    const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy, nx);
    const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy, nx);

#include "reconstruct_variable_gradient_diag.h"

    gradient_variable_diag[cell_ooo] = limited_grad_variable_diag; // 1 store
    gradient_variable_antidiag[cell_ooo] = limited_grad_variable_antidiag; // 1 store 
  }

} //end ReconstructGradientDiagBoundaryWall


void ReconstructMassFluxOrder2XCornerFluxesBoundaryWall(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dx,
							const RealType dy,
							const RealType dt,
							const RealType* RESTRICT predicted_u,
							const RealType* RESTRICT predicted_v,
							const RealType* RESTRICT volume_fluxes_tot,
							const RealType* RESTRICT volume_fluxes,
							const RealType* RESTRICT cell_density,
							const RealType* RESTRICT cell_density_gradient,
							RealType* RESTRICT reconstructed_density_faces_x,
							RealType* RESTRICT mass_flux) {
  
  //x min only //
  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
    const index_t face_sym = face_ooo + nx;
      
    mass_flux[face_ooo] = 0.;
    mass_flux[face_sym] = 0.;
    reconstructed_density_faces_x[face_ooo] = 0.;    
    reconstructed_density_faces_x[face_sym] = 0.;

  }

} //end  ReconstructMassFluxOrder2XCornerFluxesBoundaryWall


void ReconstructMassFluxOrder2YCornerFluxesBoundaryWall(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dx,
							const RealType dy,
							const RealType dt,
							const RealType* RESTRICT predicted_u,
							const RealType* RESTRICT predicted_v,
							const RealType* RESTRICT volume_fluxes_tot,
							const RealType* RESTRICT volume_fluxes,
							const RealType* RESTRICT cell_density,
							const RealType* RESTRICT cell_density_gradient,
							RealType* RESTRICT reconstructed_density_faces_y,
							RealType* RESTRICT mass_flux) {

  //y min only 
  const index_t iy = 0;
  const index_t iy_sym = ny;
  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect
       
    const index_t face_ooo = (nx * iy) + ix;
    const index_t face_sym = iy_sym * nx + ix;
      
    mass_flux[face_ooo] = 0.;
    mass_flux[face_sym] = 0.;
    reconstructed_density_faces_y[face_ooo] = 0.;
    reconstructed_density_faces_y[face_sym] = 0.;

  }
  
} //end  ReconstructMassFluxOrder2YCornerFluxesBoundaryWall



void ReconstructMassFluxOrder2CornerBoundaryWall(//in
						 index_t nx, 
						 index_t ny,
						 RealType dx,
						 RealType dy,
						 RealType dt,
						 const RealType* RESTRICT directional_lagrangian_density,
						 const RealType* RESTRICT reconstructed_density_faces_x,
						 const RealType* RESTRICT reconstructed_density_faces_y,
						 const RealType* RESTRICT volume_fluxes_corner, 
						 const RealType* RESTRICT in_vx, 
						 const RealType* RESTRICT in_vy,
						 //out
						 RealType* RESTRICT mass_corner_fluxes) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;
 
      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym] = 0.;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym] = 0.;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym_x] = 0.;
      mass_corner_fluxes[node_sym_y] = 0.;
      mass_corner_fluxes[node_sym_xy] = 0.; 
  }

} //end ReconstructMassFluxOrder2CornerBoundaryWall



void ReconstructMassFluxOrder2DiagCornerBoundaryWall(//in
						     index_t nx, 
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* RESTRICT predicted_u,
						     const RealType* RESTRICT predicted_v,
						     const RealType* RESTRICT directional_lagrangian_density,
						     const RealType* RESTRICT volume_fluxes_corner, 
						     const RealType* RESTRICT density_gradient_diag, 
						     const RealType* RESTRICT density_gradient_antidiag, 
						     //out
						     RealType* RESTRICT mass_corner_fluxes) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;
 
      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym] = 0.;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym] = 0.;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      mass_corner_fluxes[node_ooo] = 0.;
      mass_corner_fluxes[node_sym_x] = 0.;
      mass_corner_fluxes[node_sym_y] = 0.;
      mass_corner_fluxes[node_sym_xy] = 0.; 
  }

} //end ReconstructMassFluxOrder2DiagCornerBoundaryWall



void ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundaryWall(index_t nx, 
								     index_t ny, 
								     index_t halo_width,
								     const RealType dx,
								     const RealType dy,
								     const RealType dt,
								     const RealType* RESTRICT predicted_u,
								     const RealType* RESTRICT predicted_v,
								     const RealType* RESTRICT volume_fluxes,
								     const RealType* RESTRICT mass_flux,
								     const RealType* RESTRICT cell_variable,
								     const RealType* RESTRICT cell_variable_gradient,
								     RealType* RESTRICT reconstructed_variable_faces_x,
								     RealType* RESTRICT variable_flux) {
  

  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
    const index_t face_sym = face_ooo + nx;
      
    variable_flux[face_ooo] = 0.;
    variable_flux[face_sym] = 0.;
    reconstructed_variable_faces_x[face_ooo] = 0.; 
    reconstructed_variable_faces_x[face_sym] = 0.; 
    
  }

}// end ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundaryWall



void ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundaryWall(index_t nx, 
								     index_t ny, 
								     index_t halo_width,
								     const RealType dx,
								     const RealType dy,
								     const RealType dt,
								     const RealType* RESTRICT predicted_u,
								     const RealType* RESTRICT predicted_v,
								     const RealType* RESTRICT volume_fluxes,
								     const RealType* RESTRICT mass_flux,
								     const RealType* RESTRICT cell_variable,
								     const RealType* RESTRICT cell_variable_gradient,
								     RealType* RESTRICT reconstructed_variable_faces_y,
								     RealType* RESTRICT variable_flux) {


  index_t iy = 0;
  index_t iy_sym = ny;

  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect but boundary
      
    const index_t face_ooo = (iy * nx) + ix;
    const index_t face_sym = (iy_sym * nx) + ix;
      
    variable_flux[face_ooo] = 0.;
    variable_flux[face_sym] = 0.;
    reconstructed_variable_faces_y[face_ooo] = 0.; 
    reconstructed_variable_faces_y[face_sym] = 0.; 

  }
  
}// end ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundaryWall



void ReconstructIntensiveVariableFluxOrder2CornerBoundaryWall(//in
							      index_t nx, 
							      index_t ny,
							      RealType dx,
							      RealType dy,
							      RealType dt,
							      const RealType* RESTRICT cell_variable,
							      const RealType* RESTRICT reconstructed_variable_faces_x,
							      const RealType* RESTRICT reconstructed_variable_faces_y,
							      const RealType* RESTRICT mass_corner_fluxes, 
							      const RealType* RESTRICT in_vx, 
							      const RealType* RESTRICT in_vy,
							      //out
							      RealType* RESTRICT variable_flux_corner) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;
 
      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym] = 0.;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym] = 0.;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym_x] = 0.;
      variable_flux_corner[node_sym_y] = 0.;
      variable_flux_corner[node_sym_xy] = 0.; 
  }

} //end ReconstructIntensiveVariableFluxOrder2CornerBoundaryWall


void ReconstructIntensiveVariableFluxOrder2DiagCornerBoundaryWall(//in
								  index_t nx, 
								  index_t ny,
								  RealType dx,
								  RealType dy,
								  RealType dt,
								  const RealType* RESTRICT predicted_u,
								  const RealType* RESTRICT predicted_v,
								  const RealType* RESTRICT cell_variable,
								  const RealType* RESTRICT mass_corner_fluxes, 
								  const RealType* RESTRICT variable_gradient_diag, 
								  const RealType* RESTRICT variable_gradient_antidiag, 
								  //out
								  RealType* RESTRICT variable_flux_corner) {

  // X min
  for (index_t iy = 1; iy < ny; ++iy) {  
    index_t ix = 0;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = node_ooo + nx;
 
      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym] = 0.;
  
  }

  // Y min
  for (index_t ix = 1; ix < nx; ++ix) {  
    index_t iy = 0;
    index_t iy_sym = ny;

      const index_t node_ooo = (nx + 1)* iy + ix;
      const index_t node_sym = (nx + 1) * ny + ix;

      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym] = 0.;
  }

  // corners
  { index_t ix = 0; 
    index_t iy = 0;

      const index_t node_ooo = 0;
      const index_t node_sym_x = nx;
      const index_t node_sym_y = (nx + 1) * ny;
      const index_t node_sym_xy = (nx + 1) * (ny + 1) - 1;

      variable_flux_corner[node_ooo] = 0.;
      variable_flux_corner[node_sym_x] = 0.;
      variable_flux_corner[node_sym_y] = 0.;
      variable_flux_corner[node_sym_xy] = 0.; 
  }

} //end ReconstructIntensiveVariableFluxOrder2DiagCornerBoundaryWall

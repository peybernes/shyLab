// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_corner_fluxes_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


void ComputeDirectionalLagrangianQuantitiesCornerFluxesXY(index_t nx, 
							  index_t ny, 
							  RealType dt,
							  RealType dx,
							  RealType dy,
							  const RealType* RESTRICT in_vx,
							  const RealType* RESTRICT in_vy,
							  const RealType* RESTRICT cell_mass,
							  RealType* RESTRICT volume_fluxes,
							  RealType* RESTRICT volume_fluxes_y,
							  RealType* RESTRICT volume_fluxes_corner,
							  index_t* RESTRICT sign_x_corner_fluxes,
							  index_t* RESTRICT sign_y_corner_fluxes,
							  RealType* RESTRICT directional_lagrangian_volume,
							  RealType* RESTRICT directional_lagrangian_density) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    for (index_t ix = 0; ix < nx + 1; ++ix) {

      const index_t node_ooo = (nx + 1) * iy + ix;
      sign_x_corner_fluxes[node_ooo] = sign(in_vx[node_ooo]);
      sign_y_corner_fluxes[node_ooo] = sign(in_vy[node_ooo]);
      RealType vol_flux_corner = fabs(in_vx[node_ooo] * in_vy[node_ooo]) * dt * dt;
      if (vol_flux_corner / (dx * dy) < epsilon_corner) {
	vol_flux_corner = 0.0;
      }
      volume_fluxes_corner[node_ooo] = vol_flux_corner;

    }
  }

#pragma omp parallel 
  {
#pragma omp  for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianXpart1");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx + 1; ++ix) {
      
	const index_t face_ooo = ((nx + 1) * iy) + ix;
      
	const index_t prev_node = FaceXNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceXNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];                   
	const RealType vx_next = in_vx[next_node];                   
	const RealType vy_prev = in_vy[prev_node];                   
	const RealType vy_next = in_vy[next_node];

	RealType x = ix * dx;
	RealType ymin = iy * dy;
	RealType ymax = (iy + 1) * dy;

	RealType ymin_lag = ymin + vy_prev * dt;
	RealType ymax_lag = ymax + vy_next * dt;
	RealType xmin_lag = x + vx_prev * dt;
	RealType xmax_lag = x + vx_next * dt;

	// We compute x_ymin, x_ymax, ymin, ymax : the coordinates of the diagonal points of the volume flux

	RealType x_ymin, x_ymax;

	ymin = std::max(ymin, ymin_lag);
	ymax = std::min(ymax, ymax_lag);

	//Computation of coefficients of the line joining the two lagrangian vertices
	RealType a = ymax_lag - ymin_lag;
	RealType b = - (xmax_lag - xmin_lag);
	RealType c = - xmin_lag * a - ymin_lag * b;

	//Computation of x_ymin, x_ymax
	if (b == 0.0) {
	  x_ymin = xmin_lag;
	  x_ymax = xmin_lag;
	} else {
	  x_ymin = - (b * ymin + c) / a;
	  x_ymax = - (b * ymax + c) / a;
	}

	const RealType half = 0.5;
	RealType volume_flux = (ymax - ymin) * (half * (x_ymin + x_ymax) - x);
	if (fabs(volume_flux) / (dx * dy) < epsilon_face) {
	  volume_flux = 0.0;
	}

	volume_fluxes[face_ooo] = volume_flux;	

      }
      //likwid_markerStopRegion("directionnnal_lagrangianXpart1");
    }  
    
#pragma omp for
    for (index_t iy = 0; iy < ny + 1; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianYpart1");
       //#pragma unroll (UnrollFactor)
#pragma simd
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = (nx * iy) + ix;
      
	const index_t prev_node = FaceYNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceYNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];                   
	const RealType vx_next = in_vx[next_node];                   
	const RealType vy_prev = in_vy[prev_node];                   
	const RealType vy_next = in_vy[next_node];

	RealType y = iy * dy;
	RealType xmin = ix * dx;
	RealType xmax = (ix + 1) * dx;

	RealType xmin_lag = xmin + vx_prev * dt;
	RealType xmax_lag = xmax + vx_next * dt;
	RealType ymin_lag = y + vy_prev * dt;
	RealType ymax_lag = y + vy_next * dt;

	// We compute xmin, xmax, y_xmin, y_xmax : the coordinates of the diagonal points of the volume flux

	RealType y_xmin, y_xmax;

	xmin = std::max(xmin, xmin_lag);
	xmax = std::min(xmax, xmax_lag);

	//Computation of coefficients of the line joining the two lagrangian vertices
	RealType a = ymax_lag - ymin_lag;
	RealType b = - (xmax_lag - xmin_lag);
	RealType c = - xmin_lag * a - ymin_lag * b;

	//Computation of y_xmin, y_xmax
	if (a == 0.0) {
	  y_xmin = ymin_lag;
	  y_xmax = ymin_lag;
	} else {
	  y_xmin = - (a * xmin + c) / b;
	  y_xmax = - (a * xmax + c) / b;
	}

	const RealType half = 0.5;
	RealType volume_flux = (xmax - xmin) * (half * (y_xmin + y_xmax) - y);
	if (fabs(volume_flux) / (dx * dy) < epsilon_face) {
	  volume_flux = 0.0;
	}
      
	volume_fluxes_y[face_ooo] = volume_flux;     // 1 Store

     }
      //likwid_markerStopRegion("directionnnal_lagrangianYpart1");
    }  


#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianXpart2");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
	const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);

	const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

	const RealType volume_flux_prev_x = volume_fluxes[prev_face_x];    //1 Load
	const RealType volume_flux_next_x = volume_fluxes[next_face_x];     // 0 load (in cache?)
	const RealType volume_flux_prev_y = volume_fluxes_y[prev_face_y];    //1 Load
	const RealType volume_flux_next_y = volume_fluxes_y[next_face_y];
	const RealType current_mass = cell_mass[cell_ooo];

	//Compute algebraic corner fluxes
	RealType volume_flux_prev_prev = 0.0;
	RealType volume_flux_prev_next = 0.0;
	RealType volume_flux_next_prev = 0.0;
	RealType volume_flux_next_next = 0.0;

	if ((in_vx[node_m1m1] > 0.0) && (in_vy[node_m1m1] > 0.0)) {
	  volume_flux_prev_prev = - fabs(in_vx[node_m1m1]) * dt * fabs(in_vy[node_m1m1]) * dt;
	} else if ((in_vx[node_m1m1] < 0.0) && (in_vy[node_m1m1] < 0.0)) {
	  volume_flux_prev_prev = fabs(in_vx[node_m1m1]) * dt * fabs(in_vy[node_m1m1]) * dt;
	}
	if (fabs(volume_flux_prev_prev) / (dx * dy) < epsilon_corner) {
	  volume_flux_prev_prev = 0.0;
	}

	if ((in_vx[node_m1p1] > 0.0) && (in_vy[node_m1p1] < 0.0)) {
	  volume_flux_prev_next = - fabs(in_vx[node_m1p1]) * dt * fabs(in_vy[node_m1p1]) * dt;
	} else if ((in_vx[node_m1p1] < 0.0) && (in_vy[node_m1p1] > 0.0)) {
	  volume_flux_prev_next = fabs(in_vx[node_m1p1]) * dt * fabs(in_vy[node_m1p1]) * dt;
	}
	if (fabs(volume_flux_prev_next) / (dx * dy) < epsilon_corner) {
	  volume_flux_prev_next = 0.0;
	}

	if ((in_vx[node_p1m1] < 0.0) && (in_vy[node_p1m1] > 0.0)) {
	  volume_flux_next_prev = - fabs(in_vx[node_p1m1]) * dt * fabs(in_vy[node_p1m1]) * dt;
	} else if ((in_vx[node_p1m1] > 0.0) && (in_vy[node_p1m1] < 0.0)) {
	  volume_flux_next_prev = fabs(in_vx[node_p1m1]) * dt * fabs(in_vy[node_p1m1]) * dt;
	}
	if (fabs(volume_flux_next_prev) / (dx * dy) < epsilon_corner) {
	  volume_flux_next_prev = 0.0;
	}

	if ((in_vx[node_p1p1] < 0.0) && (in_vy[node_p1p1] < 0.0)) {
	  volume_flux_next_next = - fabs(in_vx[node_p1p1]) * dt * fabs(in_vy[node_p1p1]) * dt;
	} else if ((in_vx[node_p1p1] > 0.0) && (in_vy[node_p1p1] > 0.0)) {
	  volume_flux_next_next = fabs(in_vx[node_p1p1]) * dt * fabs(in_vy[node_p1p1]) * dt;
	}
	if (fabs(volume_flux_next_next) / (dx * dy) < epsilon_corner) {
	  volume_flux_next_next = 0.0;
	}

	const RealType current_directional_lagrangian_volume =         // 1 FMA   1 ADD
	  (dx * dy) - volume_flux_prev_x + volume_flux_next_x - volume_flux_prev_y + volume_flux_next_y
	  + volume_flux_prev_prev + volume_flux_prev_next + volume_flux_next_prev + volume_flux_next_next;

	const RealType current_directional_lagrangian_density =        // 1 DIV
	  current_mass / current_directional_lagrangian_volume;

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume; // 1 Store
      
	directional_lagrangian_density[cell_ooo] = current_directional_lagrangian_density; // 1 Store


	assert(0.0 < current_directional_lagrangian_volume);
	assert(0.0 < current_directional_lagrangian_density);
      
      }
      //likwid_markerStopRegion("directionnnal_lagrangianXpart2");
    }  
  }
} // total function manual optimist count --  3 Load + 3 Store    2 MUL 2 FMA 1 Add 1 DIV  



void ReconstructGradientDiag(index_t nx,
			     index_t ny,
			     RealType dx,
			     RealType dy,
			     RealType dt,
			     const RealType* RESTRICT predicted_u,
			     const RealType* RESTRICT predicted_v,
			     const RealType* RESTRICT lag_variable,
			     RealType* RESTRICT gradient_variable_diag,
			     RealType* RESTRICT gradient_variable_antidiag) {


#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
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
  }
 
}//end ReconstructGradientDiag



void ReconstructMassFluxOrder2XCornerFluxes(index_t nx, 
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
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("massFluxO2X");
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx; ++ix) {
      
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;
      const index_t node_om1 = FaceXNodeM1(face_ooo, iy, nx);
      const index_t node_op1 = FaceXNodeP1(face_ooo, iy, nx);

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell];         // 1 load
      const RealType next_cell_variable = cell_density[next_cell];         // 0 load
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell]; // 1 load
      const RealType next_cell_gradient = cell_density_gradient[next_cell]; // 0 load

      const RealType vol_flux = volume_fluxes[face_ooo];          // 1 load
      const RealType vol_flux_tot = volume_fluxes_tot[face_ooo];          // 1 load
      const RealType vol_flux_m1o = volume_fluxes_tot[face_m1o];      // 0 load
      const RealType vol_flux_p1o = volume_fluxes_tot[face_p1o];      // 0 load

      const RealType half = 0.5;
      
      const RealType y_om1 = (node_om1 / (nx + 1)) * dy;
      const RealType y_op1 = y_om1 + dy;
      const RealType y_om1_lag = y_om1 + predicted_v[node_om1] * dt;
      const RealType y_op1_lag = y_op1 + predicted_v[node_op1] * dt;
      const RealType dy_lag = std::min(y_op1,y_op1_lag) - std::max(y_om1,y_om1_lag);          
      const RealType one_over_dy_lag = 1.0 / dy_lag;                // 1 div

      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o *  one_over_dy_lag;   // 1 fma
      const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o * one_over_dy_lag;	 // 1 fma  1 add
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dx_lag_prev_corrected; // 1 mul 1 fma 
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dx_lag_next_corrected; // 1 mul 1 fma
      
      const RealType mass_flux_ooo =
	(half * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) +  // 2 abs 3 fma 1 mul
	(half * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo; // 1 store
      reconstructed_density_faces_x[face_ooo] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);

    }
    //likwid_markerStopRegion("massFluxO2X");
  }

} //end ReconstructMassFluxOrder2XCornerFluxes



void ReconstructMassFluxOrder2YCornerFluxes(index_t nx, 
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

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    //likwid_markerStartRegion("massFluxO2Y");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {
      
      
      const index_t face_ooo = (nx * iy) + ix;
      const index_t node_m1o = FaceYNodeM1(face_ooo, iy, nx);
      const index_t node_p1o = FaceYNodeP1(face_ooo, iy, nx);

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell]; // 1 load
      const RealType next_cell_variable = cell_density[next_cell];// 1 load
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell];// 1 load
      const RealType next_cell_gradient = cell_density_gradient[next_cell];// 1 load

      const RealType vol_flux = volume_fluxes[face_ooo];// 1 load
      const RealType vol_flux_tot = volume_fluxes_tot[face_ooo];// 1 load
      const RealType vol_flux_om1 = volume_fluxes_tot[face_om1];// 1 load
      const RealType vol_flux_op1 = volume_fluxes_tot[face_op1];// 1 load
      
      const RealType half = 0.5;

      const RealType x_m1o = (node_m1o - (node_m1o / (nx + 1)) * (nx + 1)) * dx;
      const RealType x_p1o = x_m1o + dx;
      const RealType x_m1o_lag = x_m1o + predicted_u[node_m1o] * dt;
      const RealType x_p1o_lag = x_p1o + predicted_u[node_p1o] * dt;
      const RealType dx_lag = std::min(x_p1o,x_p1o_lag) - std::max(x_m1o,x_m1o_lag); 
      const RealType one_over_dx_lag = 1.0 / dx_lag;                         

      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx_lag; // 1fma
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx_lag; // 1 add 1 fma	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected; // 1 mul  1 fma
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected; // 1 mul 1 fma
      
      const RealType mass_flux_ooo =
	(half * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + // 2 abs 3 fma 1 mul 
	(half * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo; // 1 store
      reconstructed_density_faces_y[face_ooo] = (half * (1 + sign(vol_flux_tot)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux_tot)) * next_cell_variable_o2);

    }
    //likwid_markerStopRegion("massFluxO2Y");
  }

} //end ReconstructMassFluxOrder2YCornerFluxes



void ReconstructMassFluxOrder2Corner(//in
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

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1)* iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);
      
      const index_t face_m1o = NodeFaceYM1OO(node_ooo,iy,nx);
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

      assert((dmass > 0.) || (dmass == 0.));

      mass_corner_fluxes[node_ooo] = dmass;  

    }
  }  

} //end ReconstructMassFluxOrder2Corner
			     


void ReconstructMassFluxOrder2DiagCorner(//in
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

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1) * iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
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

      }
    }
  }  

} //end ReconstructMassFluxOrder2DiagCorner



void ProjectMassCornerFluxes(//in
			     index_t nx, 
			     index_t ny,
			     RealType dt,
			     const RealType* RESTRICT mass_corner_fluxes,
			     const RealType* RESTRICT in_vx, 
			     const RealType* RESTRICT in_vy,
			     //out
			     RealType* RESTRICT out_cell_mass) {

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1)* iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
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
      RealType out_cell_mass_proj_giv = out_cell_mass[giving_cell] - dmass;
      RealType out_cell_mass_proj_rec = out_cell_mass[receiving_cell] + dmass;

      out_cell_mass[giving_cell] = out_cell_mass_proj_giv;
      out_cell_mass[receiving_cell] = out_cell_mass_proj_rec;    
    }
  }  

} //end ProjectMassCornerFluxes




void ReconstructIntensiveVariableFluxOrder2XCornerFluxes(index_t nx, 
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
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("fluxO2X");
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx; ++ix) {
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;
      const index_t node_om1 = FaceXNodeM1(face_ooo, iy, nx);
      const index_t node_op1 = FaceXNodeP1(face_ooo, iy, nx);
      
      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];        // 1 load
      const RealType next_cell_variable = cell_variable[next_cell];        // 0 load
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell]; // 1 load
      const RealType next_cell_gradient = cell_variable_gradient[next_cell]; // 0 load

      const RealType mass_flux_face = mass_flux[face_ooo];              // 1 load
      
      const RealType vol_flux = volume_fluxes[face_ooo];            // 1 load
      const RealType vol_flux_m1o = volume_fluxes[face_m1o];            // 1 load
      const RealType vol_flux_p1o = volume_fluxes[face_p1o];            // 0 load
      
      const RealType half = 0.5;

      const RealType y_om1 = (node_om1 / (nx + 1)) * dy;
      const RealType y_op1 = y_om1 + dy;
      const RealType y_om1_lag = y_om1 + predicted_v[node_om1] * dt;
      const RealType y_op1_lag = y_op1 + predicted_v[node_op1] * dt;
      const RealType dy_lag = std::min(y_op1,y_op1_lag) - std::max(y_om1,y_om1_lag);          
      const RealType one_over_dy_lag = 1.0 / dy_lag;                // 1 div

      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o *  one_over_dy_lag;   // 1 fma
      const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o * one_over_dy_lag;	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dx_lag_prev_corrected; // 1 fma 1 mul
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dx_lag_next_corrected; // 1 fma 1 mul
      
      const RealType variable_flux_ooo =
	(half * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) +  // 2 abs  3 fma 1 mul
	(half * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
     
      variable_flux[face_ooo] = variable_flux_ooo; // 1 store
      reconstructed_variable_faces_x[face_ooo] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2);

    }
    //likwid_markerStopRegion("fluxO2X");
 }

} //end ReconstructIntensiveVariableFluxOrder2XCornerFluxes



void ReconstructIntensiveVariableFluxOrder2YCornerFluxes(index_t nx, 
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

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    //likwid_markerStartRegion("FluxO2Y");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {
      
      const index_t face_ooo = (nx * iy) + ix;
      const index_t node_m1o = FaceYNodeM1(face_ooo, iy, nx);
      const index_t node_p1o = FaceYNodeP1(face_ooo, iy, nx);

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_variable[prev_cell];  // 1 load
      const RealType next_cell_variable = cell_variable[next_cell]; // 1 load
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];// 1 load
      const RealType next_cell_gradient = cell_variable_gradient[next_cell];// 1 load
 
      const RealType mass_flux_face = mass_flux[face_ooo];// 1 load
      
      const RealType vol_flux = volume_fluxes[face_ooo];            // 1 load
      const RealType vol_flux_om1 = volume_fluxes[face_om1];// 1 load
      const RealType vol_flux_op1 = volume_fluxes[face_op1];// 1 load
      
      const RealType half = 0.5;

      const RealType x_m1o = (node_m1o - (node_m1o / (nx + 1)) * (nx + 1)) * dx;
      const RealType x_p1o = x_m1o + dx;
      const RealType x_m1o_lag = x_m1o + predicted_u[node_m1o] * dt;
      const RealType x_p1o_lag = x_p1o + predicted_u[node_p1o] * dt;
      const RealType dx_lag = std::min(x_p1o,x_p1o_lag) - std::max(x_m1o,x_m1o_lag); 
      const RealType one_over_dx_lag = 1.0 / dx_lag;                         

      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx_lag; // 1fma
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx_lag;
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dy_lag_prev_corrected; // 1 mul  1 fma
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dy_lag_next_corrected; // 1 mul 1 fma
      
      const RealType variable_flux_ooo =
	(half * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) +  // 2 abs 3 fma 1 mul
	(half * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
      variable_flux[face_ooo] = variable_flux_ooo; // 1 store
      reconstructed_variable_faces_y[face_ooo] = (half * (1 + sign(vol_flux)) * prev_cell_variable_o2) + (half * (1 - sign(vol_flux)) * next_cell_variable_o2);  

    }
    //likwid_markerStopRegion("FluxO2Y");
  }
} // end ReconstructIntensiveVariableFluxOrder2YCornerFluxes



void ReconstructIntensiveVariableFluxOrder2Corner(//in
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

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1)* iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
      const index_t cell_p1m1 = NodeCellP1M1(node_ooo,iy,nx);
      const index_t cell_p1p1 = NodeCellP1P1(node_ooo,iy,nx);

      const index_t face_m1o = NodeFaceYM1OO(node_ooo,iy,nx);
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
    }
  }  

} //end ReconstructIntensiveVariableFluxOrder2Corner



void ReconstructIntensiveVariableFluxOrder2DiagCorner(//in
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

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1) * iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
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

      }
    }
  }  

} //end ReconstructintensiveVariableFluxOrder2DiagCorner


void MassProjectIntensiveVariableCornerFluxes(index_t nx, 
					      index_t ny,
					      RealType dt,
					      const RealType* RESTRICT in_cell_variable,
					      const RealType* RESTRICT variable_flux_corner,
					      const RealType* RESTRICT in_vx, 
					      const RealType* RESTRICT in_vy, 
					      const RealType* RESTRICT out_cell_mass,
					      RealType* RESTRICT out_cell_variable) {

  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {  
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1)* iy + ix;

      const index_t cell_m1m1 = NodeCellM1M1(node_ooo,iy,nx);
      const index_t cell_m1p1 = NodeCellM1P1(node_ooo,iy,nx);
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
  }  

} //end MassProjectIntensiveVariableCornerFluxes


void CheckMassPositivity(index_t nx,
			 index_t ny,
			 RealType* out_cell_mass) {
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      index_t cell_ooo = iy * nx + ix;
      RealType out_mass = out_cell_mass[cell_ooo];

      if (fabs(out_mass) < 1.0e-14) {
	out_cell_mass[cell_ooo] = 0.0;
      }

    }
  }
}



void ProjectNodalIntensiveVariableCornerFluxes(index_t nx, 
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

 
   //Boucle sur les noeuds
#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) { 

    const index_t node_oooo = ((nx + 1) * iy) + ix;

    const index_t node_m1oo = NodeNodeM1O(node_oooo, iy, nx);
    const index_t node_p1oo = NodeNodeP1O(node_oooo, iy, nx);
    const index_t node_oom1 = NodeNodeOM1(node_oooo, iy, nx);
    const index_t node_oop1 = NodeNodeOP1(node_oooo, iy, nx);
    const index_t node_m1m1 = NodeNodeM1M1(node_oooo, iy, nx);
    const index_t node_p1m1 = NodeNodeP1M1(node_oooo, iy, nx);
    const index_t node_m1p1 = NodeNodeM1P1(node_oooo, iy, nx);
    const index_t node_p1p1 = NodeNodeP1P1(node_oooo, iy, nx);
      
    const index_t cell_m1m1 = NodeCellM1M1(node_oooo, iy, nx);
    const index_t cell_p1m1 = NodeCellP1M1(node_oooo, iy, nx);
    const index_t cell_m1p1 = NodeCellM1P1(node_oooo, iy, nx);
    const index_t cell_p1p1 = NodeCellP1P1(node_oooo, iy, nx);

#include "direct_projection_corner_fluxes_nodal_2d_computation.h"

    out_variable[node_oooo] = out_variable[node_oooo] + total_variable_flux;
  }

}

} //end ProjectNodalIntensiveVariableCornerFluxes

#include "kernel_lagrange_2d.h" 
#include "cartesian_connectivity.h"
#include "kernel_tools.h"

#include <cassert>
#include <cstdio>
#include <algorithm>




RealType TimeStep(int nx,
		  int ny,
		  const RealType dx,
		  const RealType dy,
		  const RealType CFL,
		  const RealType* RESTRICT density,
		  const RealType* RESTRICT pressure,
		  const RealType* RESTRICT in_velocity_x,
		  const RealType* RESTRICT in_velocity_y) {
  RealType max_velocity = 0.0;
  const RealType GAMMA = 1.4;
  
  #pragma omp parallel 
  {
    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny + 1; ++iy) {
      for (int ix = 0; ix < nx + 1; ++ix) {
	max_velocity = std::max(max_velocity,in_velocity_x[iy * (nx + 1) + ix]);
      }
    }


    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny + 1; ++iy) {
      for (int ix = 0; ix < nx + 1; ++ix) {
	max_velocity =  std::max(max_velocity,in_velocity_y[iy * (nx + 1) + ix]);
      }
    }

    
    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
	const RealType	p_ooo = pressure[iy * nx + ix];
	const RealType rho_ooo =  density[iy * nx + ix];
	const RealType speed_of_sound = std::sqrt(GAMMA*p_ooo/rho_ooo);
	max_velocity =  std::max(max_velocity,speed_of_sound);
      }
    }
  }

  RealType dt = CFL * std::min(dx,dy)/max_velocity;
  return dt;
}

void LagrangePressurePredicted(int nx,
			       int ny,
			       RealType dt,
		 	       RealType dx,
		       	       RealType dy,
	       		       const RealType* RESTRICT in_mass,
			       const RealType* RESTRICT in_energy,	 
			       const RealType* RESTRICT in_velocity_x,
		       	       const RealType* RESTRICT in_velocity_y,
			       RealType* RESTRICT out_pressure,
	       	      	       RealType* RESTRICT out_predicted_pressure,
       			       RealType* RESTRICT out_pseudo_pressure) {

  const double gamma = 1.4;

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {  // sqrt fction vectorization with icpc but not with gcc

      //DATA LOAD
      const int cell_ooo = nx * iy + ix; 
   
      const int node_SW = CellNodeM1M1(cell_ooo, iy, nx);
      const int node_SE = CellNodeP1M1(cell_ooo, iy, nx);
      const int node_NW = CellNodeM1P1(cell_ooo, iy, nx);
      const int node_NE = CellNodeP1P1(cell_ooo, iy, nx);

      const RealType u_x_sw = in_velocity_x[node_SW];
      const RealType u_x_se = in_velocity_x[node_SE];
      const RealType u_x_nw = in_velocity_x[node_NW];
      const RealType u_x_ne = in_velocity_x[node_NE];

      const RealType u_y_sw = in_velocity_y[node_SW];
      const RealType u_y_se = in_velocity_y[node_SE];
      const RealType u_y_nw = in_velocity_y[node_NW];
      const RealType u_y_ne = in_velocity_y[node_NE];

      const RealType mass_ooo = in_mass[cell_ooo];

      const RealType e_ooo = in_energy[cell_ooo];
      const RealType rho_ooo = mass_ooo / (dx * dy);
      const RealType p_ooo = EquationOfState( rho_ooo, e_ooo);
                                              
      const RealType delta_vol = 0.5 * dt * 0.5 * (u_x_se + u_x_ne - u_x_sw - u_x_nw) * dy + 
	0.5 *  dt * 0.5 * (u_y_nw + u_y_ne - u_y_sw - u_y_se) * dx;

      const RealType delta_v =  0.5 * (u_x_se + u_x_ne - u_x_sw - u_x_nw)
      	+ 0.5 * (u_y_nw + u_y_ne - u_y_sw - u_y_se);

      // for perfect gas law
      const RealType cs2 = gamma * p_ooo / rho_ooo;
      const RealType cs  = std::sqrt(cs2);  
      
      // viscous pressure  -- one possible formulation with both linear 
      // and quadratic coefficient equals to 1.0 ; 
      // using negative part of delta_v --  -0.5*(delta_v-abs(delta_v)) for vectorisation*/
      const RealType delta_v_neg = -0.5 * (delta_v - fabs(delta_v));
      const RealType q_ooo   = 1.0 * rho_ooo * cs * delta_v_neg + 
	1.0 * rho_ooo * delta_v_neg * delta_v_neg; 
      
      const RealType div_u_ooo = 
	1.0 / dx * 0.5 * (u_x_se + u_x_ne - u_x_sw - u_x_nw) +
	1.0 / dy * 0.5 * (u_y_nw + u_y_ne - u_y_sw - u_y_se);
      
      const RealType e_lag_ooo = e_ooo 
	- 0.5 * dt * (p_ooo + q_ooo) * div_u_ooo / mass_ooo * dx * dy;


      const RealType out_predicted_p_ooo = EquationOfState(mass_ooo / (dx * dy + delta_vol), e_lag_ooo);
  

      out_pressure[cell_ooo]       = p_ooo ;
      out_predicted_pressure[cell_ooo] = out_predicted_p_ooo;
      out_pseudo_pressure[cell_ooo] = q_ooo;

    }
  }
}


void LagrangeVelocityPredicted(int nx,
			       int ny,
			       RealType dt,
			       RealType dx,
			       RealType dy,
			       const RealType* RESTRICT in_mass,
			       const RealType* RESTRICT in_pressure,
			       const RealType* RESTRICT in_pseudo_pressure,
			       const RealType* RESTRICT in_velocity_x,
			       const RealType* RESTRICT in_velocity_y,
			       RealType* RESTRICT out_velocity_x,
			       RealType* RESTRICT out_velocity_y)
{

#pragma omp parallel for
  for (int iy = 1; iy < ny; ++iy) {
    for (int ix = 1; ix < nx; ++ix) {

      const int node_ooo = (nx + 1) * iy + ix; 

      const int cell_SW = NodeCellM1M1(node_ooo, iy, nx);
      const int cell_SE = NodeCellP1M1(node_ooo, iy, nx);
      const int cell_NW = NodeCellM1P1(node_ooo, iy, nx);
      const int cell_NE = NodeCellP1P1(node_ooo, iy, nx);

#include "kernel_lagrange_velocity.h"

      out_velocity_x[node_ooo] = out_u_x;
      out_velocity_y[node_ooo] = out_u_y;
        
    }
  }
}



void LagrangeCorrection(int nx,
			int ny,
			RealType dt,
			RealType dx,
			RealType dy,
			const RealType* RESTRICT in_mass,
			const RealType* RESTRICT in_energy,
			const RealType* RESTRICT in_pressure,
			const RealType* RESTRICT in_pseudo_pressure,
			const RealType* RESTRICT in_velocity_x,
			const RealType* RESTRICT in_velocity_y,
			RealType* RESTRICT out_energy) {

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx ; ++ix) {

      //DATA LOAD
      const int cell_ooo = nx * iy + ix;

      const int node_SW = CellNodeM1M1(cell_ooo,iy,nx);
      const int node_SE = CellNodeP1M1(cell_ooo,iy,nx);
      const int node_NW = CellNodeM1P1(cell_ooo,iy,nx);
      const int node_NE = CellNodeP1P1(cell_ooo,iy,nx);
      
      const RealType u_x_sw = in_velocity_x[node_SW];
      const RealType u_x_se = in_velocity_x[node_SE];
      const RealType u_x_nw = in_velocity_x[node_NW];
      const RealType u_x_ne = in_velocity_x[node_NE];

      const RealType u_y_sw = in_velocity_y[node_SW];
      const RealType u_y_se = in_velocity_y[node_SE];
      const RealType u_y_nw = in_velocity_y[node_NW];
      const RealType u_y_ne = in_velocity_y[node_NE];

      const RealType mass_ooo = in_mass[cell_ooo];

      const RealType e_ooo = in_energy[cell_ooo];
      const RealType p_ooo = in_pressure[cell_ooo];
      const RealType q_ooo = in_pseudo_pressure[cell_ooo];
     
      const RealType div_u_ooo = 
	1.0 / dx * (0.5 * (u_x_se + u_x_ne - u_x_sw - u_x_nw)) +
	1.0 / dy * (0.5 * (u_y_nw + u_y_ne - u_y_sw - u_y_se));
      
      const RealType e_lag_ooo = e_ooo 
	- dt * (p_ooo + q_ooo) * div_u_ooo / mass_ooo * dx * dy;
      

      out_energy[cell_ooo] = e_lag_ooo;

    }
  }
}


void LagrangeVelocityCorrection(int nx,
				int ny,
				const RealType* RESTRICT in_velocity_x,
				const RealType* RESTRICT in_velocity_y,
				const RealType* RESTRICT predicted_velocity_x,
				const RealType* RESTRICT predicted_velocity_y,
				RealType* RESTRICT lagrangian_velocity_x,
				RealType* RESTRICT lagrangian_velocity_y)
{
#pragma omp parallel for
  for (int iy = 0; iy < ny+1; ++iy) {
    for (int ix = 0; ix < nx+1 ; ++ix) {
      const int node_ooo = iy * (nx + 1) + ix;
      RealType in_u = in_velocity_x[node_ooo];
      RealType in_v = in_velocity_y[node_ooo];
      RealType u_predicted = predicted_velocity_x[node_ooo];
      RealType v_predicted = predicted_velocity_y[node_ooo];
      RealType u_lag = 2.0 * u_predicted - in_u;
      RealType v_lag = 2.0 * v_predicted - in_v;

      lagrangian_velocity_x[node_ooo] = u_lag;
      lagrangian_velocity_y[node_ooo] = v_lag;
    }
  }
}


/*============================================*/
/* Specific functions for boundary conditions */ 
/*============================================*/

void PeriodicBoundaryCopy(int nx,
			  int ny,
			  RealType* RESTRICT in_velocity_x,
			  RealType* RESTRICT in_velocity_y)
{
  for (int ix = 0; ix < nx + 1; ++ix){  // not vectorized 
    
    const int pos_bottom_border = ix;
    const int pos_top_border    = ny * (nx + 1) + ix;
    in_velocity_x[pos_top_border] =  in_velocity_x[pos_bottom_border];
    in_velocity_y[pos_top_border] =  in_velocity_y[pos_bottom_border];
  }				   

  for (int iy = 0; iy < ny + 1; ++iy){ // not vectorized 

    const int pos_lefft_border = iy * (nx + 1);
    const int pos_right_border = (iy + 1) * (nx + 1) - 1;
    in_velocity_x[pos_right_border] =  in_velocity_x[pos_lefft_border];
    in_velocity_y[pos_right_border] =  in_velocity_y[pos_lefft_border];
  }
}



void PeriodicBoundaryVelocityPrediction(int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					const RealType* RESTRICT in_mass,
					const RealType* RESTRICT in_pressure,
					const RealType* RESTRICT in_pseudo_pressure,
					const RealType* RESTRICT in_velocity_x,
					const RealType* RESTRICT in_velocity_y,
					RealType* RESTRICT out_velocity_x,
					RealType* RESTRICT out_velocity_y)
{

  for (int ix = 1 ; ix < nx; ++ix){ // not verctorized
    //bottom -- top computation
    const int pos_bottom_border = ix;
    const int pos_top_border    = ny * (nx + 1) + ix;
    const int node_ooo = pos_bottom_border;

    const int cell_SW  = NodeCellM1M1(pos_top_border, ny, nx);
    const int cell_SE  = NodeCellP1M1(pos_top_border, ny, nx);
    const int cell_NW  = NodeCellM1P1(pos_bottom_border, 0, nx);
    const int cell_NE  = NodeCellP1P1(pos_bottom_border, 0, nx);

#include "kernel_lagrange_velocity.h"

    out_velocity_x[pos_bottom_border] = out_u_x; 
    out_velocity_x[pos_top_border] = out_u_x;
    out_velocity_y[pos_bottom_border] = out_u_y;
    out_velocity_y[pos_top_border] = out_u_y;
  
  }

  for (int iy = 1; iy < ny; ++iy){  // not verctorized
    //left -- right computation 
    const int pos_lefft_border = iy * (nx + 1);
    const int pos_right_border = (iy + 1) * (nx + 1) - 1;
    const int node_ooo = pos_lefft_border;

    const int cell_SW  = NodeCellM1M1(pos_right_border, iy , nx);
    const int cell_SE  = NodeCellP1M1(pos_lefft_border, iy , nx);
    const int cell_NW  = NodeCellM1P1(pos_right_border, iy, nx);
    const int cell_NE  = NodeCellP1P1(pos_lefft_border, iy, nx);
 
#include "kernel_lagrange_velocity.h"

    out_velocity_x [pos_lefft_border] = out_u_x; 
    out_velocity_x [pos_right_border] = out_u_x;
    out_velocity_y [pos_lefft_border] = out_u_y;
    out_velocity_y [pos_right_border] = out_u_y;
 }

 // corners
 { 
   const int node_ooo = 0;
 
   const int cell_SW  = ny * nx -1;
   const int cell_SE  = (ny - 1) * nx ;
   const int cell_NW  = nx - 1;
   const int cell_NE  = 0;

#include "kernel_lagrange_velocity.h"

   out_velocity_x[0]  = out_u_x; 
   out_velocity_y[0]  = out_u_y;
   out_velocity_x[nx]  = out_u_x; 
   out_velocity_y[nx]  = out_u_y;
   out_velocity_x[ny * (nx + 1)]  = out_u_x; 
   out_velocity_y[ny * (ny + 1)]  = out_u_y;
   out_velocity_x[(ny + 1) * (nx + 1) - 1]  = out_u_x; 
   out_velocity_y[(ny + 1) * (ny + 1) - 1 ]  = out_u_y; 
 }
 

}

void RtBoundaryVelocityPrediction(int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					const RealType* RESTRICT in_mass,
					const RealType* RESTRICT in_pressure,
					const RealType* RESTRICT in_pseudo_pressure,
					const RealType* RESTRICT in_velocity_x,
					const RealType* RESTRICT in_velocity_y,
					RealType* RESTRICT out_velocity_x,
					RealType* RESTRICT out_velocity_y)
{

  for (int ix = 1 ; ix < nx; ++ix){ // not verctorized
    //bottom -- top computation
    const int pos_bottom_border = ix;
    const int pos_top_border    = ny * (nx + 1) + ix;

    out_velocity_x[pos_bottom_border] = 0.; 
    out_velocity_x[pos_top_border] = 0.;
    out_velocity_y[pos_bottom_border] = 0.;
    out_velocity_y[pos_top_border] = 0.;
  
  }

  for (int iy = 1; iy < ny; ++iy){  // not verctorized
    //left -- right computation 
    const int pos_lefft_border = iy * (nx + 1);
    const int pos_right_border = (iy + 1) * (nx + 1) - 1;
    const int node_ooo = pos_lefft_border;

    const int cell_SW  = NodeCellM1M1(pos_right_border, iy , nx);
    const int cell_SE  = NodeCellP1M1(pos_lefft_border, iy , nx);
    const int cell_NW  = NodeCellM1P1(pos_right_border, iy, nx);
    const int cell_NE  = NodeCellP1P1(pos_lefft_border, iy, nx);
 
#include "kernel_lagrange_velocity.h"

    out_velocity_x [pos_lefft_border] = out_u_x; 
    out_velocity_x [pos_right_border] = out_u_x;
    out_velocity_y [pos_lefft_border] = out_u_y;
    out_velocity_y [pos_right_border] = out_u_y;
 }

 // corners
 { 

   out_velocity_x[0]  = 0.; 
   out_velocity_y[0]  = 0.;
   out_velocity_x[nx]  = 0.; 
   out_velocity_y[nx]  = 0.;
   out_velocity_x[ny * (nx + 1)]  = 0.; 
   out_velocity_y[ny * (ny + 1)]  = 0.;
   out_velocity_x[(ny + 1) * (nx + 1) - 1]  = 0.; 
   out_velocity_y[(ny + 1) * (ny + 1) - 1 ]  = 0.; 
 }
 

}

void CopyVelocity(int nx,
		  int ny,
		  const RealType* RESTRICT in_u,
		  RealType* RESTRICT out_u) {
  for (int iy = 0; iy < ny + 1; ++iy) { //not vectorised (seems unefficent)
    for (int ix = 0; ix < nx + 1; ++ix) {
      const int node_o = iy * (nx + 1) + ix;
      out_u[node_o] = in_u[node_o];
    }
  }
}

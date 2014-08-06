// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_2d_boundary_conditions.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"




  // Periodic conditions


void ProjectNodalIntensiveVariablePeriodicBoundaryDirect(index_t nx, 
				                 index_t ny, 
				                 index_t halo_width,
				                 const RealType* RESTRICT lag_cell_mass,
				                 const RealType* RESTRICT out_cell_mass,
				                 const RealType* RESTRICT in_vx,
				                 const RealType* RESTRICT mass_flux_x,
				                 const RealType* RESTRICT mass_flux_y,
				                 RealType* RESTRICT out_vx) {

  // xmin (xmax is deduced in periodic conditions)

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t ix = 0; ix < halo_width; ++ix) {
    for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t node_m1o = NodeNodeM1O(node_sym, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_sym, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t facexm1m1 = NodeFaceXM1M1(node_sym, iy, nx);
      const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
      const index_t facexm1p1 = NodeFaceXM1P1(node_sym, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
      const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
      const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_sym, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_sym, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "direct_nodal_projection_2d_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
        const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
	assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux_x = 0.25 * (mass_flux_x_m1m1 + mass_flux_x_oom1 + mass_flux_x_m1p1 + mass_flux_x_oop1);
      
      const RealType next_dual_mass_flux_x = 0.25 * (mass_flux_x_p1m1 + mass_flux_x_oom1 + mass_flux_x_p1p1 + mass_flux_x_oop1); 
 
      const RealType prev_dual_mass_flux_y = 0.25 * (mass_flux_y_m1m1 + mass_flux_y_m1oo + mass_flux_y_p1m1 + mass_flux_y_p1oo);
      
      const RealType next_dual_mass_flux_y = 0.25 * (mass_flux_y_m1p1 + mass_flux_y_m1oo + mass_flux_y_p1p1 + mass_flux_y_p1oo);
      
#include "reconstruct_dual_variable_xy_computation.h"   
   
      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;

      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 

      out_vx[node_sym] = out_moment_x_ooo / out_node_mass_ooo; 
 
    }
  }  

  // ymin (ymax is deduced in periodic conditions)

  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + (nx + 1) * ny;
      const int iy_sym = ny; 
	
      index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      index_t node_om1 = NodeNodeOM1(node_sym, iy, nx);
      index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);    
      
      index_t cellm1m1 = NodeCellM1M1(node_sym, iy_sym, nx);
      index_t cellp1m1 = NodeCellP1M1(node_sym, iy_sym, nx);
      index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

 
      const index_t facexm1m1 = NodeFaceXM1M1(node_sym, iy_sym, nx);
      const index_t facexoom1 = NodeFaceXOOM1(node_sym, iy_sym, nx);
      const index_t facexm1p1 = NodeFaceXM1P1(node_ooo, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
      const index_t facexp1m1 = NodeFaceXP1M1(node_sym, iy_sym, nx);
      const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy_sym, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_sym, iy_sym, nx);

#include "direct_nodal_projection_2d_data_load.h"

       const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
        const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
	assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux_x = 0.25 * (mass_flux_x_m1m1 + mass_flux_x_oom1 + mass_flux_x_m1p1 + mass_flux_x_oop1);
      
      const RealType next_dual_mass_flux_x = 0.25 * (mass_flux_x_p1m1 + mass_flux_x_oom1 + mass_flux_x_p1p1 + mass_flux_x_oop1); 
 
      const RealType prev_dual_mass_flux_y = 0.25 * (mass_flux_y_m1m1 + mass_flux_y_m1oo + mass_flux_y_p1m1 + mass_flux_y_p1oo);
      
      const RealType next_dual_mass_flux_y = 0.25 * (mass_flux_y_m1p1 + mass_flux_y_m1oo + mass_flux_y_p1p1 + mass_flux_y_p1oo);

#include "reconstruct_dual_variable_xy_computation.h"   
   
      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;


      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 

      out_vx[node_sym] = out_moment_x_ooo / out_node_mass_ooo; 

    }

  }  

  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    const index_t node_ooo = node_cmm;
	
    index_t node_m1o = nx - 1;
    index_t node_p1o = 1;
    index_t node_om1 = (ny - 1) * (nx + 1);
    index_t node_op1 = nx + 1;    
      
    index_t cellm1m1 = nx * ny - 1;
    index_t cellp1m1 = (ny - 1) * nx;
    index_t cellm1p1 = nx - 1;
    index_t cellp1p1 = 0 ;

 
    const index_t facexm1m1 = ny * (nx + 1) - 2;
    const index_t facexoom1 = (ny - 1) * (nx + 1);
    const index_t facexm1p1 = nx - 1;
    const index_t facexoop1 = 0;
    const index_t facexp1p1 = 1;
    const index_t facexp1m1 = (ny - 1) * (nx + 1) + 1 ;
    const index_t faceym1m1 = ny * nx - 1;
    const index_t faceym1oo = nx - 1;
    const index_t faceym1p1 = nx + nx - 1;
    const index_t faceyp1p1 = nx;
    const index_t faceyp1oo = 0;
    const index_t faceyp1m1 = (ny - 1) * nx;

#include "direct_nodal_projection_2d_data_load.h"

    const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
    assert(0.0 < lag_node_mass_ooo);
    const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
    assert(0.0 < out_node_mass_ooo);

    const RealType prev_dual_mass_flux_x = 0.25 * (mass_flux_x_m1m1 + mass_flux_x_oom1 + mass_flux_x_m1p1 + mass_flux_x_oop1);
      
    const RealType next_dual_mass_flux_x = 0.25 * (mass_flux_x_p1m1 + mass_flux_x_oom1 + mass_flux_x_p1p1 + mass_flux_x_oop1); 
 
    const RealType prev_dual_mass_flux_y = 0.25 * (mass_flux_y_m1m1 + mass_flux_y_m1oo + mass_flux_y_p1m1 + mass_flux_y_p1oo);
      
    const RealType next_dual_mass_flux_y = 0.25 * (mass_flux_y_m1p1 + mass_flux_y_m1oo + mass_flux_y_p1p1 + mass_flux_y_p1oo);

#include "reconstruct_dual_variable_xy_computation.h"   
   
    const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

    const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;


    const RealType out_v_x = out_moment_x_ooo / out_node_mass_ooo; 

    out_vx[node_cmm] = out_v_x;
    out_vx[node_cpm] = out_v_x;
    out_vx[node_cmp] = out_v_x;
    out_vx[node_cpp] = out_v_x;

  }
} // end ProjectNodalIntensiveVariablePeriodicBoundaryDirect



void ProjectNodalIntensiveVariableOrder2PeriodicBoundaryDirect(index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable_x,
						  const RealType* gradient_variable_y,
						  const RealType* mass_flux_x,
						  const RealType* mass_flux_y,
						  RealType* out_variable) {

 // x min and x max
   
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect but boundary
    index_t ix = 0;

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + nx;

    const index_t node_m1o = NodeNodeM1O(node_sym, iy, nx);
    const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
    const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
    const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
    const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
    const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
    const index_t cellm1p1 = NodeCellM1P1(node_sym, iy, nx);
    const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

    const index_t facexm1m1 = NodeFaceXM1M1(node_sym, iy, nx);
    const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
    const index_t facexm1p1 = NodeFaceXM1P1(node_sym, iy, nx);
    const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
    const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
    const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
    const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
    const index_t faceym1oo = NodeFaceYM1OO(node_sym, iy, nx);
    const index_t faceym1p1 = NodeFaceYM1P1(node_sym, iy, nx);
    const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
    const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
    const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "direct_projection_nodal_o2_2d_computation.h"

    out_variable[node_ooo] = out_nodal_variable ;
    out_variable[node_sym] = out_nodal_variable;
  }


  //ymin and y max

  index_t iy = 0;
  index_t iy_sym = ny;
  for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + ny * (nx + 1);

    const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
    const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
    const index_t node_om1 = NodeNodeOM1(node_sym, iy_sym, nx);
    const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
     
    const index_t cellm1m1 = NodeCellM1M1(node_sym, iy_sym, nx);
    const index_t cellp1m1 = NodeCellP1M1(node_sym, iy_sym, nx);
    const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
    const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

    const index_t facexm1m1 = NodeFaceXM1M1(node_sym, iy_sym, nx);
    const index_t facexoom1 = NodeFaceXOOM1(node_sym, iy_sym, nx);
    const index_t facexm1p1 = NodeFaceXM1P1(node_ooo, iy, nx);
    const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
    const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
    const index_t facexp1m1 = NodeFaceXP1M1(node_sym, iy_sym, nx);
    const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy_sym, nx);
    const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
    const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
    const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
    const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
    const index_t faceyp1m1 = NodeFaceYP1M1(node_sym, iy_sym, nx);

#include "direct_projection_nodal_o2_2d_computation.h"
  
   
    out_variable[node_ooo] = out_nodal_variable ;
    out_variable[node_sym] = out_nodal_variable;
  } // end Y boundary

  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;

    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    const index_t node_ooo = node_cmm;
	
    const index_t node_m1o = nx - 1;
    const index_t node_p1o = 1;
    const index_t node_om1 = (ny - 1) * (nx + 1) ;
    const index_t node_op1 = nx + 1;
 
    const index_t cellm1m1 = nx * ny - 1;
    const index_t cellp1m1 = (ny - 1) * nx;
    const index_t cellm1p1 = nx - 1;
    const index_t cellp1p1 = 0;
     
 
    const index_t facexm1m1 = ny * (nx + 1) - 2;
    const index_t facexoom1 = (ny - 1) * (nx + 1);
    const index_t facexm1p1 = nx - 1;
    const index_t facexoop1 = 0;
    const index_t facexp1p1 = 1;
    const index_t facexp1m1 = (ny - 1) * (nx + 1) + 1 ;
    const index_t faceym1m1 = ny * nx - 1;
    const index_t faceym1oo = nx - 1;
    const index_t faceym1p1 = 2 * nx - 1;
    const index_t faceyp1oo = 0;
    const index_t faceyp1p1 = nx;
    const index_t faceyp1m1 = (ny - 1) * nx;

#include "direct_projection_nodal_o2_2d_computation.h" 

    out_variable[node_cmm] = out_nodal_variable ;
    out_variable[node_cmp] = out_nodal_variable ;
    out_variable[node_cpm] = out_nodal_variable ;
    out_variable[node_cpp] = out_nodal_variable ;
  } // end corners

} // end ProjectNodalIntensiveVariableOrder2PeriodicBoundaryDirect





  //Rayleigh Taylor conditions


void ProjectNodalIntensiveVariableOrder2BoundaryRtDirect(index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable_x,
						  const RealType* gradient_variable_y,
						  const RealType* mass_flux_x,
						  const RealType* mass_flux_y,
						  RealType* out_variable) {

 // x min and x max
   
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect but boundary
    index_t ix = 0;

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + nx;

    const index_t node_m1o = NodeNodeM1O(node_sym, iy, nx);
    const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
    const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
    const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
    const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
    const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
    const index_t cellm1p1 = NodeCellM1P1(node_sym, iy, nx);
    const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

    const index_t facexm1m1 = NodeFaceXM1M1(node_sym, iy, nx);
    const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
    const index_t facexm1p1 = NodeFaceXM1P1(node_sym, iy, nx);
    const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
    const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
    const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
    const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
    const index_t faceym1oo = NodeFaceYM1OO(node_sym, iy, nx);
    const index_t faceym1p1 = NodeFaceYM1P1(node_sym, iy, nx);
    const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
    const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
    const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "direct_projection_nodal_o2_2d_computation.h"

    out_variable[node_ooo] = out_nodal_variable ;
    out_variable[node_sym] = out_nodal_variable;
  }


  //ymin and y max

  index_t iy = 0;
  index_t iy_sym = ny;
  for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + ny * (nx + 1);
   
    out_variable[node_ooo] = 0.;
    out_variable[node_sym] = 0.;
  }


  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    const index_t node_ooo = node_cmm;

    out_variable[node_cmm] = 0.;
    out_variable[node_cmp] = 0.;
    out_variable[node_cpm] = 0.;
    out_variable[node_cpp] = 0.;
  } // end corners

} // end ProjectNodalIntensiveVariableOrder2BoundaryRtDirect




  //Wall conditions


void ProjectNodalIntensiveVariableUxWallBoundaryDirect(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
					    const RealType* mass_flux_x,
					    const RealType* mass_flux_y,
					    RealType* out_vx) {


  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect but boundary
  // x min 
    for (index_t ix = 0; ix < halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      out_vx[node_ooo] = 0.; 
 
    }
  // x max
    for (index_t ix = nx; ix < nx + 1; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      out_vx[node_ooo] = 0.; 

    }
  } // end X boundary 


  //ymin
  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
	
      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);   
      
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

 
      const index_t facexm1p1 = NodeFaceXM1P1(node_ooo, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);

      const RealType in_vx_ooo = in_vx[node_ooo];
      const RealType in_vx_m1o = in_vx[node_m1o];
      const RealType in_vx_p1o = in_vx[node_p1o];
      const RealType in_vx_om1 = in_vx_ooo;
      const RealType in_vx_op1 = in_vx[node_op1];

      const RealType cell_mass_m1p1 = lag_cell_mass[cellm1p1];
      const RealType cell_mass_p1p1 = lag_cell_mass[cellp1p1];

      const RealType out_mass_m1p1 = out_cell_mass[cellm1p1];
      const RealType out_mass_p1p1 = out_cell_mass[cellp1p1];

      const RealType mass_flux_x_m1p1 = mass_flux_x[facexm1p1];
      const RealType mass_flux_x_oop1 = mass_flux_x[facexoop1];
      const RealType mass_flux_x_p1p1 = mass_flux_x[facexp1p1];
      const RealType mass_flux_y_m1oo = mass_flux_y[faceym1oo];
      const RealType mass_flux_y_m1p1 = mass_flux_y[faceym1p1];
      const RealType mass_flux_y_p1p1 = mass_flux_y[faceyp1p1];
      const RealType mass_flux_y_p1oo = mass_flux_y[faceyp1oo];

      const RealType lag_node_mass_ooo = 0.5 * (cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.5 * (out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux_x = 0.5 * (mass_flux_x_m1p1 + mass_flux_x_oop1);
      
      const RealType next_dual_mass_flux_x = 0.5 * (mass_flux_x_p1p1 + mass_flux_x_oop1);  

      const RealType prev_dual_mass_flux_y = 0.5 * (mass_flux_y_m1oo + mass_flux_y_p1oo);
      
      const RealType next_dual_mass_flux_y = 0.25 * (mass_flux_y_p1oo + mass_flux_y_m1oo + mass_flux_y_p1p1 + mass_flux_y_m1p1);
      
#include "reconstruct_dual_variable_xy_computation.h"      

      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x - moment_x_flux_next_x  + moment_x_flux_prev_y - moment_x_flux_next_y ;


      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 

    }
  }

  //ymax
  for (index_t iy = ny; iy < ny + 1; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
	
      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);   
      
      const index_t cellm1m1 = NodeCellM1M1(node_ooo, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
 
      const index_t facexm1m1 = NodeFaceXM1M1(node_ooo, iy, nx);
      const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
      const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1m1 = NodeFaceYM1M1(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);

      const RealType in_vx_ooo = in_vx[node_ooo];
      const RealType in_vx_m1o = in_vx[node_m1o];
      const RealType in_vx_p1o = in_vx[node_p1o];
      const RealType in_vx_op1 = in_vx_ooo;
      const RealType in_vx_om1 = in_vx[node_om1];

      const RealType cell_mass_m1m1 = lag_cell_mass[cellm1m1];
      const RealType cell_mass_p1m1 = lag_cell_mass[cellp1m1];

      const RealType out_mass_m1m1 = out_cell_mass[cellm1m1];
      const RealType out_mass_p1m1 = out_cell_mass[cellp1m1];

      const RealType mass_flux_x_m1m1 = mass_flux_x[facexm1m1];
      const RealType mass_flux_x_oom1 = mass_flux_x[facexoom1];
      const RealType mass_flux_x_p1m1 = mass_flux_x[facexp1m1];
      const RealType mass_flux_y_m1oo = mass_flux_y[faceym1oo];
      const RealType mass_flux_y_m1m1 = mass_flux_y[faceym1m1];
      const RealType mass_flux_y_p1m1 = mass_flux_y[faceyp1m1];
      const RealType mass_flux_y_p1oo = mass_flux_y[faceyp1oo];

      const RealType lag_node_mass_ooo = 0.5 * (cell_mass_m1m1 +  cell_mass_p1m1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.5 * (out_mass_m1m1 + out_mass_p1m1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux_x = 0.5 * (mass_flux_x_m1m1 + mass_flux_x_oom1);
      
      const RealType next_dual_mass_flux_x = 0.5 * (mass_flux_x_p1m1 + mass_flux_x_oom1);

      const RealType prev_dual_mass_flux_y = 0.25 * (mass_flux_y_p1m1 + mass_flux_y_m1m1 + mass_flux_y_m1oo + mass_flux_y_p1oo);
      
      const RealType next_dual_mass_flux_y = 0.5 * (mass_flux_y_p1oo + mass_flux_y_m1oo);   
      
#include "reconstruct_dual_variable_xy_computation.h"      

      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x - moment_x_flux_next_x + moment_x_flux_prev_y - moment_x_flux_next_y ;

      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 

    } // end Y boundary
  }

  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    out_vx[node_cmm] = 0.;
    out_vx[node_cpm] = 0.;
    out_vx[node_cpp] = 0.;
    out_vx[node_cmp] = 0.;
  }


} // end ProjectNodalIntensiveVariableUxWallBoundaryDirect



void ProjectNodalIntensiveVariableUyWallBoundaryDirect(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
					    const RealType* mass_flux_x,
					    const RealType* mass_flux_y,
					    RealType* out_vx) {

  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect but boundary
  // x min 
    for (index_t ix = 0; ix < halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);     
      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);

      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
      const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

      const RealType in_vx_ooo = in_vx[node_ooo];
      const RealType in_vx_p1o = in_vx[node_p1o];
      const RealType in_vx_m1o = in_vx_ooo;
      const RealType in_vx_om1 = in_vx[node_om1];
      const RealType in_vx_op1 = in_vx[node_op1];

      const RealType cell_mass_p1m1 = lag_cell_mass[cellp1m1];
      const RealType cell_mass_p1p1 = lag_cell_mass[cellp1p1];

      const RealType out_mass_p1m1 = out_cell_mass[cellp1m1];
      const RealType out_mass_p1p1 = out_cell_mass[cellp1p1];

      const RealType mass_flux_x_oom1 = mass_flux_x[facexoom1];
      const RealType mass_flux_x_oop1 = mass_flux_x[facexoop1];
      const RealType mass_flux_x_p1p1 = mass_flux_x[facexp1p1];
      const RealType mass_flux_x_p1m1 = mass_flux_x[facexp1m1];
      const RealType mass_flux_y_p1p1 = mass_flux_y[faceyp1p1];
      const RealType mass_flux_y_p1oo = mass_flux_y[faceyp1oo];
      const RealType mass_flux_y_p1m1 = mass_flux_y[faceyp1m1];

      const RealType lag_node_mass_ooo = 0.5 * (cell_mass_p1m1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.5 * (out_mass_p1m1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux_x = 0.5 * (mass_flux_x_oom1 + mass_flux_x_oop1);
      
      const RealType next_dual_mass_flux_x = 0.25 * (mass_flux_x_p1m1 + mass_flux_x_oom1 + mass_flux_x_p1p1 + mass_flux_x_oop1);  
     
      const RealType prev_dual_mass_flux_y = 0.5 * (mass_flux_y_p1m1 + mass_flux_y_p1oo);
      
      const RealType next_dual_mass_flux_y = 0.5 * (mass_flux_y_p1oo + mass_flux_y_p1p1);  
 
#include "reconstruct_dual_variable_xy_computation.h"      

      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x - moment_x_flux_next_x + moment_x_flux_prev_y - moment_x_flux_next_y;
   
      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 
 
    }
  // x max
    for (index_t ix = nx; ix < nx + 1; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);

      const index_t cellm1m1 = NodeCellM1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);

      const index_t facexm1m1 = NodeFaceXM1M1(node_ooo, iy, nx);
      const index_t facexm1p1 = NodeFaceXM1P1(node_ooo, iy, nx);
      const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1m1 = NodeFaceYM1M1(node_ooo, iy, nx);

      const RealType in_vx_ooo = in_vx[node_ooo];
      const RealType in_vx_m1o = in_vx[node_m1o];
      const RealType in_vx_p1o = in_vx_ooo;
      const RealType in_vx_om1 = in_vx[node_om1];
      const RealType in_vx_op1 = in_vx[node_op1];

      const RealType cell_mass_m1m1 = lag_cell_mass[cellm1m1];
      const RealType cell_mass_m1p1 = lag_cell_mass[cellm1p1];

      const RealType out_mass_m1m1 = out_cell_mass[cellm1m1];
      const RealType out_mass_m1p1 = out_cell_mass[cellm1p1];

      const RealType mass_flux_x_oom1 = mass_flux_x[facexoom1];
      const RealType mass_flux_x_oop1 = mass_flux_x[facexoop1];
      const RealType mass_flux_x_m1p1 = mass_flux_x[facexm1p1];
      const RealType mass_flux_x_m1m1 = mass_flux_x[facexm1m1];
      const RealType mass_flux_y_m1p1 = mass_flux_y[faceym1p1];
      const RealType mass_flux_y_m1oo = mass_flux_y[faceym1oo];
      const RealType mass_flux_y_m1m1 = mass_flux_y[faceym1m1];  

      const RealType lag_node_mass_ooo = 0.5 * (cell_mass_m1m1 +  cell_mass_m1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.5 * (out_mass_m1m1 + out_mass_m1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType next_dual_mass_flux_x = 0.5 * (mass_flux_x_oom1 + mass_flux_x_oop1);
      
      const RealType prev_dual_mass_flux_x = 0.25 * (mass_flux_x_m1m1 + mass_flux_x_oom1 + mass_flux_x_m1p1 + mass_flux_x_oop1);

      const RealType prev_dual_mass_flux_y = 0.5 * (mass_flux_y_m1m1 + mass_flux_y_m1oo);
      
      const RealType next_dual_mass_flux_y = 0.5 * (mass_flux_y_m1oo + mass_flux_y_m1p1);
      
#include "reconstruct_dual_variable_xy_computation.h"      

      const RealType in_moment_x_ooo = lag_node_mass_ooo * in_vx_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x - moment_x_flux_next_x + moment_x_flux_prev_y - moment_x_flux_next_y;
   
      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 

    }
  } // end X boundary 


  //ymin
  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
	
      out_vx[node_ooo] = 0.; 

    }
  }
  //ymax
  for (index_t iy = ny; iy < ny + 1; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
	
      out_vx[node_ooo] = 0.; 

    } // end Y boundary
  }

  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    out_vx[node_cmm] = 0.;
    out_vx[node_cpm] = 0.;
    out_vx[node_cpp] = 0.;
    out_vx[node_cmp] = 0.;
  }

} // end ProjectNodalIntensiveVariableUyWallBoundaryDirect




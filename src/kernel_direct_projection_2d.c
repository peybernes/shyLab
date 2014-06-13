// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"

void ProjectMassDirect(index_t nx, 
			       index_t ny, 
			       const RealType* in_cell_mass,
			       const RealType* mass_flux_x, 
			       const RealType* mass_flux_y, 
		               RealType* out_cell_mass) {
   

  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);

      // dmass
      const RealType mass_flux_prev_x = mass_flux_x[prev_face_x];
      const RealType mass_flux_next_x = mass_flux_x[next_face_x];
      const RealType mass_flux_prev_y = mass_flux_y[prev_face_y];
      const RealType mass_flux_next_y = mass_flux_y[next_face_y];

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      
      const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev_x + mass_flux_prev_y - mass_flux_next_x - mass_flux_next_y;

      out_cell_mass[cell_ooo] = out_cell_mass_ooo; 
      
    }

  }  

}

void MassProjectIntensiveVariableDirect(index_t nx, 
				        index_t ny, 
				        const RealType* RESTRICT in_cell_mass,
				        //const RealType* RESTRICT mass_flux_x,
				        //const RealType* RESTRICT mass_flux_y,
				        const RealType* RESTRICT in_cell_variable,
				        const RealType* RESTRICT in_variable_flux_x,
				        const RealType* RESTRICT in_variable_flux_y,
				        const RealType* RESTRICT out_cell_mass,
				        RealType* RESTRICT out_cell_variable) {
  
  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;

      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);      

      // dmass
      //      const RealType mass_flux_prev_x = mass_flux_x[prev_face_x];
      //      const RealType mass_flux_next_x = mass_flux_x[next_face_x];
      //      const RealType mass_flux_prev_y = mass_flux_y[prev_face_y];
      //      const RealType mass_flux_next_y = mass_flux_y[next_face_y];

      // reconstruct face variable
      const RealType face_variable_prev_x = in_variable_flux_x[prev_face_x];
      const RealType face_variable_next_x = in_variable_flux_x[next_face_x];
      const RealType face_variable_prev_y = in_variable_flux_y[prev_face_y];
      const RealType face_variable_next_y = in_variable_flux_y[next_face_y];

      /// mass and variable
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      const RealType out_cell_mass_ooo = out_cell_mass[cell_ooo];
      
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo];
      
      const RealType out_cell_variable_ooo = 
	in_cell_mass_ooo * in_cell_variable_ooo +
	//	mass_flux_prev_x * face_variable_prev_x + mass_flux_prev_y * face_variable_prev_y -
	//        mass_flux_next_x * face_variable_next_x - mass_flux_next_y * face_variable_next_y;
	face_variable_prev_x + face_variable_prev_y -
        face_variable_next_x - face_variable_next_y;

      assert(0.0 < out_cell_mass_ooo);
      out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo; 
      
    }

  }  

}

void ProjectNodalIntensiveVariableDirect(index_t nx, 
				         index_t ny, 
				         index_t halo_width,
				         const RealType* RESTRICT lag_cell_mass,
				         const RealType* RESTRICT out_cell_mass,
				         const RealType* RESTRICT in_vx,
				         const RealType* RESTRICT in_vy,
				         const RealType* RESTRICT mass_flux_x,
				         const RealType* RESTRICT mass_flux_y,
				         RealType* RESTRICT out_vx,
				         RealType* RESTRICT out_vy) {

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_ooo, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t facexm1m1 = NodeFaceXM1M1(node_ooo, iy, nx);
      const index_t facexoom1 = NodeFaceXOOM1(node_ooo, iy, nx);
      const index_t facexm1p1 = NodeFaceXM1P1(node_ooo, iy, nx);
      const index_t facexoop1 = NodeFaceXOOP1(node_ooo, iy, nx);
      const index_t facexp1p1 = NodeFaceXP1P1(node_ooo, iy, nx);
      const index_t facexp1m1 = NodeFaceXP1M1(node_ooo, iy, nx);
      const index_t faceym1m1 = NodeFaceYM1M1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
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
      const RealType in_moment_y_ooo = lag_node_mass_ooo * in_vy_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;
     
      const RealType out_moment_y_ooo = 
	in_moment_y_ooo + moment_y_flux_prev_x + moment_y_flux_prev_y - moment_y_flux_next_x - moment_y_flux_next_y;


      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 
      out_vy[node_ooo] = out_moment_y_ooo / out_node_mass_ooo;  

     // printf("INFO: nx=%d, ny=%d, iy=%d, prev_node=%d, next_node=%d, vx_prev=%lf, vx_next=%lf, volume_flux=%lf\n", 
      // 	     nx, ny, iy, prev_node, next_node, vx_prev, vx_next, volume_flux);

    }
  }  

}

void ProjectNodalIntensiveVariablePeriodicBoundaryDirect(index_t nx, 
				                 index_t ny, 
				                 index_t halo_width,
				                 const RealType* RESTRICT lag_cell_mass,
				                 const RealType* RESTRICT out_cell_mass,
				                 const RealType* RESTRICT in_vx,
				                 const RealType* RESTRICT in_vy,
				                 const RealType* RESTRICT mass_flux_x,
				                 const RealType* RESTRICT mass_flux_y,
				                 RealType* RESTRICT out_vx,
				                 RealType* RESTRICT out_vy) {

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
      const RealType in_moment_y_ooo = lag_node_mass_ooo * in_vy_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;
     
      const RealType out_moment_y_ooo = 
	in_moment_y_ooo + moment_y_flux_prev_x + moment_y_flux_prev_y - moment_y_flux_next_x - moment_y_flux_next_y;


      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 
      out_vy[node_ooo] = out_moment_y_ooo / out_node_mass_ooo;

      out_vx[node_sym] = out_moment_x_ooo / out_node_mass_ooo; 
      out_vy[node_sym] = out_moment_y_ooo / out_node_mass_ooo;  
 
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
      const RealType in_moment_y_ooo = lag_node_mass_ooo * in_vy_ooo; 

      const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;
     
      const RealType out_moment_y_ooo = 
	in_moment_y_ooo + moment_y_flux_prev_x + moment_y_flux_prev_y - moment_y_flux_next_x - moment_y_flux_next_y;


      out_vx[node_ooo] = out_moment_x_ooo / out_node_mass_ooo; 
      out_vy[node_ooo] = out_moment_y_ooo / out_node_mass_ooo;

      out_vx[node_sym] = out_moment_x_ooo / out_node_mass_ooo; 
      out_vy[node_sym] = out_moment_y_ooo / out_node_mass_ooo;  

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
    const RealType in_moment_y_ooo = lag_node_mass_ooo * in_vy_ooo; 

    const RealType out_moment_x_ooo = 
	in_moment_x_ooo + moment_x_flux_prev_x + moment_x_flux_prev_y - moment_x_flux_next_x - moment_x_flux_next_y;
     
    const RealType out_moment_y_ooo = 
	in_moment_y_ooo + moment_y_flux_prev_x + moment_y_flux_prev_y - moment_y_flux_next_x - moment_y_flux_next_y;



    const RealType out_v_x = out_moment_x_ooo / out_node_mass_ooo; 
    const RealType out_v_y = out_moment_y_ooo / out_node_mass_ooo; 

    out_vx[node_cmm] = out_v_x;
    out_vx[node_cpm] = out_v_x;
    out_vx[node_cmp] = out_v_x;
    out_vx[node_cpp] = out_v_x;

    out_vy[node_cmm] = out_v_y;
    out_vy[node_cpm] = out_v_y;
    out_vy[node_cmp] = out_v_y;
    out_vy[node_cpp] = out_v_y;

  }
}



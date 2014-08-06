// -*- c++ -*- (for emacs users)
//#include <likwid.h> // perf profiling

#include "kernel_ad_projection_2d.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"

void ComputeDirectionalLagrangianQuantitiesX(index_t nx, 
					     index_t ny, 
					     RealType dt,
					     RealType dx,
					     RealType dy,
					     const RealType* RESTRICT in_vx,
					     const RealType* RESTRICT cell_mass,
					     RealType* RESTRICT volume_fluxes,
					     RealType* RESTRICT directional_lagrangian_volume,
					     RealType* RESTRICT directional_lagrangian_density) {
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

	const RealType vx_prev = in_vx[prev_node];                    //1 Load
	const RealType vx_next = in_vx[next_node];                    // 0 load (in cache?)

	const RealType half = 0.5;

	const RealType face_velocity = half * (vx_prev + vx_next);    // 1 FMA
      
	const RealType volume_flux = dt * face_velocity * dy;         // 2 MUL
      
	volume_fluxes[face_ooo] = volume_flux;                        // 1 store                   

      }
      //likwid_markerStopRegion("directionnnal_lagrangianXpart1");
    }  
    
#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianXpart2");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

	const RealType volume_flux_prev = volume_fluxes[prev_face];    //1 Load
	const RealType volume_flux_next = volume_fluxes[next_face];     // 0 load (in cache?)
	const RealType current_mass = cell_mass[cell_ooo];             //1 Load

	const RealType current_directional_lagrangian_volume =         // 1 FMA   1 ADD
	  (dx * dy) - volume_flux_prev + volume_flux_next;

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


void ComputeDirectionalLagrangianQuantitiesY(index_t nx, 
					     index_t ny, 
					     RealType dt,
					     RealType dx,
					     RealType dy,
					     const RealType* RESTRICT in_vy,
					     const RealType* RESTRICT cell_mass,
					     RealType* RESTRICT volume_fluxes,
					     RealType* RESTRICT directional_lagrangian_volume,
					     RealType* RESTRICT directional_lagrangian_density) {
#pragma omp parallel
  {
#pragma omp for
    for (index_t iy = 0; iy < ny + 1; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianYpart1");
       //#pragma unroll (UnrollFactor)
#pragma simd
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = (nx * iy) + ix;
      
	const index_t prev_node = FaceYNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceYNodeP1(face_ooo, iy, nx);

	const RealType vy_prev = in_vy[prev_node];  // 1 Load
	const RealType vy_next = in_vy[next_node];  // 1 Load

	const RealType half = 0.5;
	
	const RealType face_velocity = half * (vy_prev + vy_next); //1 FMA
      
	const RealType volume_flux = dt * face_velocity * dx;      // 2 MUL
      
	volume_fluxes[face_ooo] = volume_flux;     // 1 Store

     }
      //likwid_markerStopRegion("directionnnal_lagrangianYpart1");
    }  

#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianYpart2");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	const RealType volume_flux_prev = volume_fluxes[prev_face]; // 1 load
	const RealType volume_flux_next = volume_fluxes[next_face]; // 1 load
	const RealType current_mass = cell_mass[cell_ooo];         // 1 load

	const RealType current_directional_lagrangian_volume =      // 1 fma  1 add
	  (dx * dy) - volume_flux_prev + volume_flux_next;
      
	const RealType current_directional_lagrangian_density =     // 1 Div
	  current_mass / current_directional_lagrangian_volume;

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume; //1 store
      
	directional_lagrangian_density[cell_ooo] = current_directional_lagrangian_density; // 1 store

	assert(0.0 < current_directional_lagrangian_volume);
	assert(0.0 < current_directional_lagrangian_density);
      
      }
      //likwid_markerStopRegion("directionnnal_lagrangianYpart2");
    }  
  }
}// total function manual optimist count  5 Load + 3 Store    2 MUL 2 FMA 1 Add 1 DIV  


void CheckFluxPeriodicalPropertyX(index_t nx,
				  index_t ny,
				  const RealType* RESTRICT flux) {

  for (index_t iy = 0; iy < ny; ++iy) {
    
    const index_t first_face = (nx + 1) * iy + 0;
    const index_t lasst_face = (nx + 1) * iy + nx;

    const RealType first_flux = flux[first_face];
    const RealType lasst_flux = flux[lasst_face];

    assert(flux[first_face] == flux[lasst_face]);
      
    
  }

}


void CheckFluxPeriodicalPropertyY(index_t nx,
				  index_t ny,
				  const RealType* RESTRICT flux) {

  for (index_t ix = 0; ix < nx; ++ix) {
    
    const index_t first_face = ix;
    const index_t lasst_face = (nx * ny) + ix;

    const RealType first_flux = flux[first_face];
    const RealType lasst_flux = flux[lasst_face];
    
    if (first_flux != lasst_flux) {

      // printf("ERROR: nx=%d, ny=%d, ix=%d, first_flux=%lf, lasst_flux=%lf\n", nx, ny, ix, first_flux, lasst_flux);
      assert(flux[first_face] == flux[lasst_face]);
    
    }
  }

}


void ReconstructMassFluxOrder1X(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType* RESTRICT volume_fluxes,
				const RealType* RESTRICT cell_variable,
				RealType* RESTRICT mass_flux) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      
#include "reconstruct_mass_flux_computation.h"
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }    

}


void ReconstructMassFluxOrder1Y(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType* RESTRICT volume_fluxes,
				const RealType* RESTRICT cell_variable,
				RealType* RESTRICT mass_flux) {

#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      
#include "reconstruct_mass_flux_computation.h"
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }    
}


void ReconstructIntensiveVariableFluxOrder1X(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* RESTRICT mass_flux,
					     const RealType* RESTRICT cell_variable,
					     RealType* RESTRICT variable_flux) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      
#include "reconstruct_variable_computation.h"
      
      variable_flux[face_ooo] = variable_flux_ooo;

    }
  }  
  
}


void ReconstructIntensiveVariableFluxOrder1Y(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* RESTRICT mass_flux,
					     const RealType* RESTRICT cell_variable,
					     RealType* RESTRICT variable_flux) {

#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      
#include "reconstruct_variable_computation.h"
      
      variable_flux[face_ooo] = variable_flux_ooo;

    }
  }    

}


void ProjectMassX(index_t nx, 
		  index_t ny, 
		  const RealType* RESTRICT in_cell_mass,
		  const RealType* RESTRICT mass_flux,
		  RealType* RESTRICT out_cell_mass) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
#pragma simd //force inner loop vectorization  no real perf difference when running more than 2 threads 
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

      // dmass
      const RealType mass_flux_prev = mass_flux[prev_face]; // 1 load
      const RealType mass_flux_next = mass_flux[next_face]; // 0 load (in cache)

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo]; // 1 load
      
      const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next; // 2 add

      out_cell_mass[cell_ooo] = out_cell_mass_ooo;     //1 store
      
    }
  }  

} // total function manual optimist count  2 Load   1 Store    2 Add


void ProjectMassY(index_t nx, 
		  index_t ny, 
		  const RealType* RESTRICT in_cell_mass,
		  const RealType* RESTRICT mass_flux,
		  RealType* RESTRICT out_cell_mass) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {  
    #pragma simd //force inner loop vectorization  no real perf difference when running more than 2 threads  
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

      // dmass
      const RealType mass_flux_prev = mass_flux[prev_face];  // 1 Load
      const RealType mass_flux_next = mass_flux[next_face];  // 1 Load

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo]; // 1 Load
      
      const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next; // 2 Add

      out_cell_mass[cell_ooo] = out_cell_mass_ooo;       // 1 Store
      
    }

  }  

} //total function  manual optimist count  3 load 1 store    2 add


void MassProjectIntensiveVariableX(index_t nx, 
				   index_t ny, 
				   const RealType* RESTRICT in_cell_mass,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT in_cell_variable,
				   const RealType* RESTRICT in_variable_flux,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT out_cell_mass,
				   RealType* RESTRICT out_cell_variable) {
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("massProjectIntensiveX");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

      // dmass
      //const RealType mass_flux_prev = mass_flux[prev_face]; // 1 load
      //const RealType mass_flux_next = mass_flux[next_face]; // 0 load (cache)

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo]; // 1 load
      
      /*const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next; // 2 add

        out_cell_mass[cell_ooo] = out_cell_mass_ooo; // 1 Store
      */
      const RealType out_cell_mass_ooo = out_cell_mass[cell_ooo];

      // reconstruct face variable
      const RealType face_variable_prev = in_variable_flux[prev_face]; // 1 load
      const RealType face_variable_next = in_variable_flux[next_face]; // 0 load (cache)

      /// variable     
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo]; // 1 load
      
      const RealType out_cell_variable_ooo = 
	in_cell_mass_ooo * in_cell_variable_ooo +              //1 FMA   1 Add
	//	mass_flux_prev * face_variable_prev -
	//        mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      assert( (0.0 < out_cell_mass_ooo) || (0.0 == out_cell_mass_ooo) );

      if (out_cell_mass_ooo == 0.0) {					
	out_cell_variable[cell_ooo] = 0.0;
      } else {
	out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo;
      }
      //      assert(0.0 < out_cell_mass_ooo);
      //      out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo;  // 1 Div
                                                                                // 1 Store
    }
    //likwid_markerStopRegion("massProjectIntensiveX");
  }  

}// Total function manual optimist count  4 load  2 store     1 Fma  3 add  1 div  


void MassProjectIntensiveVariableY(index_t nx, 
				   index_t ny, 
				   const RealType* RESTRICT in_cell_mass,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT in_cell_variable,
				   const RealType* RESTRICT in_variable_flux,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT out_cell_mass,
				   RealType* RESTRICT out_cell_variable) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("massProjectIntensiveY");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

      // dmass
      //const RealType mass_flux_prev = mass_flux[prev_face]; // 1 load
      //const RealType mass_flux_next = mass_flux[next_face]; // 1 load

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo]; // 1 load
      /*
      const RealType out_cell_mass_ooo = // 2 add
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next;

      out_cell_mass[cell_ooo] = out_cell_mass_ooo; // 1 store
      */
      const RealType out_cell_mass_ooo = out_cell_mass[cell_ooo];

      // reconstruct face variable
      const RealType face_variable_prev = in_variable_flux[prev_face]; // 1 load
      const RealType face_variable_next = in_variable_flux[next_face]; // 1 load

      /// variable
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo]; // 1 load
      
      const RealType out_cell_variable_ooo = in_cell_mass_ooo * in_cell_variable_ooo +   // 1 FMA   1 ADD
	//mass_flux_prev * face_variable_prev -
        //mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      assert( (0.0 < out_cell_mass_ooo) || (0.0 == out_cell_mass_ooo) );

      if (out_cell_mass_ooo == 0.0) {					
	out_cell_variable[cell_ooo] = 0.0;
      } else {
	out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo;
      }

      //      assert(0.0 < out_cell_mass_ooo);
      //      out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo;    // 1 Div   
                                                                                 //1 Store
      
    }
    //likwid_markerStopRegion("massProjectIntensiveY");
  }  

} //total function  manual optimist count - 6 load  2  store    1 fma  3 add 1 div

void ProjectNodalIntensiveVariableX(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT lag_cell_mass,
				    const RealType* RESTRICT out_cell_mass,
				    const RealType* RESTRICT in_vx,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_vx) {

  //Boucle sur les noeuds
#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    //likwid_markerStartRegion("massProjectNodalIntensiveX");
#pragma unroll (UnrollFactor)
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
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

      const RealType one_quarter = 0.25;

#include "ad_nodal_projection_2d_X_data_load.h"

      const RealType lag_node_mass_ooo = one_quarter * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);

      const RealType out_node_mass_ooo = one_quarter * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = one_quarter * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = one_quarter * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;



      out_vx[node_ooo] = out_moment_ooo / out_node_mass_ooo; 

    }
    //likwid_markerStopRegion("massProjectNodalIntensiveX");
 }  

}


void ProjectNodalIntensiveVariableY(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT lag_cell_mass,
				    const RealType* RESTRICT out_cell_mass, 
				    const RealType* RESTRICT in_vy,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_vy) {

  //Boucle sur les noeuds
#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    //likwid_markerStartRegion("massProjectNodalIntensiveY");
#pragma unroll (UnrollFactor)
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_ooo, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t faceym1m1 = NodeFaceYM1M1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

      const RealType one_quarter = 0.25;

#include "ad_nodal_projection_2d_Y_data_load.h"

      const RealType lag_node_mass_ooo = one_quarter * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);

      const RealType out_node_mass_ooo = one_quarter * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = one_quarter * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = one_quarter * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;


      out_vy[node_ooo] = out_moment_ooo / out_node_mass_ooo; 
    
    }
    //likwid_markerStopRegion("massProjectNodalIntensiveY");
 }  

}


void ReconstructGradientX(index_t nx,
			  index_t ny,
			  RealType dx,
			  RealType dy,
			  const RealType* RESTRICT volume_fluxes,
			  const RealType* RESTRICT lag_variable,
			  RealType* RESTRICT gradient_variable) {

  
 
  
  

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("gradX"); // 0 over Head perf measurment   
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx - 1; ++ix) {
    
     


      const int cell_ooo = iy * nx + ix;
      const int cell_m1o = CellCellM1O(cell_ooo, nx);
      const int cell_p1o = CellCellP1O(cell_ooo, nx);
      
      const int face_m2o = CellFaceM1O(cell_m1o, iy, nx);
      const int face_m1o = CellFaceM1O(cell_ooo, iy, nx);
      const int face_p1o = CellFaceP1O(cell_ooo, iy, nx);
      const int face_p2o = CellFaceP1O(cell_p1o, iy, nx);
      //could be put in --.h file
      //data load
      const RealType variable_m1o = lag_variable[cell_m1o]; // 1 load
      const RealType variable_ooo = lag_variable[cell_ooo]; // 0 load
      const RealType variable_p1o = lag_variable[cell_p1o]; // 0 load
      const RealType d_vol_m2o = volume_fluxes[face_m2o]; // 1 load
      const RealType d_vol_m1o = volume_fluxes[face_m1o]; // 0 load
      const RealType d_vol_p1o = volume_fluxes[face_p1o]; // 0 load
      const RealType d_vol_p2o = volume_fluxes[face_p2o]; // 0 load

      const RealType one_over_dy = 1.0 / dy;    // 1div
      const RealType grad_m1o = (variable_ooo - variable_m1o ) /   //2 add  1 FMA 1 div
	(dx + (d_vol_p1o - d_vol_m2o) * one_over_dy);

      const RealType grad_p1o = (variable_p1o - variable_ooo) /     //2 add 1 FMA  1 div
	(dx + (d_vol_p2o - d_vol_m1o) * one_over_dy);

      const RealType limited_grad_variable = VanAlbadaLimiter(grad_m1o, grad_p1o); // VanAlbada  : 4 FMA  3 Mul  1 logic  1 Div     || ( minmod 2 logic  3 mul  1 min  2 abs)
      gradient_variable[cell_ooo] = limited_grad_variable; // 1 store
      

     
    }//inner loop
    //likwid_markerStopRegion("gradX"); // perf counter
  } // outer loop
 
}// function manual optimist count - 2 load 1 store    6 fma  3 mul  4 add   4 div  (1 logic) 


void ReconstructMassFluxOrder2X(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType dx,
				const RealType dy,
				const RealType* RESTRICT volume_fluxes,
				const RealType* RESTRICT cell_density,
				const RealType* RESTRICT cell_density_gradient,
				RealType* RESTRICT mass_flux) {
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("massFluxO2X");
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx; ++ix) {
      
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell];         // 1 load
      const RealType next_cell_variable = cell_density[next_cell];         // 0 load
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell]; // 1 load
      const RealType next_cell_gradient = cell_density_gradient[next_cell]; // 0 load

      const RealType vol_flux = volume_fluxes[face_ooo];          // 1 load
      const RealType vol_flux_m1o = volume_fluxes[face_m1o];      // 0 load
      const RealType vol_flux_p1o = volume_fluxes[face_p1o];      // 0 load

      const RealType half = 0.5;
      
      const RealType one_over_dy = 1.0 / dy;                // 1 div

      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o *  one_over_dy;   // 1 fma
      const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o * one_over_dy;	 // 1 fma  1 add
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dx_lag_prev_corrected; // 1 mul 1 fma 
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dx_lag_next_corrected; // 1 mul 1 fma
      
      const RealType mass_flux_ooo =
	(half * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) +  // 2 abs 3 fma 1 mul
	(half * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo; // 1 store

    }
    //likwid_markerStopRegion("massFluxO2X");
  }

} // total function manual optimist count -  3 load  1 store       7 fma  3 mul 1add   2 abs  1 div


void ReconstructIntensiveVariableFluxOrder2X(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType dx,
					     const RealType dy,
					     const RealType* RESTRICT volume_fluxes,
					     const RealType* RESTRICT mass_flux,
					     const RealType* RESTRICT cell_variable,
					     const RealType* RESTRICT cell_variable_gradient,
					     RealType* RESTRICT variable_flux) {
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    //likwid_markerStartRegion("fluxO2X");
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx; ++ix) {
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;
      
      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];        // 1 load
      const RealType next_cell_variable = cell_variable[next_cell];        // 0 load
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell]; // 1 load
      const RealType next_cell_gradient = cell_variable_gradient[next_cell]; // 0 load

      const RealType mass_flux_face = mass_flux[face_ooo];              // 1 load
      
      const RealType vol_flux_m1o = volume_fluxes[face_m1o];            // 1 load
      const RealType vol_flux_p1o = volume_fluxes[face_p1o];            // 0 load
      
      const RealType half = 0.5;

      const RealType one_over_dy = 1.0 / dy;                  // 1 div
 
      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o * one_over_dy ;  // 1 fma
      const RealType dx_lag_next_corrected = - dx - vol_flux_p1o * one_over_dy ; // 1 fma  1 add	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dx_lag_prev_corrected; // 1 fma 1 mul
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dx_lag_next_corrected; // 1 fma 1 mul
      
      const RealType variable_flux_ooo =
	(half * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) +  // 2 abs  3 fma 1 mul
	(half * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
     
      variable_flux[face_ooo] = variable_flux_ooo; // 1 store

    }
    //likwid_markerStopRegion("fluxO2X");
 }

} //total function manual optimist count  4 load  1 store   7 fma  3 mul  1 add  2abs 1 div 



void ReconstructGradientNodalX(index_t nx,
			       index_t ny,
			       RealType dx,
			       RealType dt,
			       const RealType* RESTRICT predicted_velocity,
			       const RealType* RESTRICT lag_nodal_variable,
			       RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    //likwid_markerStartRegion("gradientNodalX");
#pragma unroll (UnrollFactor)
    for (index_t ix = 1; ix < nx; ++ix) {
    
      const int node_ooo = iy * (nx + 1) + ix;
      const int node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const int node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
      //data load
      const RealType variable_m1o = lag_nodal_variable[node_m1o];   // 1 load
      const RealType variable_ooo = lag_nodal_variable[node_ooo];   // 0 load
      const RealType variable_p1o = lag_nodal_variable[node_p1o];   // 0 load

      const RealType velocity_m1o = predicted_velocity[node_m1o];   // 1 load
      const RealType velocity_ooo = predicted_velocity[node_ooo];   // 0 load
      const RealType velocity_p1o = predicted_velocity[node_p1o];   // 0 load
     

      const RealType grad_m1o = (variable_ooo - variable_m1o ) /    //2 add 1 fma  1 div 
       	(dx + dt * (velocity_ooo - velocity_m1o));

      const RealType grad_p1o = (variable_p1o - variable_ooo ) / //2 add 1 fma  1 div 
       	(dx + dt * (velocity_p1o - velocity_ooo));

      const RealType limited_grad_variable = VanAlbadaLimiter(grad_m1o,grad_p1o) ; // VanAlbada  : 4 FMA  3 Mul  1 logic  1 Div   
      gradient_variable[node_ooo] = limited_grad_variable; // 1 store
      
    }
    //likwid_markerStopRegion("gradientNodalX");
 } 

} //total function manual optimist count  --   2 load  1 store   6 fma  3 mul 4 add  3 div  (1 logic)


void ProjectNodalIntensiveVariableOrder2X(index_t nx, 
					  index_t ny, 
					  index_t halo_width,
					  const RealType dx,
					  const RealType dt,
					  const RealType* RESTRICT lag_cell_mass,
					  const RealType* RESTRICT out_cell_mass,
					  const RealType* RESTRICT u_velocity_pred,
					  const RealType* RESTRICT in_variable,
					  const RealType* RESTRICT gradient_variable,
					  const RealType* RESTRICT mass_flux,
					  RealType* RESTRICT out_variable) {

  //Boucle sur les noeuds
#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    //likwid_markerStartRegion("projectVariableNodalX");
#pragma unroll (UnrollFactor)
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
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

#include "ad_projection_nodal_o2_2d_computation_X.h"
      // 9 load   8 fma  8 add 2 mul   
      //   if statement but similare branch
      //   4 fma ( or 4 fma  and 1-2 add) // end if
      //3 mul  2 add 1 div
      
      out_variable[node_ooo] = out_nodal_variable ; // 1 store
    }
    //likwid_markerStopRegion("projectVariableNodalX");
  }

} //total function manual optimist count -- 9 load 1 strore   12 fma  10 add  5 mul 1 div

void ReconstructGradientY(index_t nx,
			  index_t ny,
			  RealType dx,
			  RealType dy,
			  const RealType* RESTRICT volume_fluxes,
			  const RealType* RESTRICT lag_variable,
			  RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    //likwid_markerStartRegion("gradY");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {

      const int cell_ooo = iy * nx + ix;
      const int cell_om1 = CellCellOM1(cell_ooo, nx);
      const int cell_op1 = CellCellOP1(cell_ooo, nx);
      
      const int face_om2 = CellFaceOM1(cell_om1, iy - 1, nx);
      const int face_om1 = CellFaceOM1(cell_ooo, iy, nx);
      const int face_op1 = CellFaceOP1(cell_ooo, iy, nx);
      const int face_op2 = CellFaceOP1(cell_op1, iy + 1, nx);

      const RealType variable_om1 = lag_variable[cell_om1]; // 1 load
      const RealType variable_ooo = lag_variable[cell_ooo]; // 1 load
      const RealType variable_op1 = lag_variable[cell_op1]; // 1 load
      const RealType d_vol_om2 = volume_fluxes[face_om2];   // 1 load
      const RealType d_vol_om1 = volume_fluxes[face_om1];   // 1 load
      const RealType d_vol_op1 = volume_fluxes[face_op1];   // 1 load
      const RealType d_vol_op2 = volume_fluxes[face_op2];   // 1 load

      const RealType one_over_dx = 1.0 / dx;    // 1 div
	
      const RealType grad_om1 = (variable_ooo - variable_om1 ) / //2 add 1 fma 1 div
	(dy + (d_vol_op1 - d_vol_om2) * one_over_dx);

      const RealType grad_op1 = (variable_op1 - variable_ooo )  /  //2 add 1 fma 1 div
	(dy + (d_vol_op2 - d_vol_om1) *one_over_dx);

      const RealType limited_grad_variable = VanAlbadaLimiter(grad_om1,grad_op1) ;// VanAlbada  : 4 FMA  3 Mul  1 logic  1 Div 
      gradient_variable[cell_ooo] = limited_grad_variable; // 1 store
    }
    //likwid_markerStopRegion("gradY");
  }

} // end ReconstructGradientY
//total function optimist manual count - 7 load 1 store  6 fma  4add 3mul  4div  (logic 1)





void ReconstructMassFluxOrder2Y(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType dx,
				const RealType dy,
				const RealType* RESTRICT volume_fluxes,
				const RealType* RESTRICT cell_density,
				const RealType* RESTRICT cell_density_gradient,
				RealType* RESTRICT mass_flux) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    //likwid_markerStartRegion("massFluxO2Y");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {
      
      
      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell]; // 1 load
      const RealType next_cell_variable = cell_density[next_cell];// 1 load
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell];// 1 load
      const RealType next_cell_gradient = cell_density_gradient[next_cell];// 1 load

      const RealType vol_flux = volume_fluxes[face_ooo];// 1 load
      const RealType vol_flux_om1 = volume_fluxes[face_om1];// 1 load
      const RealType vol_flux_op1 = volume_fluxes[face_op1];// 1 load
      
      const RealType half = 0.5;

      const RealType one_over_dx = 1.0 / dx;                          // 1 div
      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx; // 1fma
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx; // 1 add 1 fma	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected; // 1 mul  1 fma
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected; // 1 mul 1 fma
      
      const RealType mass_flux_ooo =
	(half * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + // 2 abs 3 fma 1 mul 
	(half * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo; // 1 store

    }
    //likwid_markerStopRegion("massFluxO2Y");
  }

} //end ReconstructMassFluxOrder2Y
//total function manual optimist count -- 7 load 1 store -- 7fma   1 add  3 mul  1 div  2 abs



void ReconstructIntensiveVariableFluxOrder2Y(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType dx,
					     const RealType dy,
					     const RealType* RESTRICT volume_fluxes,
					     const RealType* RESTRICT mass_flux,
					     const RealType* RESTRICT cell_variable,
					     const RealType* RESTRICT cell_variable_gradient,
					     RealType* RESTRICT variable_flux) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    //likwid_markerStartRegion("FluxO2Y");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx; ++ix) {
      
      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_variable[prev_cell];  // 1 load
      const RealType next_cell_variable = cell_variable[next_cell]; // 1 load
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];// 1 load
      const RealType next_cell_gradient = cell_variable_gradient[next_cell];// 1 load
 
      const RealType mass_flux_face = mass_flux[face_ooo];// 1 load
      
      const RealType vol_flux_om1 = volume_fluxes[face_om1];// 1 load
      const RealType vol_flux_op1 = volume_fluxes[face_op1];// 1 load
      
      const RealType half = 0.5;
      const RealType one_over_dx = 1.0 / dx;                    //1 div
      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 * one_over_dx; // 1 fma
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 * one_over_dx; // 1add 1 fma	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + half * prev_cell_gradient * dy_lag_prev_corrected; // 1 mul  1 fma
      const RealType next_cell_variable_o2 =  next_cell_variable + half * next_cell_gradient * dy_lag_next_corrected; // 1 mul 1 fma
      
      const RealType variable_flux_ooo =
	(half * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) +  // 2 abs 3 fma 1 mul
	(half * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
      variable_flux[face_ooo] = variable_flux_ooo; // 1 store
  

    }
    //likwid_markerStopRegion("FluxO2Y");
  }
} // end ReconstructIntensiveVariableFluxOrder2Y
// total function manual optimist count - 7 load 1 store-- 7 fma  1add  3mul   1 div  2 abs


void ReconstructGradientNodalY(index_t nx,
			       index_t ny,
			       RealType dy,
			       RealType dt,
			       const RealType* RESTRICT predicted_velocity,
			       const RealType* RESTRICT lag_nodal_variable,
			       RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    //likwid_markerStartRegion("gradientNodalY");
#pragma unroll (UnrollFactor)
    for (index_t ix = 0; ix < nx + 1; ++ix) {
    
      const int node_ooo = iy * (nx + 1) + ix;
      const int node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const int node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      //data load
      const RealType variable_om1 = lag_nodal_variable[node_om1]; // 1 load
      const RealType variable_ooo = lag_nodal_variable[node_ooo]; // 1 load
      const RealType variable_op1 = lag_nodal_variable[node_op1]; // 1 load

      const RealType velocity_om1 = predicted_velocity[node_om1]; // 1 load
      const RealType velocity_ooo = predicted_velocity[node_ooo]; // 1 load
      const RealType velocity_op1 = predicted_velocity[node_op1]; // 1 load
     

      const RealType grad_om1 = (variable_ooo - variable_om1 ) / //1 fma 2 add 1 div
       	(dy + dt * (velocity_ooo - velocity_om1));

      const RealType grad_op1 = (variable_op1 - variable_ooo ) / //1 fma 2 add 1 div
       	(dy + dt * (velocity_op1 - velocity_ooo));

      const RealType limited_grad_variable = VanAlbadaLimiter(grad_om1,grad_op1) ;// VanAlbada  : 4 FMA  3 Mul  (1 logic)  1 Div 
      gradient_variable[node_ooo] = limited_grad_variable; // 1 store
      
    }
    //likwid_markerStopRegion("gradientNodalY");
 } 

} // end ReconstructGradientNodalY
// total function optimist manual count -- 6 load 1 store -- 6 fma  4 add 3 mul 1 div (1logic)


void ProjectNodalIntensiveVariableOrder2Y(index_t nx, 
					  index_t ny, 
					  index_t halo_width,
					  const RealType dy,
					  const RealType dt,
					  const RealType* RESTRICT lag_cell_mass,
					  const RealType* RESTRICT out_cell_mass,
					  const RealType* RESTRICT v_velocity_pred,
					  const RealType* RESTRICT in_variable,
					  const RealType* RESTRICT gradient_variable,
					  const RealType* RESTRICT mass_flux,
					  RealType* RESTRICT out_variable) {

  //Boucle sur les noeuds
#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    //likwid_markerStartRegion("projectVariableNodalY");
#pragma unroll (UnrollFactor)
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;

      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_ooo, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t faceym1m1 = NodeFaceYM1M1(node_ooo, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "ad_projection_nodal_o2_2d_computation_Y.h"
      //16 load  8 fma  8 add   2mul 
      // if - 4 fma (or fma 1-2 add)
      //3mul  2 add 1 div

      out_variable[node_ooo] = out_nodal_variable ; // 1 store
    }
    //likwid_markerStopRegion("projectVariableNodalY");
  }

} //end ProjectNodalIntensiveVariableOrder2Y
//total function manual optimist count -- 16 load 1 store --12 fma  10 add 5mul 1 div



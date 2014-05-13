// -*- c++ -*- (for emacs users)
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
      for (index_t ix = 0; ix < nx + 1; ++ix) {
      
	const index_t face_ooo = ((nx + 1) * iy) + ix;
      
	const index_t prev_node = FaceXNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceXNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];
	const RealType vx_next = in_vx[next_node];

	const RealType face_velocity = 0.5 * (vx_prev + vx_next);
      
	const RealType volume_flux = dt * face_velocity * dy;
      
	volume_fluxes[face_ooo] = volume_flux;

      }
    }  
    
#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

	const RealType volume_flux_prev = volume_fluxes[prev_face];
	const RealType volume_flux_next = volume_fluxes[next_face];
	const RealType current_mass = cell_mass[cell_ooo];

	const RealType current_directional_lagrangian_volume =
	  (dx * dy) - volume_flux_prev + volume_flux_next;

	const RealType current_directional_lagrangian_density =
	  current_mass / current_directional_lagrangian_volume;

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume;
      
	directional_lagrangian_density[cell_ooo] = current_directional_lagrangian_density;


	assert(0.0 < current_directional_lagrangian_volume);
	assert(0.0 < current_directional_lagrangian_density);
      
      }
    }  
  }
}


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
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = (nx * iy) + ix;
      
	const index_t prev_node = FaceYNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceYNodeP1(face_ooo, iy, nx);

	const RealType vy_prev = in_vy[prev_node];
	const RealType vy_next = in_vy[next_node];

	const RealType face_velocity = 0.5 * (vy_prev + vy_next);
      
	const RealType volume_flux = dt * face_velocity * dx;
      
	volume_fluxes[face_ooo] = volume_flux;

      }
    }  

#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	const RealType volume_flux_prev = volume_fluxes[prev_face];
	const RealType volume_flux_next = volume_fluxes[next_face];
	const RealType current_mass = cell_mass[cell_ooo];

	const RealType current_directional_lagrangian_volume =
	  (dx * dy) - volume_flux_prev + volume_flux_next;
      
	const RealType current_directional_lagrangian_density =
	  current_mass / current_directional_lagrangian_volume;

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume;
      
	directional_lagrangian_density[cell_ooo] = current_directional_lagrangian_density;

	assert(0.0 < current_directional_lagrangian_volume);
	assert(0.0 < current_directional_lagrangian_density);
      
      }
    }  
  }
}


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


void ReconstructMassFluxOrder1BoundaryX(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_variable,
					RealType* RESTRICT mass_flux) {

  // xmin
  for (index_t iy = 0; iy < ny; ++iy) {

    for (index_t ix = 0; ix < halo_width; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t prev_cell = nx * (iy + 1) - 1;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_mass_flux_computation.h"
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }

  // xmax
  for (index_t iy = 0; iy < ny; ++iy) {

    for (index_t ix = nx + 1 - halo_width; ix < nx + 1; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);    
      const index_t next_cell = nx * iy;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_mass_flux_computation.h"
      
      mass_flux[face_ooo] = mass_flux_ooo;
      
    }
  }

}


void ReconstructMassFluxOrder1BoundaryY(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_variable,
					RealType* RESTRICT mass_flux) {

  // ymin
  for (index_t iy = 0; iy < halo_width; ++iy) {

    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t prev_cell = (nx * (ny - 1)) + ix;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_mass_flux_computation.h"
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }

  // ymax
  for (index_t iy = ny + 1 - halo_width; iy < ny + 1; ++iy) {

    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);    
      const index_t next_cell = ix;

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


void ReconstructIntensiveVariableFluxOrder1BoundaryX(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     RealType* RESTRICT variable_flux) {

  // xmin
  for (index_t iy = 0; iy < ny; ++iy) {

    for (index_t ix = 0; ix < halo_width; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t prev_cell = nx * (iy + 1) - 1;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_variable_computation.h"
      
      variable_flux[face_ooo] = variable_flux_ooo;

    }
  }

  // xmax
  for (index_t iy = 0; iy < ny; ++iy) {

    for (index_t ix = nx + 1 - halo_width; ix < nx + 1; ++ix) {

      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);    
      const index_t next_cell = nx * iy;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_variable_computation.h"
      
      variable_flux[face_ooo] = variable_flux_ooo;
      
    }
  }
}


void ReconstructIntensiveVariableFluxOrder1BoundaryY(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     RealType* RESTRICT variable_flux) {

  // ymin
  for (index_t iy = 0; iy < halo_width; ++iy) {

    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t prev_cell = (nx * (ny - 1)) + ix;

      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];

#include "reconstruct_variable_computation.h"
      
      variable_flux[face_ooo] = variable_flux_ooo;

    }
  }

  // ymax
  for (index_t iy = ny + 1 - halo_width; iy < ny + 1; ++iy) {

    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);    
      const index_t next_cell = ix;

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
      const RealType mass_flux_prev = mass_flux[prev_face];
      const RealType mass_flux_next = mass_flux[next_face];

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      
      const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next;

      out_cell_mass[cell_ooo] = out_cell_mass_ooo; 
      
    }
  }  

}


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
      const RealType mass_flux_prev = mass_flux[prev_face];
      const RealType mass_flux_next = mass_flux[next_face];

      /// mass
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      
      const RealType out_cell_mass_ooo = 
	in_cell_mass_ooo +  mass_flux_prev - mass_flux_next;

      out_cell_mass[cell_ooo] = out_cell_mass_ooo; 
      
    }

  }  

}


void MassProjectIntensiveVariableX(index_t nx, 
				   index_t ny, 
				   const RealType* RESTRICT in_cell_mass,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT in_cell_variable,
				   const RealType* RESTRICT in_variable_flux,
				   const RealType* RESTRICT out_cell_mass,
				   RealType* RESTRICT out_cell_variable) {
  
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

      // dmass
      //      const RealType mass_flux_prev = mass_flux[prev_face];
      //      const RealType mass_flux_next = mass_flux[next_face];

      // reconstruct face variable
      const RealType face_variable_prev = in_variable_flux[prev_face];
      const RealType face_variable_next = in_variable_flux[next_face];

      /// mass and variable
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      const RealType out_cell_mass_ooo = out_cell_mass[cell_ooo];
      
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo];
      
      const RealType out_cell_variable_ooo = 
	in_cell_mass_ooo * in_cell_variable_ooo +
	//	mass_flux_prev * face_variable_prev -
	//        mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      assert(0.0 < out_cell_mass_ooo);
      out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo; 
      
    }

  }  

}


void MassProjectIntensiveVariableY(index_t nx, 
				   index_t ny, 
				   const RealType* RESTRICT in_cell_mass,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT in_cell_variable,
				   const RealType* RESTRICT in_variable_flux,
				   const RealType* RESTRICT out_cell_mass,
				   RealType* RESTRICT out_cell_variable) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;
      
      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

      // dmass
      // const RealType mass_flux_prev = mass_flux[prev_face];
      // const RealType mass_flux_next = mass_flux[next_face];

      // reconstruct face variable
      const RealType face_variable_prev = in_variable_flux[prev_face];
      const RealType face_variable_next = in_variable_flux[next_face];

      /// mass and variable
      const RealType in_cell_mass_ooo = in_cell_mass[cell_ooo];
      const RealType out_cell_mass_ooo = out_cell_mass[cell_ooo];
      
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo];
      
      const RealType out_cell_variable_ooo = in_cell_mass_ooo * in_cell_variable_ooo +
	//mass_flux_prev * face_variable_prev -
        //mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      assert(0.0 < out_cell_mass_ooo);
      out_cell_variable[cell_ooo] = out_cell_variable_ooo / out_cell_mass_ooo; 
      
    }

  }  

}


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

#include "ad_nodal_projection_2d_X_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);

      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;



      out_vx[node_ooo] = out_moment_ooo / out_node_mass_ooo; 

    }
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

#include "ad_nodal_projection_2d_Y_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);

      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;


      out_vy[node_ooo] = out_moment_ooo / out_node_mass_ooo; 
    
    }
  }  

}


void ProjectNodalIntensiveVariableBoundaryX(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* RESTRICT lag_cell_mass,
					    const RealType* RESTRICT out_cell_mass,  
					    const RealType* RESTRICT in_vx,
					    const RealType* RESTRICT mass_flux,
					    RealType* RESTRICT out_vx) {

  // xmin (xmax is deduced in periodic conditions)

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t ix = 0; ix < halo_width; ++ix) {
    for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) { //problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t node_m1o = NodeNodeM1O(node_sym, iy, nx);
      const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
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

#include "ad_nodal_projection_2d_X_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;
   
      out_vx[node_ooo] = out_moment_ooo / out_node_mass_ooo; 

      out_vx[node_sym] = out_moment_ooo / out_node_mass_ooo;
 
    }
  }  

  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary 

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + (nx + 1) * ny;
      const int iy_sym = ny; 
	
      index_t node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
    
      
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

#include "ad_nodal_projection_2d_X_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);

      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo =lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;


      out_vx[node_ooo] = out_moment_ooo / out_node_mass_ooo; 

      out_vx[node_sym] = out_moment_ooo / out_node_mass_ooo;
 
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


#include "ad_nodal_projection_2d_X_data_load.h"

    const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
    assert(0.0 < lag_node_mass_ooo);
    const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
    assert(0.0 < out_node_mass_ooo);

    const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
    const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

    const RealType in_moment_ooo =lag_node_mass_ooo * in_vx_oo; 

    const RealType out_moment_ooo = 
      in_moment_ooo + moment_flux_prev - moment_flux_next;


    const RealType out_v = out_moment_ooo / out_node_mass_ooo; 


    out_vx[node_cmm] = out_v;
    out_vx[node_cpm] = out_v;
    out_vx[node_cmp] = out_v;
    out_vx[node_cpp] = out_v;
  }
}


void ProjectNodalIntensiveVariableBoundaryY(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* RESTRICT lag_cell_mass,
					    const RealType* RESTRICT out_cell_mass,
					    const RealType* RESTRICT in_vy,
					    const RealType* RESTRICT mass_flux,
					    RealType* RESTRICT out_vy) {

  // ymin (ymax is deduced in periodic conditions)

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect but boundary

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + (nx + 1) * ny;
      const index_t iy_sym = ny;

      const index_t node_om1 = NodeNodeOM1(node_sym, iy_sym, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_sym, iy_sym, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_sym, iy_sym, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy_sym, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_sym, iy_sym, nx);

#include "ad_nodal_projection_2d_Y_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);


      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;

      assert(0.0 < out_node_mass_ooo);

      out_vy[node_ooo] = out_moment_ooo / out_node_mass_ooo; 
 
      out_vy[node_sym] = out_moment_ooo / out_node_mass_ooo;  
 
    }
  }  

  for (index_t ix = 0; ix < halo_width; ++ix) {
    for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect
      
      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + nx;

      const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_sym, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_sym, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_sym, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "ad_nodal_projection_2d_Y_data_load.h"

      const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      assert(0.0 < lag_node_mass_ooo);
      const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
      assert(0.0 < out_node_mass_ooo);
     
      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo =lag_node_mass_ooo * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;


      out_vy[node_ooo] = out_moment_ooo / out_node_mass_ooo; 
 
      out_vy[node_sym] = out_moment_ooo / out_node_mass_ooo;
    }
  }  

  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    const index_t node_ooo = node_cmm;
	
    const index_t node_om1 = (ny - 1) * (nx + 1);
    const index_t node_op1 = nx + 1;
    
      
    const index_t cellm1m1 = nx * ny - 1;
    const index_t cellp1m1 = (ny - 1) * nx;
    const index_t cellm1p1 = nx - 1;
    const index_t cellp1p1 = 0 ;

    const index_t faceym1m1 = ny * nx - 1;
    const index_t faceym1oo = nx - 1;
    const index_t faceym1p1 = nx + nx - 1;
    const index_t faceyp1p1 = nx;
    const index_t faceyp1oo = 0;
    const index_t faceyp1m1 = (ny - 1) * nx;


#include "ad_nodal_projection_2d_Y_data_load.h"

    const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
    assert(0.0 < lag_node_mass_ooo);
    const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
    assert(0.0 < out_node_mass_ooo);
     
    const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
    const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

    const RealType in_moment_ooo = lag_node_mass_ooo * in_vx_oo; 

    const RealType out_moment_ooo = 
      in_moment_ooo + moment_flux_prev - moment_flux_next;

    const RealType out_v = out_moment_ooo / out_node_mass_ooo; 
    
    
    out_vy[node_cmm] = out_v;
    out_vy[node_cpm] = out_v;
    out_vy[node_cmp] = out_v;
    out_vy[node_cpp] = out_v;
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
      const RealType variable_m1o = lag_variable[cell_m1o];
      const RealType variable_ooo = lag_variable[cell_ooo];
      const RealType variable_p1o = lag_variable[cell_p1o];
      const RealType d_vol_m2o = volume_fluxes[face_m2o];
      const RealType d_vol_m1o = volume_fluxes[face_m1o];
      const RealType d_vol_p1o = volume_fluxes[face_p1o];
      const RealType d_vol_p2o = volume_fluxes[face_p2o];


      const RealType grad_m1o = (variable_ooo - variable_m1o ) /
	(dx + (d_vol_p1o - d_vol_m2o) / dy);

      const RealType grad_p1o = (variable_p1o - variable_ooo ) /
	(dx + (d_vol_p2o - d_vol_m1o) / dy);

      const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
      gradient_variable[cell_ooo] = limited_grad_variable;
      
    }
  } 

}


void ReconstructGradientXBoundary(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT volume_fluxes,
				  const RealType* RESTRICT lag_variable,
				  RealType* RESTRICT gradient_variable) {
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
 
    //could be put in --.h file
    //data load
    const RealType variable_m1o = lag_variable[cell_m1o];
    const RealType variable_ooo = lag_variable[cell_ooo];
    const RealType variable_p1o = lag_variable[cell_p1o];
    const RealType d_vol_m2o = volume_fluxes[face_m2o];
    const RealType d_vol_m1o = volume_fluxes[face_m1o];
    const RealType d_vol_p1o = volume_fluxes[face_p1o];
    const RealType d_vol_p2o = volume_fluxes[face_p2o];


    const RealType grad_m1o = (variable_ooo - variable_m1o ) /
      (dx + (d_vol_p1o - d_vol_m2o) / dy);

    const RealType grad_p1o = (variable_p1o - variable_ooo ) /
      (dx + (d_vol_p2o - d_vol_m1o) / dy);

    const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
    gradient_variable[cell_ooo] = limited_grad_variable;   
    
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
 
    //could be put in --.h file
    //data load
    const RealType variable_m1o = lag_variable[cell_m1o];
    const RealType variable_ooo = lag_variable[cell_ooo];
    const RealType variable_p1o = lag_variable[cell_p1o];
    const RealType d_vol_m2o = volume_fluxes[face_m2o];
    const RealType d_vol_m1o = volume_fluxes[face_m1o];
    const RealType d_vol_p1o = volume_fluxes[face_p1o];
    const RealType d_vol_p2o = volume_fluxes[face_p2o];


    const RealType grad_m1o = (variable_ooo - variable_m1o ) /
      (dx + (d_vol_p1o - d_vol_m2o) / dy);

    const RealType grad_p1o = (variable_p1o - variable_ooo ) /
      (dx + (d_vol_p2o - d_vol_m1o) / dy);

    const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
    gradient_variable[cell_ooo] = limited_grad_variable; 
  }

}


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
    for (index_t ix = 1; ix < nx; ++ix) {
      
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;

      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell];
      const RealType next_cell_variable = cell_density[next_cell];
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell];
      const RealType next_cell_gradient = cell_density_gradient[next_cell];

      const RealType vol_flux = volume_fluxes[face_ooo];
      const RealType vol_flux_m1o = volume_fluxes[face_m1o];
      const RealType vol_flux_p1o = volume_fluxes[face_p1o];
      
      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o / dy;
      const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o / dy;	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
      const RealType mass_flux_ooo =
	(0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
	(0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }

}


void ReconstructMassFluxOrder2XBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_density,
					const RealType* RESTRICT cell_density_gradient,
					RealType* RESTRICT mass_flux) {
  
  //x min only //
  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
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
    const RealType vol_flux_m1o = volume_fluxes[face_m1o];
    const RealType vol_flux_p1o = volume_fluxes[face_p1o];
      
    const RealType dx_lag_prev_corrected = dx  - vol_flux_m1o / dy;
    const RealType dx_lag_next_corrected = - dx  - vol_flux_p1o / dy;	
    const RealType prev_cell_variable_o2 = prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
    const RealType mass_flux_ooo =
      (0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
      (0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
    mass_flux[face_ooo] = mass_flux_ooo;
    mass_flux[face_sym] = mass_flux_ooo;
    
  }

}


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
    for (index_t ix = 1; ix < nx; ++ix) {
      
      const index_t face_ooo = ((nx + 1) * iy) + ix;
      
      const index_t prev_cell = FaceXCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceXCellP1(face_ooo, iy, nx);
      const index_t face_m1o = CellFaceM1O( prev_cell, iy, nx);
      const index_t face_p1o = CellFaceP1O( next_cell, iy, nx);
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];
      const RealType next_cell_gradient = cell_variable_gradient[next_cell];

      const RealType mass_flux_face = mass_flux[face_ooo];
      
      const RealType vol_flux_m1o = volume_fluxes[face_m1o];
      const RealType vol_flux_p1o = volume_fluxes[face_p1o];
      
      const RealType dx_lag_prev_corrected = dx - vol_flux_m1o / dy;
      const RealType dx_lag_next_corrected = - dx - vol_flux_p1o / dy;	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
      const RealType variable_flux_ooo =
	(0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
	(0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
      variable_flux[face_ooo] = variable_flux_ooo;

    }
  }

}



void ReconstructIntensiveVariableFluxOrder2XBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType* RESTRICT volume_fluxes,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     const RealType* RESTRICT cell_variable_gradient,
						     RealType* RESTRICT variable_flux) {
  

  for (index_t iy = 0; iy < ny; ++iy) {//problem not vect but boundary
    index_t ix = 0;
      
    const index_t face_ooo = ((nx + 1) * iy) + ix;
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
      
    const RealType vol_flux_m1o = volume_fluxes[face_m1o];
    const RealType vol_flux_p1o = volume_fluxes[face_p1o];
      
    const RealType dx_lag_prev_corrected = dx - vol_flux_m1o / dy;
    const RealType dx_lag_next_corrected = - dx - vol_flux_p1o / dy;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dx_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dx_lag_next_corrected;
      
    const RealType variable_flux_ooo =
      (0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
      (0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
    variable_flux[face_ooo] = variable_flux_ooo;
    variable_flux[face_sym] = variable_flux_ooo;
    
  }

}



void ReconstructGradientNodalX(index_t nx,
			       index_t ny,
			       RealType dx,
			       RealType dt,
			       const RealType* RESTRICT predicted_velocity,
			       const RealType* RESTRICT lag_nodal_variable,
			       RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    for (index_t ix = 1; ix < nx; ++ix) {
    
      const int node_ooo = iy * (nx + 1) + ix;
      const int node_m1o = NodeNodeM1O(node_ooo, iy, nx);
      const int node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
      //data load
      const RealType variable_m1o = lag_nodal_variable[node_m1o];
      const RealType variable_ooo = lag_nodal_variable[node_ooo];
      const RealType variable_p1o = lag_nodal_variable[node_p1o];

      const RealType velocity_m1o = predicted_velocity[node_m1o];
      const RealType velocity_ooo = predicted_velocity[node_ooo];
      const RealType velocity_p1o = predicted_velocity[node_p1o];
     

      const RealType grad_m1o = (variable_ooo - variable_m1o ) /
       	(dx + dt * (velocity_ooo - velocity_m1o));

      const RealType grad_p1o = (variable_p1o - variable_ooo ) /
       	(dx + dt * (velocity_p1o - velocity_ooo));

      const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
      gradient_variable[node_ooo] = limited_grad_variable;
      
    }
  } 

}



void ReconstructGradientNodalXBoundary(index_t nx,
				       index_t ny,
				       RealType dx,
				       RealType dt,
				       const RealType* RESTRICT predicted_velocity,
				       const RealType* RESTRICT lag_nodal_variable,
				       RealType* RESTRICT gradient_variable) {
  //x min
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    index_t ix = 0;
    
    const int node_ooo = iy * (nx + 1) + ix;
    const int node_m1o = node_ooo + nx - 1;
    const int node_p1o = NodeNodeP1O(node_ooo, iy, nx);
    //data load
    const RealType variable_m1o = lag_nodal_variable[node_m1o];
    const RealType variable_ooo = lag_nodal_variable[node_ooo];
    const RealType variable_p1o = lag_nodal_variable[node_p1o];

    const RealType velocity_m1o = predicted_velocity[node_m1o];
    const RealType velocity_ooo = predicted_velocity[node_ooo];
    const RealType velocity_p1o = predicted_velocity[node_p1o];
     

    const RealType grad_m1o = (variable_ooo - variable_m1o ) /
      (dx + dt * (velocity_ooo - velocity_m1o));

    const RealType grad_p1o = (variable_p1o - variable_ooo ) /
      (dx + dt * (velocity_p1o - velocity_ooo));

    const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
    gradient_variable[node_ooo] = limited_grad_variable;    
    
  } 

  //x max
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    index_t ix = nx;
    const int node_ooo = iy * (nx + 1) + ix;
    const int node_m1o = NodeNodeM1O(node_ooo, iy, nx);
    const int node_p1o = iy * (nx + 1) + 1;
      
    //data load
    const RealType variable_m1o = lag_nodal_variable[node_m1o];
    const RealType variable_ooo = lag_nodal_variable[node_ooo];
    const RealType variable_p1o = lag_nodal_variable[node_p1o];

    const RealType velocity_m1o = predicted_velocity[node_m1o];
    const RealType velocity_ooo = predicted_velocity[node_ooo];
    const RealType velocity_p1o = predicted_velocity[node_p1o];
     

    const RealType grad_m1o = (variable_ooo - variable_m1o ) /
      (dx + dt * (velocity_ooo - velocity_m1o));

    const RealType grad_p1o = (variable_p1o - variable_ooo ) /
      (dx + dt * (velocity_p1o - velocity_ooo));

    const RealType limited_grad_variable = LimiterMinmod(grad_m1o,grad_p1o) ;
    gradient_variable[node_ooo] = limited_grad_variable;    
    
  }  

}



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

      out_variable[node_ooo] = out_nodal_variable ;
    }
  }

}



void ProjectNodalIntensiveVariableOrder2XBoundary(index_t nx, 
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

  // x min and x max
   
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect but boundary
    index_t ix = 0;

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + nx;
    const index_t node_m1o = NodeNodeM1O(node_sym, iy, nx);
    const index_t node_p1o = NodeNodeP1O(node_ooo, iy, nx);
      
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

#include "ad_projection_nodal_o2_2d_computation_X.h"

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

#include "ad_projection_nodal_o2_2d_computation_X.h"
  
   
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
	
    index_t node_m1o = nx - 1;
    index_t node_p1o = 1;
 
    index_t cellm1m1 = nx * ny - 1;
    index_t cellp1m1 = (ny - 1) * nx;
    index_t cellm1p1 = nx - 1;
    index_t cellp1p1 = 0;
     
 
    const index_t facexm1m1 = ny * (nx + 1) - 2;
    const index_t facexoom1 = (ny - 1) * (nx + 1);
    const index_t facexm1p1 = nx - 1;
    const index_t facexoop1 = 0;
    const index_t facexp1p1 = 1;
    const index_t facexp1m1 = (ny - 1) * (nx + 1) + 1 ;

#include "ad_projection_nodal_o2_2d_computation_X.h" 

    out_variable[node_cmm] = out_nodal_variable ;
    out_variable[node_cmp] = out_nodal_variable ;
    out_variable[node_cpm] = out_nodal_variable ;
    out_variable[node_cpp] = out_nodal_variable ;
  } // end corners

} // end  ProjectNodalIntensiveVariableOrder2XBoundary



void ReconstructGradientY(index_t nx,
			  index_t ny,
			  RealType dx,
			  RealType dy,
			  const RealType* RESTRICT volume_fluxes,
			  const RealType* RESTRICT lag_variable,
			  RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const int cell_ooo = iy * nx + ix;
      const int cell_om1 = CellCellOM1(cell_ooo, nx);
      const int cell_op1 = CellCellOP1(cell_ooo, nx);
      
      const int face_om2 = CellFaceOM1(cell_om1, iy - 1, nx);
      const int face_om1 = CellFaceOM1(cell_ooo, iy, nx);
      const int face_op1 = CellFaceOP1(cell_ooo, iy, nx);
      const int face_op2 = CellFaceOP1(cell_op1, iy + 1, nx);

      const RealType variable_om1 = lag_variable[cell_om1];
      const RealType variable_ooo = lag_variable[cell_ooo];
      const RealType variable_op1 = lag_variable[cell_op1];
      const RealType d_vol_om2 = volume_fluxes[face_om2];
      const RealType d_vol_om1 = volume_fluxes[face_om1];
      const RealType d_vol_op1 = volume_fluxes[face_op1];
      const RealType d_vol_op2 = volume_fluxes[face_op2];


      const RealType grad_om1 = (variable_ooo - variable_om1 ) /
	(dy + (d_vol_op1 - d_vol_om2) / dx);

      const RealType grad_op1 = (variable_op1 - variable_ooo ) /
	(dy + (d_vol_op2 - d_vol_om1) / dx);

      const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
      gradient_variable[cell_ooo] = limited_grad_variable;
    }
  }

} // end ReconstructGradientY



void ReconstructGradientYBoundary(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT volume_fluxes,
				  const RealType* RESTRICT lag_variable,
				  RealType* RESTRICT gradient_variable) {
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

    const RealType variable_om1 = lag_variable[cell_om1];
    const RealType variable_ooo = lag_variable[cell_ooo];
    const RealType variable_op1 = lag_variable[cell_op1];
    const RealType d_vol_om2 = volume_fluxes[face_om2];
    const RealType d_vol_om1 = volume_fluxes[face_om1];
    const RealType d_vol_op1 = volume_fluxes[face_op1];
    const RealType d_vol_op2 = volume_fluxes[face_op2];


    const RealType grad_om1 = (variable_ooo - variable_om1 ) /
      (dy + (d_vol_op1 - d_vol_om2) / dx);

    const RealType grad_op1 = (variable_op1 - variable_ooo ) /
      (dy + (d_vol_op2 - d_vol_om1) / dx);

    const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
    gradient_variable[cell_ooo] = limited_grad_variable;
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

    const RealType variable_om1 = lag_variable[cell_om1];
    const RealType variable_ooo = lag_variable[cell_ooo];
    const RealType variable_op1 = lag_variable[cell_op1];
    const RealType d_vol_om2 = volume_fluxes[face_om2];
    const RealType d_vol_om1 = volume_fluxes[face_om1];
    const RealType d_vol_op1 = volume_fluxes[face_op1];
    const RealType d_vol_op2 = volume_fluxes[face_op2];


    const RealType grad_om1 = (variable_ooo - variable_om1 ) /
      (dy + (d_vol_op1 - d_vol_om2) / dx);

    const RealType grad_op1 = (variable_op1 - variable_ooo ) /
      (dy + (d_vol_op2 - d_vol_om1) / dx);

    const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
    gradient_variable[cell_ooo] = limited_grad_variable;

  }

} // end ReconstructGradientYBoundary




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
    for (index_t ix = 0; ix < nx; ++ix) {
      
      
      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_density[prev_cell];
      const RealType next_cell_variable = cell_density[next_cell];
      const RealType prev_cell_gradient = cell_density_gradient[prev_cell];
      const RealType next_cell_gradient = cell_density_gradient[next_cell];

      const RealType vol_flux = volume_fluxes[face_ooo];
      const RealType vol_flux_om1 = volume_fluxes[face_om1];
      const RealType vol_flux_op1 = volume_fluxes[face_op1];
      
      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 / dx;
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 / dx;	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
      const RealType mass_flux_ooo =
	(0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
	(0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
      mass_flux[face_ooo] = mass_flux_ooo;

    }
  }

} //end ReconstructMassFluxOrder2Y



void ReconstructMassFluxOrder2YBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_density,
					const RealType* RESTRICT cell_density_gradient,
					RealType* RESTRICT mass_flux) {

  //y min only 
  const index_t iy = 0;
  const index_t iy_sym = ny;
  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect
       
    const index_t face_ooo = (nx * iy) + ix;
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
    const RealType vol_flux_om1 = volume_fluxes[face_om1];
    const RealType vol_flux_op1 = volume_fluxes[face_op1];
      
    const RealType dy_lag_prev_corrected = dy - vol_flux_om1 / dx;
    const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 / dx;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
    const RealType mass_flux_ooo =
      (0.5 * (vol_flux + fabs(vol_flux)) * prev_cell_variable_o2) + 
      (0.5 * (vol_flux - fabs(vol_flux)) * next_cell_variable_o2);
      
    mass_flux[face_ooo] = mass_flux_ooo;
    mass_flux[face_sym] = mass_flux_ooo;

  }
  
} //end  ReconstructMassFluxOrder2YBoundary





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
    for (index_t ix = 0; ix < nx; ++ix) {
      
      const index_t face_ooo = (nx * iy) + ix;

      const index_t prev_cell = FaceYCellM1(face_ooo, iy, nx);
      const index_t next_cell = FaceYCellP1(face_ooo, iy, nx);
      const index_t face_om1 = CellFaceOM1( prev_cell, iy, nx);
      const index_t face_op1 = CellFaceOP1( next_cell, iy, nx);;
      
      const RealType prev_cell_variable = cell_variable[prev_cell];
      const RealType next_cell_variable = cell_variable[next_cell];
      const RealType prev_cell_gradient = cell_variable_gradient[prev_cell];
      const RealType next_cell_gradient = cell_variable_gradient[next_cell];
 
      const RealType mass_flux_face = mass_flux[face_ooo];
      
      const RealType vol_flux_om1 = volume_fluxes[face_om1];
      const RealType vol_flux_op1 = volume_fluxes[face_op1];
      
      const RealType dy_lag_prev_corrected = dy - vol_flux_om1 / dx;
      const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 / dx;	
      const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
      const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
      const RealType variable_flux_ooo =
	(0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
	(0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
      variable_flux[face_ooo] = variable_flux_ooo;
  

    }
  }
} // end ReconstructIntensiveVariableFluxOrder2Y




void ReconstructIntensiveVariableFluxOrder2YBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType* RESTRICT volume_fluxes,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     const RealType* RESTRICT cell_variable_gradient,
						     RealType* RESTRICT variable_flux) {


  index_t iy = 0;
  index_t iy_sym = ny;

  for (index_t ix = 0; ix < nx; ++ix) {//problem not vect but boundary
      
    const index_t face_ooo = (iy * nx) + ix;
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
      
    const RealType vol_flux_om1 = volume_fluxes[face_om1];
    const RealType vol_flux_op1 = volume_fluxes[face_op1];
      
    const RealType dy_lag_prev_corrected = dy - vol_flux_om1 / dx;
    const RealType dy_lag_next_corrected = - dy  - vol_flux_op1 / dx;	
    const RealType prev_cell_variable_o2 =  prev_cell_variable + 0.5 * prev_cell_gradient * dy_lag_prev_corrected;
    const RealType next_cell_variable_o2 =  next_cell_variable + 0.5 * next_cell_gradient * dy_lag_next_corrected;
      
    const RealType variable_flux_ooo =
      (0.5 * (mass_flux_face + fabs(mass_flux_face)) * prev_cell_variable_o2) + 
      (0.5 * (mass_flux_face - fabs(mass_flux_face)) * next_cell_variable_o2);
      
    variable_flux[face_ooo] = variable_flux_ooo;
    variable_flux[face_sym] = variable_flux_ooo;

  }
  
}// end ReconstructIntensiveVariableFluxOrder2YBoundary



void ReconstructGradientNodalY(index_t nx,
			       index_t ny,
			       RealType dy,
			       RealType dt,
			       const RealType* RESTRICT predicted_velocity,
			       const RealType* RESTRICT lag_nodal_variable,
			       RealType* RESTRICT gradient_variable) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx + 1; ++ix) {
    
      const int node_ooo = iy * (nx + 1) + ix;
      const int node_om1 = NodeNodeOM1(node_ooo, iy, nx);
      const int node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      //data load
      const RealType variable_om1 = lag_nodal_variable[node_om1];
      const RealType variable_ooo = lag_nodal_variable[node_ooo];
      const RealType variable_op1 = lag_nodal_variable[node_op1];

      const RealType velocity_om1 = predicted_velocity[node_om1];
      const RealType velocity_ooo = predicted_velocity[node_ooo];
      const RealType velocity_op1 = predicted_velocity[node_op1];
     

      const RealType grad_om1 = (variable_ooo - variable_om1 ) /
       	(dy + dt * (velocity_ooo - velocity_om1));

      const RealType grad_op1 = (variable_op1 - variable_ooo ) /
       	(dy + dt * (velocity_op1 - velocity_ooo));

      const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
      gradient_variable[node_ooo] = limited_grad_variable;
      
    }
  } 

} // end ReconstructGradientNodalY


void ReconstructGradientNodalYBoundary(index_t nx,
				       index_t ny,
				       RealType dy,
				       RealType dt,
				       const RealType* RESTRICT predicted_velocity,
				       const RealType* RESTRICT lag_nodal_variable,
				       RealType* RESTRICT gradient_variable) {
  //ymin
  for (index_t ix = 0; ix < nx + 1; ++ix) {
    index_t iy = 0;

    const int node_ooo = iy * (nx + 1) + ix;
    const int node_om1 = (ny - 1) * (nx + 1) + ix;
    const int node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
    //data load
    const RealType variable_om1 = lag_nodal_variable[node_om1];
    const RealType variable_ooo = lag_nodal_variable[node_ooo];
    const RealType variable_op1 = lag_nodal_variable[node_op1];

    const RealType velocity_om1 = predicted_velocity[node_om1];
    const RealType velocity_ooo = predicted_velocity[node_ooo];
    const RealType velocity_op1 = predicted_velocity[node_op1];
     

    const RealType grad_om1 = (variable_ooo - variable_om1 ) /
      (dy + dt * (velocity_ooo - velocity_om1));

    const RealType grad_op1 = (variable_op1 - variable_ooo ) /
      (dy + dt * (velocity_op1 - velocity_ooo));

    const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
    gradient_variable[node_ooo] = limited_grad_variable;
      
  }

  //ymax
  for (index_t ix = 0; ix < nx + 1; ++ix) {
    index_t iy = ny;
    const int node_ooo = iy * (nx + 1) + ix;
    const int node_om1 = NodeNodeOM1(node_ooo, iy, nx);
    const int node_op1 = ix + nx + 1;;
      
    //data load
    const RealType variable_om1 = lag_nodal_variable[node_om1];
    const RealType variable_ooo = lag_nodal_variable[node_ooo];
    const RealType variable_op1 = lag_nodal_variable[node_op1];

    const RealType velocity_om1 = predicted_velocity[node_om1];
    const RealType velocity_ooo = predicted_velocity[node_ooo];
    const RealType velocity_op1 = predicted_velocity[node_op1];
     

    const RealType grad_om1 = (variable_ooo - variable_om1 ) /
      (dy + dt * (velocity_ooo - velocity_om1));

    const RealType grad_op1 = (variable_op1 - variable_ooo ) /
      (dy + dt * (velocity_op1 - velocity_ooo));

    const RealType limited_grad_variable = LimiterMinmod(grad_om1,grad_op1) ;
    gradient_variable[node_ooo] = limited_grad_variable;
      
  }
} // end ReconstructGradientNodalYBoundary



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

      out_variable[node_ooo] = out_nodal_variable ;
    }
  }

} //end ProjectNodalIntensiveVariableOrder2Y



void ProjectNodalIntensiveVariableOrder2YBoundary(index_t nx, 
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

  // xmin and x max
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {//problem not vect

    index_t ix = 0;
    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + nx;
 
    const index_t node_om1 = NodeNodeOM1(node_ooo, iy, nx);
    const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
    const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
    const index_t cellp1m1 = NodeCellP1M1(node_ooo, iy, nx);
    const index_t cellm1p1 = NodeCellM1P1(node_sym, iy, nx);
    const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

    const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
    const index_t faceym1oo = NodeFaceYM1OO(node_sym, iy, nx);
    const index_t faceym1p1 = NodeFaceYM1P1(node_sym, iy, nx);
    const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
    const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
    const index_t faceyp1m1 = NodeFaceYP1M1(node_ooo, iy, nx);

#include "ad_projection_nodal_o2_2d_computation_Y.h"
    out_variable[node_ooo] = out_nodal_variable ;
    out_variable[node_sym] = out_nodal_variable;
  }  // end X boundary


  //ymin and y max
  for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {//problem not vect
    index_t iy = 0;
    index_t iy_sym = ny;

    const index_t node_ooo = ((nx + 1) * iy) + ix;
    const index_t node_sym = node_ooo + ny * (nx + 1);
 
    const index_t node_om1 = NodeNodeOM1(node_sym, iy_sym, nx);
    const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
    const index_t cellm1m1 = NodeCellM1M1(node_sym, iy_sym, nx);
    const index_t cellp1m1 = NodeCellP1M1(node_sym, iy_sym, nx);
    const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
    const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

    const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy_sym, nx);
    const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
    const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
    const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
    const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
    const index_t faceyp1m1 = NodeFaceYP1M1(node_sym, iy_sym, nx);

#include "ad_projection_nodal_o2_2d_computation_Y.h"
    
    out_variable[node_ooo] = out_nodal_variable ;
    out_variable[node_sym] = out_nodal_variable;
  }//end  Y boundary


  //corners
  {
    const index_t node_cmm = 0;
    const index_t node_cpm = nx;
    const index_t node_cmp = ny * (nx + 1) ;
    const index_t node_cpp = (nx + 1) * (ny + 1) - 1;
    const index_t node_ooo = node_cmm;

    const index_t node_om1 = (ny - 1) * (nx + 1) ;
    const index_t node_op1 = nx + 1;
      
    const index_t cellm1m1 = nx * ny - 1 ;
    const index_t cellp1m1 = (ny - 1) * nx;
    const index_t cellm1p1 = nx - 1;
    const index_t cellp1p1 = 0;

    const index_t faceym1m1 = ny * nx - 1;
    const index_t faceym1oo = nx - 1;
    const index_t faceym1p1 = 2 * nx - 1;
    const index_t faceyp1oo = 0;
    const index_t faceyp1p1 = nx;
    const index_t faceyp1m1 = (ny - 1) * nx;	
   

#include "ad_projection_nodal_o2_2d_computation_Y.h" 

    out_variable[node_cmm] = out_nodal_variable ;
    out_variable[node_cmp] = out_nodal_variable ;
    out_variable[node_cpm] = out_nodal_variable ;
    out_variable[node_cpp] = out_nodal_variable ;
  } // end corners

}// end ProjectNodalIntensiveVariableOrder2YBoundar

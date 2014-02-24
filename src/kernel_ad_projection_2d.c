// -*- c++ -*- (for emacs users)
#include "kernel_ad_projection_2d.h"

//#define NDEBUG
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

  //#pragma omp parallel for
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

      // printf("INFO: nx=%d, ny=%d, iy=%d, prev_node=%d, next_node=%d, vx_prev=%lf, vx_next=%lf, volume_flux=%lf\n", 
      // 	     nx, ny, iy, prev_node, next_node, vx_prev, vx_next, volume_flux);

    }
  }  

  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;

      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

      const RealType volume_flux_prev = volume_fluxes[prev_face];
      const RealType volume_flux_next = volume_fluxes[next_face];
      const RealType current_mass = cell_mass[cell_ooo];

      const RealType current_directional_lagrangian_volume =
	(dx * dy) + volume_flux_prev - volume_flux_next;

      const RealType current_directional_lagrangian_density =
	current_mass / current_directional_lagrangian_volume;

      directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume;
      
      directional_lagrangian_density[cell_ooo] = current_directional_lagrangian_density;

      //printf("%lf\n", current_directional_lagrangian_density);

      assert(0.0 < current_directional_lagrangian_volume);
      assert(0.0 < current_directional_lagrangian_density);
      
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

  //#pragma omp parallel for
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

      // printf("INFO: nx=%d, ny=%d, iy=%d, prev_node=%d, next_node=%d, vx_prev=%lf, vx_next=%lf, volume_flux=%lf\n", 
      // 	     nx, ny, iy, prev_node, next_node, vx_prev, vx_next, volume_flux);

    }
  }  

  //#pragma omp parallel for
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

void CheckFluxPeriodicalPropertyX(index_t nx,
				  index_t ny,
				  const RealType* RESTRICT flux) {

  for (index_t iy = 0; iy < ny; ++iy) {
    
    const index_t first_face = (nx + 1) * iy + 0;
    const index_t lasst_face = (nx + 1) * iy + nx;

    const RealType first_flux = flux[first_face];
    const RealType lasst_flux = flux[lasst_face];
    
    if (first_flux != lasst_flux) {

      printf("ERROR: nx=%d, ny=%d, iy=%d, first_flux=%lf, lasst_flux=%lf\n", nx, ny, iy, first_flux, lasst_flux);
      assert(flux[first_face] == flux[lasst_face]);
    
    }
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

      printf("ERROR: nx=%d, ny=%d, ix=%d, first_flux=%lf, lasst_flux=%lf\n", nx, ny, ix, first_flux, lasst_flux);
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

  //#pragma omp parallel for
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

  //#pragma omp parallel for
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

  //#pragma omp parallel for
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

  //#pragma omp parallel for
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

  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    
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

  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    
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
				   RealType* RESTRICT out_cell_variable) {
  
  //#pragma omp parallel for
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
      
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo];
      
      const RealType out_cell_variable_ooo = 
	in_cell_mass_ooo * in_cell_variable_ooo +
	//	mass_flux_prev * face_variable_prev -
	//        mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      out_cell_variable[cell_ooo] = out_cell_variable_ooo; 
      
    }

  }  

}

void MassProjectIntensiveVariableY(index_t nx, 
				   index_t ny, 
				   const RealType* RESTRICT in_cell_mass,
				   //const RealType* RESTRICT mass_flux,
				   const RealType* RESTRICT in_cell_variable,
				   const RealType* RESTRICT in_variable_flux,
				   RealType* RESTRICT out_cell_variable) {

  //#pragma omp parallel for
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
      
      const RealType in_cell_variable_ooo = in_cell_variable[cell_ooo];
      
      const RealType out_cell_variable_ooo = in_cell_mass_ooo * in_cell_variable_ooo +
	//mass_flux_prev * face_variable_prev -
        //mass_flux_next * face_variable_next;
	face_variable_prev -
        face_variable_next;

      out_cell_variable[cell_ooo] = out_cell_variable_ooo; 
      
    }

  }  

}


void ProjectNodalIntensiveVariableX(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT in_cell_mass,
				    const RealType* RESTRICT in_vx,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_moment) {

  //Boucle sur les noeuds
  //#pragma omp parallel for
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

      const RealType node_mass = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      
      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = node_mass * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;

      out_moment[node_ooo] = out_moment_ooo; 
 

     // printf("INFO: nx=%d, ny=%d, iy=%d, prev_node=%d, next_node=%d, vx_prev=%lf, vx_next=%lf, volume_flux=%lf\n", 
      // 	     nx, ny, iy, prev_node, next_node, vx_prev, vx_next, volume_flux);

    }
  }  

}

void ProjectNodalIntensiveVariableY(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT in_cell_mass,
				    const RealType* RESTRICT in_vx,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_moment) {

  //Boucle sur les noeuds
  //#pragma omp parallel for
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

      const RealType node_mass = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      
      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = node_mass * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;

      out_moment[node_ooo] = out_moment_ooo; 
 
    }
  }  

}

void ProjectNodalIntensiveVariableBoundaryX(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT in_cell_mass,
				    const RealType* RESTRICT in_vx,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_moment) {

  // xmin (xmax is deduced in periodic conditions)

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t iy = halo_width; iy < ny + 1 - halo_width; ++iy) {
    for (index_t ix = 0; ix < halo_width; ++ix) {

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

      const RealType node_mass = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      
      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = node_mass * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;

      out_moment[node_ooo] = out_moment_ooo; 

      out_moment[node_sym] = out_moment_ooo;
 
    }
  }  

}

void ProjectNodalIntensiveVariableBoundaryY(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* RESTRICT in_cell_mass,
				    const RealType* RESTRICT in_vx,
				    const RealType* RESTRICT mass_flux,
				    RealType* RESTRICT out_moment) {

  // ymin (ymax is deduced in periodic conditions)

  //Boucle sur les noeuds
  //#pragma omp parallel for
  for (index_t iy = 0; iy < halo_width; ++iy) {
    for (index_t ix = halo_width; ix < nx + 1 - halo_width; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;
      const index_t node_sym = node_ooo + (nx + 1) * ny;

      const index_t node_om1 = NodeNodeOM1(node_sym, iy, nx);
      const index_t node_op1 = NodeNodeOP1(node_ooo, iy, nx);
      
      const index_t cellm1m1 = NodeCellM1M1(node_sym, iy, nx);
      const index_t cellp1m1 = NodeCellP1M1(node_sym, iy, nx);
      const index_t cellm1p1 = NodeCellM1P1(node_ooo, iy, nx);
      const index_t cellp1p1 = NodeCellP1P1(node_ooo, iy, nx);

      const index_t faceym1m1 = NodeFaceYM1M1(node_sym, iy, nx);
      const index_t faceym1oo = NodeFaceYM1OO(node_ooo, iy, nx);
      const index_t faceym1p1 = NodeFaceYM1P1(node_ooo, iy, nx);
      const index_t faceyp1p1 = NodeFaceYP1P1(node_ooo, iy, nx);
      const index_t faceyp1oo = NodeFaceYP1OO(node_ooo, iy, nx);
      const index_t faceyp1m1 = NodeFaceYP1M1(node_sym, iy, nx);

#include "ad_nodal_projection_2d_Y_data_load.h"

      const RealType node_mass = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
      
      const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_p1m1 + mass_flux_m1oo + mass_flux_p1oo);
      
      const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1oo + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_m1p1);  
      
#include "reconstruct_dual_variable_computation.h"      

      const RealType in_moment_ooo = node_mass * in_vx_oo; 

      const RealType out_moment_ooo = 
	in_moment_ooo + moment_flux_prev - moment_flux_next;

      out_moment[node_ooo] = out_moment_ooo; 
 
      out_moment[node_sym] = out_moment_ooo;
    }
  }  

}

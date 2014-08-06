#ifndef KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H

#include "kernel_multimat_interface_2d.h"

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif


void ReconstructGradientMultimatXPeriodicBoundary(index_t nx,
						  index_t ny,
						  RealType dx,
						  RealType dy,
						  const RealType* in_cell_volumic_fraction,
						  const RealType* volume_fluxes_1,
						  const RealType* volume_fluxes_2,
						  const RealType* lag_variable_1,
						  const RealType* lag_variable_2,
						  RealType* gradient_variable_1,
						  RealType* gradient_variable_2);


void ReconstructGradientMultimatYPeriodicBoundary(index_t nx,
						  index_t ny,
						  RealType dx,
						  RealType dy,
						  const RealType* in_cell_volumic_fraction,
						  const RealType* volume_fluxes_1,
						  const RealType* volume_fluxes_2,
						  const RealType* lag_variable_1,
						  const RealType* lag_variable_2,
						  RealType* gradient_variable_1,
						  RealType* gradient_variable_2);


void ReconstructNormalInterfacePeriodicBoundary(// in
						index_t nx,
						index_t ny,
						RealType dx,
						RealType dy,
						const RealType* in_variable,
						// out
						RealType* interface_normal_x,
						RealType* interface_normal_y);



void ReconstructNormalInterfaceBoundaryRt(// in
					  index_t nx,
					  index_t ny,
					  RealType dx,
					  RealType dy,
					  const RealType* in_variable,
					  // out
					  RealType* interface_normal_x,
					  RealType* interface_normal_y);



void ReconstructNormalInterfaceBoundaryWall(// in
					    index_t nx,
					    index_t ny,
					    RealType dx,
					    RealType dy,
					    const RealType* in_variable,
					    // out
					    RealType* interface_normal_x,
					    RealType* interface_normal_y);


void ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundaryWall(// in
								     index_t nx, 
								     index_t ny, 
								     RealType dx,
								     RealType dy,
								     const RealType* volume_fluxes,
								     const RealType* directional_lagrangian_volume,
								     const RealType* cell_volumic_fraction,
								     const RealType* interface_normal_x,
								     const RealType* interface_normal_y,
								     // out
								     RealType* volume_fluxes_1,
								     RealType* volume_fluxes_2);


void ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundaryWall(// in
								     index_t nx, 
								     index_t ny, 
								     RealType dx,
								     RealType dy,
								     const RealType* volume_fluxes_y,
								     const RealType* directional_lagrangian_volume,
								     const RealType* cell_volumic_fraction,
								     const RealType* interface_normal_x,
								     const RealType* interface_normal_y,
								     // out
								     RealType* volume_fluxes_1_y,
								     RealType* volume_fluxes_2_y);


void ComputeDirectionalLagrangianFractionalVolumeFluxesXPeriodicBoundary(// in
									 index_t nx, 
									 index_t ny, 
									 RealType dx,
									 RealType dy,
									 const RealType* volume_fluxes,
									 const RealType* directional_lagrangian_volume,
									 const RealType* cell_volumic_fraction,
									 const RealType* interface_normal_x,
									 const RealType* interface_normal_y,
									 // out
									 RealType* volume_fluxes_1,
									 RealType* volume_fluxes_2);


void ComputeDirectionalLagrangianFractionalVolumeFluxesYPeriodicBoundary(// in
									 index_t nx, 
									 index_t ny, 
									 RealType dx,
									 RealType dy,
									 const RealType* volume_fluxes_y,
									 const RealType* directional_lagrangian_volume,
									 const RealType* cell_volumic_fraction,
									 const RealType* interface_normal_x,
									 const RealType* interface_normal_y,
									 // out
									 RealType* volume_fluxes_1_y,
									 RealType* volume_fluxes_2_y);



#ifdef __cplusplus
}
#endif


#endif // KERNEL_PROJECTION_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H

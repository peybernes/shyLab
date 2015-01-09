#ifndef KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H

#include "kernel_multimat_interface_2d.h"

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif


  // Periodic

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


void ReconstructGradientMultimatDiagPeriodicBoundary(index_t nx,
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* predicted_u,
						     const RealType* predicted_v,
						     const RealType* in_cell_volumic_fraction,
						     const RealType* volume_fluxes_1_corner,
						     const RealType* volume_fluxes_2_corner,
						     const RealType* lag_variable_1,
						     const RealType* lag_variable_2,
						     RealType* gradient_variable_1_diag,
						     RealType* gradient_variable_2_diag,
						     RealType* gradient_variable_1_antidiag,
						     RealType* gradient_variable_2_antidiag);


void ReconstructNormalInterfacePeriodicBoundary(// in
						index_t nx,
						index_t ny,
						RealType dx,
						RealType dy,
						const RealType* in_variable,
						// out
						RealType* interface_normal_x,
						RealType* interface_normal_y);


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

void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerPeriodicBoundary(// in
										index_t nx, 
										index_t ny, 
										RealType dx,
										RealType dy,
										RealType dt,
										const RealType* volume_fluxes,
										const RealType* volume_fluxes_y,
										const RealType* directional_lagrangian_volume,
										const RealType* cell_volumic_fraction,
										const RealType* interface_normal_x,
										const RealType* interface_normal_y,
										const RealType* in_vx,
										const RealType* in_vy,
										// out
										RealType* volume_fluxes_1,
										RealType* volume_fluxes_1_y,
										RealType* volume_fluxes_1_corner,
										RealType* volume_fluxes_2,
										RealType* volume_fluxes_2_y,
										RealType* volume_fluxes_2_corner);

  // Rayleigh Taylor

void ReconstructNormalInterfaceBoundaryRt(// in
					  index_t nx,
					  index_t ny,
					  RealType dx,
					  RealType dy,
					  const RealType* in_variable,
					  // out
					  RealType* interface_normal_x,
					  RealType* interface_normal_y);



  // Wall

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

void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerBoundaryWall(// in
									    index_t nx, 
									    index_t ny, 
									    RealType dx,
									    RealType dy,
									    RealType dt,
									    const RealType* volume_fluxes,
									    const RealType* volume_fluxes_y,
									    const RealType* directional_lagrangian_volume,
									    const RealType* cell_volumic_fraction,
									    const RealType* interface_normal_x,
									    const RealType* interface_normal_y,
									    const RealType* in_vx,
									    const RealType* in_vy,
									    // out
									    RealType* volume_fluxes_1,
									    RealType* volume_fluxes_1_y,
									    RealType* volume_fluxes_1_corner,
									    RealType* volume_fluxes_2,
									    RealType* volume_fluxes_2_y,
									    RealType* volume_fluxes_2_corner);

#ifdef __cplusplus
}
#endif


#endif // KERNEL_PROJECTION_MULTIMAT_INTERFACE_2D_BOUNDARY_CONDITIONS_H

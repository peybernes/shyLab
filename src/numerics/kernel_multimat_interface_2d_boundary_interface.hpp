#ifndef KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_INTERFACE_H
#define KERNEL_MULTIMAT_INTERFACE_2D_BOUNDARY_INTERFACE_H
#include <string>
typedef int index_t;

#include "kernel_multimat_interface_2d.h"
#include "kernel_multimat_interface_2d_boundary_conditions.h"

void ReconstructGradientMultimatXBoundary(std::string BoundaryConditions,
					  index_t nx,
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

void ReconstructGradientMultimatYBoundary(std::string BoundaryConditions,
					  index_t nx,
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


void ReconstructNormalInterfaceBoundary(std::string BoundaryConditions,
					index_t nx,
					index_t ny,
					RealType dx,
					RealType dy,
					const RealType* in_cell_volumic_fraction,
					// out
					RealType* interface_normal_x,
					RealType* interface_normal_y);


void ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 RealType dx,
								 RealType dy,
								 const RealType* volume_fluxes_x,
								 const RealType* directional_lagrangian_volume,
								 const RealType* in_cell_volumic_fraction,
								 const RealType* interface_normal_x,
								 const RealType* interface_normal_y,
								 // out
								 RealType* volume_fluxes_1_x,
								 RealType* volume_fluxes_2_x);


void ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 RealType dx,
								 RealType dy,
								 const RealType* volume_fluxes_y,
								 const RealType* directional_lagrangian_volume_y,
								 const RealType* in_cell_volumic_fraction,
								 const RealType* interface_normal_x,
								 const RealType* interface_normal_y,
								 // out
								 RealType* volume_fluxes_1_y,
								 RealType* volume_fluxes_2_y);


#endif 

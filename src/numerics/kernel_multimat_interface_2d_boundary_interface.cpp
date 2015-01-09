// -*- c++ -*- (for emacs users)
#include "kernel_multimat_interface_2d_boundary_interface.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"


void ReconstructGradientMultimatXBoundary(std::string BoundaryConditions,
					  index_t nx,
					  index_t ny,
					  RealType dx,
					  RealType dy,
					  const RealType* RESTRICT in_cell_volumic_fraction,
					  const RealType* RESTRICT volume_fluxes_1,
					  const RealType* RESTRICT volume_fluxes_2,
					  const RealType* RESTRICT lag_variable_1,
					  const RealType* RESTRICT lag_variable_2,
					  RealType* RESTRICT gradient_variable_1,
					  RealType* RESTRICT gradient_variable_2) {

  if (BoundaryConditions == "periodic") {
    ReconstructGradientMultimatXPeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1,volume_fluxes_2, lag_variable_1, lag_variable_2, gradient_variable_1, gradient_variable_2);
  }
  else if (BoundaryConditions == "Wall") {
  }
  else{
    ReconstructGradientMultimatXPeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1,volume_fluxes_2, lag_variable_1, lag_variable_2, gradient_variable_1, gradient_variable_2);
  }
}

void ReconstructGradientMultimatYBoundary(std::string BoundaryConditions,
					  index_t nx,
					  index_t ny,
					  RealType dx,
					  RealType dy,
					  const RealType* RESTRICT in_cell_volumic_fraction,
					  const RealType* RESTRICT volume_fluxes_1,
					  const RealType* RESTRICT volume_fluxes_2,
					  const RealType* RESTRICT lag_variable_1,
					  const RealType* RESTRICT lag_variable_2,
					  RealType* RESTRICT gradient_variable_1,
					  RealType* RESTRICT gradient_variable_2) {

  if (BoundaryConditions == "periodic") {
    ReconstructGradientMultimatYPeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1,volume_fluxes_2, lag_variable_1, lag_variable_2, gradient_variable_1, gradient_variable_2);
  }
  else if (BoundaryConditions == "Wall") {
  }
  else{
    ReconstructGradientMultimatYPeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1,volume_fluxes_2, lag_variable_1, lag_variable_2, gradient_variable_1, gradient_variable_2);
  }
}


void ReconstructGradientMultimatDiagBoundary(std::string BoundaryConditions,
					     index_t nx,
					     index_t ny,
					     RealType dx,
					     RealType dy,
					     RealType dt,
					     const RealType* RESTRICT predicted_u,
					     const RealType* RESTRICT predicted_v,
					     const RealType* RESTRICT in_cell_volumic_fraction,
					     const RealType* RESTRICT volume_fluxes_1_corner,
					     const RealType* RESTRICT volume_fluxes_2_corner,
					     const RealType* RESTRICT lag_variable_1,
					     const RealType* RESTRICT lag_variable_2,
					     RealType* RESTRICT gradient_variable_1_diag,
					     RealType* RESTRICT gradient_variable_2_diag,
					     RealType* RESTRICT gradient_variable_1_antidiag,
					     RealType* RESTRICT gradient_variable_2_antidiag) {

  if (BoundaryConditions == "periodic") {

    ReconstructGradientMultimatDiagPeriodicBoundary(nx, ny, dx, dy, dt,predicted_u, predicted_v, in_cell_volumic_fraction, volume_fluxes_1_corner, volume_fluxes_2_corner, lag_variable_1, lag_variable_2, gradient_variable_1_diag, gradient_variable_2_diag, gradient_variable_1_antidiag, gradient_variable_2_antidiag);

  } else if (BoundaryConditions == "Wall") {

  } else {

    ReconstructGradientMultimatDiagPeriodicBoundary(nx, ny, dx, dy, dt,predicted_u, predicted_v, in_cell_volumic_fraction, volume_fluxes_1_corner, volume_fluxes_2_corner, lag_variable_1, lag_variable_2, gradient_variable_1_diag, gradient_variable_2_diag, gradient_variable_1_antidiag, gradient_variable_2_antidiag);

  }
}


void ReconstructNormalInterfaceBoundary(std::string BoundaryConditions,
					index_t nx,
					index_t ny,
					RealType dx,
					RealType dy,
					const RealType* RESTRICT in_cell_volumic_fraction,
					// out
					RealType* RESTRICT interface_normal_x,
					RealType* RESTRICT interface_normal_y) {

  
  if (BoundaryConditions == "periodic") {
    ReconstructNormalInterfacePeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);

  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructNormalInterfaceBoundaryRt(nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);  
  
  } 
  else if (BoundaryConditions == "Wall") {
    ReconstructNormalInterfaceBoundaryWall(nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);  
  
  }  
  else{
    ReconstructNormalInterfacePeriodicBoundary(nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);

  }

}



void ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 RealType dx,
								 RealType dy,
								 const RealType* RESTRICT volume_fluxes_x,
								 const RealType* RESTRICT directional_lagrangian_volume,
								 const RealType* RESTRICT in_cell_volumic_fraction,
								 const RealType* RESTRICT interface_normal_x,
								 const RealType* RESTRICT interface_normal_y,
								 // out
								 RealType* RESTRICT volume_fluxes_1_x,
								 RealType* RESTRICT volume_fluxes_2_x) {
  

  if (BoundaryConditions == "periodic") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesXPeriodicBoundary(nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_volume, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_x, volume_fluxes_2_x);

  } 
  else if (BoundaryConditions == "Wall") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundaryWall(nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_volume, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_x, volume_fluxes_2_x);

  }  
  else{
    ComputeDirectionalLagrangianFractionalVolumeFluxesXPeriodicBoundary(nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_volume, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_x, volume_fluxes_2_x);

  }

}



void ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 RealType dx,
								 RealType dy,
								 const RealType* RESTRICT volume_fluxes_y,
								 const RealType* RESTRICT directional_lagrangian_volume_y,
								 const RealType* RESTRICT in_cell_volumic_fraction,
								 const RealType* RESTRICT interface_normal_x,
								 const RealType* RESTRICT interface_normal_y,
								 // out
								 RealType* RESTRICT volume_fluxes_1_y,
								 RealType* RESTRICT volume_fluxes_2_y) {
  

  if (BoundaryConditions == "periodic") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesYPeriodicBoundary(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_volume_y, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_y, volume_fluxes_2_y);

  } 
  else if (BoundaryConditions == "Wall") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundaryWall(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_volume_y, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_y, volume_fluxes_2_y);

  }  
  else {
    ComputeDirectionalLagrangianFractionalVolumeFluxesYPeriodicBoundary(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_volume_y, in_cell_volumic_fraction, interface_normal_x, interface_normal_y, volume_fluxes_1_y, volume_fluxes_2_y);

  }

}


void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerBoundary(// in
									std::string BoundaryConditions,
									index_t nx, 
									index_t ny, 
									RealType dx,
									RealType dy,
									RealType dt,
									const RealType* RESTRICT volume_fluxes,
									const RealType* RESTRICT volume_fluxes_y,
									const RealType* RESTRICT directional_lagrangian_volume,
									const RealType* RESTRICT cell_volumic_fraction,
									const RealType* RESTRICT interface_normal_x,
									const RealType* RESTRICT interface_normal_y,
									const RealType* RESTRICT in_vx,
									const RealType* RESTRICT in_vy,
									// out
									RealType* RESTRICT volume_fluxes_1,
									RealType* RESTRICT volume_fluxes_1_y,
									RealType* RESTRICT volume_fluxes_1_corner,
									RealType* RESTRICT volume_fluxes_2,
									RealType* RESTRICT volume_fluxes_2_y,
									RealType* RESTRICT volume_fluxes_2_corner) {

  if (BoundaryConditions == "periodic") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerPeriodicBoundary(nx, ny, dx, dy, dt, volume_fluxes, volume_fluxes_y, directional_lagrangian_volume, cell_volumic_fraction, interface_normal_x, interface_normal_y, in_vx, in_vy, volume_fluxes_1, volume_fluxes_1_y, volume_fluxes_1_corner, volume_fluxes_2, volume_fluxes_2_y, volume_fluxes_2_corner);

  } 
  else if (BoundaryConditions == "Wall") {
    ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerBoundaryWall(nx, ny, dx, dy, dt, volume_fluxes, volume_fluxes_y, directional_lagrangian_volume, cell_volumic_fraction, interface_normal_x, interface_normal_y, in_vx, in_vy, volume_fluxes_1, volume_fluxes_1_y, volume_fluxes_1_corner, volume_fluxes_2, volume_fluxes_2_y, volume_fluxes_2_corner);

  }  
  else {
    ComputeDirectionalLagrangianFractionalVolumeFluxesXYCornerPeriodicBoundary(nx, ny, dx, dy, dt, volume_fluxes, volume_fluxes_y, directional_lagrangian_volume, cell_volumic_fraction, interface_normal_x, interface_normal_y, in_vx, in_vy, volume_fluxes_1, volume_fluxes_1_y, volume_fluxes_1_corner, volume_fluxes_2, volume_fluxes_2_y, volume_fluxes_2_corner);

  }

}

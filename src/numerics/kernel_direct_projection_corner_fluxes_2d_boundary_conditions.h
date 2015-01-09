#ifndef KERNEL_DIRECT_PROJECTION_CORNER_FLUXES_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_DIRECT_PROJECTION_CORNER_FLUXES_2D_BOUNDARY_CONDITIONS_H

typedef int index_t;

#include"kernels.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <string.h>


  //Periodic boundary

void ReconstructGradientDiagPeriodicBoundary(index_t nx,
					     index_t ny,
					     RealType dx,
					     RealType dy,
					     RealType dt,
					     const RealType* predicted_u,
					     const RealType* predicted_v,
					     const RealType* lag_variable,
					     RealType* gradient_variable_diag,
					     RealType* gradient_variable_antidiag);


void ReconstructMassFluxOrder2XCornerFluxesPeriodicBoundary(index_t nx, 
							    index_t ny, 
							    index_t halo_width,
							    const RealType dx,
							    const RealType dy,
							    const RealType dt,
							    const RealType* predicted_u,
							    const RealType* predicted_v,
							    const RealType* volume_fluxes_tot,
							    const RealType* volume_fluxes,
							    const RealType* cell_density,
							    const RealType* cell_density_gradient,
							    RealType* reconstructed_density_faces_x,
							    RealType* mass_flux);


void ReconstructMassFluxOrder2YCornerFluxesPeriodicBoundary(index_t nx, 
							    index_t ny, 
							    index_t halo_width,
							    const RealType dx,
							    const RealType dy,
							    const RealType dt,
							    const RealType* predicted_u,
							    const RealType* predicted_v,
							    const RealType* volume_fluxes_tot,
							    const RealType* volume_fluxes,
							    const RealType* cell_density,
							    const RealType* cell_density_gradient,
							    RealType* reconstructed_density_faces_y,
							    RealType* mass_flux);

void ReconstructMassFluxOrder2CornerPeriodicBoundary(index_t nx, 
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* directional_lagrangian_density,
						     const RealType* reconstructed_density_faces_x,
						     const RealType* reconstructed_density_faces_y,
						     const RealType* volume_fluxes_corner, 
						     const RealType* in_vx, 
						     const RealType* in_vy,
						     RealType* mass_corner_fluxes); 


void ReconstructMassFluxOrder2DiagCornerPeriodicBoundary(index_t nx, 
							 index_t ny,
							 RealType dx,
							 RealType dy,
							 RealType dt,
							 const RealType* predicted_u,
							 const RealType* predicted_v,
							 const RealType* directional_lagrangian_density,
							 const RealType* volume_fluxes_corner, 
							 const RealType* density_gradient_diag, 
							 const RealType* density_gradient_antidiag, 
							 RealType* mass_corner_fluxes);


void ProjectMassCornerFluxesPeriodicBoundary(index_t nx, 
					     index_t ny,
					     RealType dt, 
					     const RealType* mass_corner_fluxes,
					     const RealType* in_vx, 
					     const RealType* in_vy, 
					     RealType* out_cell_mass);


void ReconstructIntensiveVariableFluxOrder2XCornerFluxesPeriodicBoundary(index_t nx, 
									 index_t ny, 
									 index_t halo_width,
									 const RealType dx,
									 const RealType dy,
									 const RealType dt,
									 const RealType* predicted_u,
									 const RealType* predicted_v,
									 const RealType* volume_fluxes,
									 const RealType* mass_flux,
									 const RealType* cell_variable,
									 const RealType* cell_variable_gradient,
									 RealType* reconstructed_variable_faces_x,
									 RealType* variable_flux);


void ReconstructIntensiveVariableFluxOrder2YCornerFluxesPeriodicBoundary(index_t nx, 
									 index_t ny, 
									 index_t halo_width,
									 const RealType dx,
									 const RealType dy,
									 const RealType dt,
									 const RealType* predicted_u,
									 const RealType* predicted_v,
									 const RealType* volume_fluxes,
									 const RealType* mass_flux,
									 const RealType* cell_variable,
									 const RealType* cell_variable_gradient,
									 RealType* reconstructed_variable_faces_y,
									 RealType* variable_flux);


void ReconstructIntensiveVariableFluxOrder2CornerPeriodicBoundary(index_t nx, 
								  index_t ny,
								  RealType dx,
								  RealType dy,
								  RealType dt,
								  const RealType* cell_variable,
								  const RealType* reconstructed_variable_faces_x,
								  const RealType* reconstructed_variable_faces_y,
								  const RealType* mass_corner_fluxes, 
								  const RealType* in_vx, 
								  const RealType* in_vy,
								  RealType* variable_flux_corner);



void ReconstructIntensiveVariableFluxOrder2DiagCornerPeriodicBoundary(index_t nx, 
								      index_t ny,
								      RealType dx,
								      RealType dy,
								      RealType dt,
								      const RealType* predicted_u,
								      const RealType* predicted_v,
								      const RealType* cell_variable,
								      const RealType* mass_corner_fluxes, 
								      const RealType* variable_gradient_diag, 
								      const RealType* variable_gradient_antidiag, 
								      RealType* variable_flux_corner);

void MassProjectIntensiveVariableCornerFluxesPeriodicBoundary(index_t nx, 
							      index_t ny,
							      RealType dt,
							      const RealType* in_cell_variable,
							      const RealType* variable_flux_corner,
							      const RealType* in_vx, 
							      const RealType* in_vy, 
							      const RealType* out_cell_mass,
							      RealType* out_cell_variable);


void ProjectNodalIntensiveVariableCornerFluxesPeriodicBoundary(index_t nx, 
							       index_t ny, 
							       RealType dt,
							       index_t halo_width,
							       const RealType* out_cell_mass,
							       const RealType* in_variable,
							       const RealType* in_u,
							       const RealType* in_v,
							       const RealType* mass_corner_fluxes,
							       const index_t* sign_x_corner_fluxes,
							       const index_t* sign_y_corner_fluxes,
							       RealType* out_variable);



  // Boundary Wall

void ReconstructGradientDiagBoundaryWall(index_t nx,
					 index_t ny,
					 RealType dx,
					 RealType dy,
					 RealType dt,
					 const RealType* predicted_u,
					 const RealType* predicted_v,
					 const RealType* lag_variable,
					 RealType* gradient_variable_diag,
					 RealType* gradient_variable_antidiag);

void ReconstructMassFluxOrder2XCornerFluxesBoundaryWall(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dx,
							const RealType dy,
							const RealType dt,
							const RealType* predicted_u,
							const RealType* predicted_v,
							const RealType* volume_fluxes_tot,
							const RealType* volume_fluxes,
							const RealType* cell_density,
							const RealType* cell_density_gradient,
							RealType* reconstructed_density_faces_x,
							RealType* mass_flux);


void ReconstructMassFluxOrder2YCornerFluxesBoundaryWall(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dx,
							const RealType dy,
							const RealType dt,
							const RealType* predicted_u,
							const RealType* predicted_v,
							const RealType* volume_fluxes_tot,
							const RealType* volume_fluxes,
							const RealType* cell_density,
							const RealType* cell_density_gradient,
							RealType* reconstructed_density_faces_y,
							RealType* mass_flux);


void ReconstructMassFluxOrder2CornerBoundaryWall(index_t nx, 
						 index_t ny,
						 RealType dx,
						 RealType dy,
						 RealType dt,
						 const RealType* directional_lagrangian_density,
						 const RealType* reconstructed_density_faces_x,
						 const RealType* reconstructed_density_faces_y,
						 const RealType* volume_fluxes_corner, 
						 const RealType* in_vx, 
						 const RealType* in_vy,
						 RealType* mass_corner_fluxes);

 
void ReconstructMassFluxOrder2DiagCornerBoundaryWall(index_t nx, 
						     index_t ny,
						     RealType dx,
						     RealType dy,
						     RealType dt,
						     const RealType* predicted_u,
						     const RealType* predicted_v,
						     const RealType* directional_lagrangian_density,
						     const RealType* volume_fluxes_corner, 
						     const RealType* density_gradient_diag, 
						     const RealType* density_gradient_antidiag, 
						     RealType* mass_corner_fluxes);


void ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundaryWall(index_t nx, 
								     index_t ny, 
								     index_t halo_width,
								     const RealType dx,
								     const RealType dy,
								     const RealType dt,
								     const RealType* predicted_u,
								     const RealType* predicted_v,
								     const RealType* volume_fluxes,
								     const RealType* mass_flux,
								     const RealType* cell_variable,
								     const RealType* cell_variable_gradient,
								     RealType* reconstructed_variable_faces_x,
								     RealType* variable_flux);


void ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundaryWall(index_t nx, 
								     index_t ny, 
								     index_t halo_width,
								     const RealType dx,
								     const RealType dy,
								     const RealType dt,
								     const RealType* predicted_u,
								     const RealType* predicted_v,
								     const RealType* volume_fluxes,
								     const RealType* mass_flux,
								     const RealType* cell_variable,
								     const RealType* cell_variable_gradient,
								     RealType* reconstructed_variable_faces_y,
								     RealType* variable_flux);


void ReconstructIntensiveVariableFluxOrder2CornerBoundaryWall(index_t nx, 
							      index_t ny,
							      RealType dx,
							      RealType dy,
							      RealType dt,
							      const RealType* cell_variable,
							      const RealType* reconstructed_variable_faces_x,
							      const RealType* reconstructed_variable_faces_y,
							      const RealType* mass_corner_fluxes, 
							      const RealType* in_vx, 
							      const RealType* in_vy,
							      RealType* variable_flux_corner);


void ReconstructIntensiveVariableFluxOrder2DiagCornerBoundaryWall(index_t nx, 
								  index_t ny,
								  RealType dx,
								  RealType dy,
								  RealType dt,
								  const RealType* predicted_u,
								  const RealType* predicted_v,
								  const RealType* cell_variable,
								  const RealType* mass_corner_fluxes, 
								  const RealType* variable_gradient_diag, 
								  const RealType* variable_gradient_antidiag, 
								  RealType* variable_flux_corner);


#ifdef __cplusplus
}
#endif


#endif



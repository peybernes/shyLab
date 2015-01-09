#ifndef KERNEL_DIRECT_PROJECTION_CORNER_FLUXES_2D_BOUNDARY_INTERFACE_H
#define KERNEL_DIRECT_PROJECTION_CORNER_FLUXES_2D_BOUNDARY_INTERFACE_H
#include <string>
typedef int index_t;

#include "kernel_direct_projection_corner_fluxes_2d_boundary_conditions.h"
#include "kernel_direct_projection_2d_boundary_conditions.h"
#include "kernel_lagrange_2d.h"



void ReconstructGradientDiagBoundary(std::string BoundaryConditions,
				     index_t nx,
				     index_t ny,
				     RealType dx,
				     RealType dy,
				     RealType dt,
				     const RealType* predicted_u,
				     const RealType* predicted_v,
				     const RealType* lag_variable,
				     RealType* gradient_variable_diag,
				     RealType* gradient_variable_antidiag);


void ReconstructMassFluxOrder2XCornerFluxesBoundary(std::string BoundaryConditions,
						    index_t nx, 
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


void ReconstructMassFluxOrder2YCornerFluxesBoundary(std::string BoundaryConditions,
						    index_t nx, 
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


void ReconstructMassFluxOrder2CornerBoundary(std::string BoundaryConditions,
					     index_t nx, 
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


void ReconstructMassFluxOrder2DiagCornerBoundary(std::string BoundaryConditions,
						 index_t nx, 
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


void ProjectMassCornerFluxesBoundary(std::string BoundaryConditions,
				     index_t nx, 
				     index_t ny,
				     RealType dt,
				     const RealType* mass_corner_fluxes,
				     const RealType* in_vx, 
				     const RealType* in_vy, 
				     RealType* out_cell_mass);


void ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary(std::string BoundaryConditions,
								 index_t nx, 
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


void ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary(std::string BoundaryConditions,
								 index_t nx, 
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


void ReconstructIntensiveVariableFluxOrder2CornerBoundary(std::string BoundaryConditions,
							  index_t nx, 
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


void ReconstructIntensiveVariableFluxOrder2DiagCornerBoundary(std::string BoundaryConditions,
							      index_t nx, 
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


void MassProjectIntensiveVariableCornerFluxesBoundary(std::string BoundaryConditions,
						      index_t nx, 
						      index_t ny,
						      RealType dt,
						      const RealType* in_cell_variable,
						      const RealType* variable_flux_corner, 
						      const RealType* in_vx, 
						      const RealType* in_vy, 
						      const RealType* out_cell_mass,
						      RealType* out_cell_variable);


void ProjectNodalIntensiveVariableCornerFluxesBoundary(std::string BoundaryConditions,
						       index_t nx, 
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
#endif 

#ifndef KERNEL_AD_PROJECTION_2D_BOUNDARY_INTERFACE_H
#define KERNEL_AD_PROJECTION_2D_BOUNDARY_INTERFACE_H
#include <string>
typedef int index_t;

#include "kernel_ad_projection_2d_boundary_conditions.h"
#include "kernel_lagrange_2d.h"

void BoundaryVelocityPrediction(std::string BoundaryConditions,
					int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					const RealType* in_mass,
					const RealType* in_pressure,
					const RealType* in_pseudo_pressure,
					const RealType* in_velocity_x,
					const RealType* in_velocity_y,
					RealType* out_velocity_x,
					RealType* out_velocity_y);

void BoundaryCopy(std::string BoundaryConditions,
		  int nx,
		  int ny,
		  RealType* in_velocity_x,
		  RealType* in_velocity_y);

void ReconstructGradientXBoundary(std::string BoundaryConditions, index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);


void ReconstructMassFluxOrder2XBoundary(std::string BoundaryConditions,
					index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2XBoundary(std::string BoundaryConditions,
						     index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);

void ReconstructGradientNodalXBoundary(std::string BoundaryConditions,
				       index_t nx,
				       index_t ny,
				       RealType dx,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2XBoundary(std::string BoundaryConditions,
						  std::string variable_projected,
						  index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable,
						  const RealType* mass_flux,
						  RealType* out_variable);

void ReconstructGradientYBoundary(std::string BoundaryConditions,
				  index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);

void ReconstructMassFluxOrder2YBoundary(std::string BoundaryConditions,
					index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2YBoundary(std::string BoundaryConditions,
						     index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);


void ReconstructGradientNodalYBoundary(std::string BoundaryConditions,
				       index_t nx,
				       index_t ny,
				       RealType dy,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2YBoundary(std::string BoundaryConditions,
						  std::string variable_projected,
						  index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable,
						  const RealType* mass_flux,
						  RealType* out_variable);

#endif 

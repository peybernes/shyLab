#ifndef KERNEL_AD_PROJECTION_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_AD_PROJECTION_2D_BOUNDARY_CONDITIONS_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif
#include <string.h>

  // Periodic conditions

void ReconstructMassFluxOrder1PeriodicBoundaryX(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* volume_fluxes_x,
					const RealType* cell_variable,
					RealType* mass_flux);

void ReconstructMassFluxOrder1PeriodicBoundaryY(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* volume_fluxes_y,
					const RealType* cell_variable,
					RealType* mass_flux);

void ReconstructIntensiveVariableFluxOrder1PeriodicBoundaryX(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ReconstructIntensiveVariableFluxOrder1PeriodicBoundaryY(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ProjectNodalIntensiveVariablePeriodicBoundaryX(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,
					    const RealType* in_vx,
					    const RealType*  mass_flux,
					    RealType* out_vx);

void ProjectNodalIntensiveVariablePeriodicBoundaryY(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,
					    const RealType* in_vy,
					    const RealType*  mass_flux,
					    RealType* out_vy);

void ReconstructGradientXPeriodicBoundary(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);

void ReconstructMassFluxOrder2XPeriodicBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2XPeriodicBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);

void ReconstructGradientNodalXPeriodicBoundary(index_t nx,
				       index_t ny,
				       RealType dx,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2XPeriodicBoundary(index_t nx, 
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


void ReconstructGradientYPeriodicBoundary(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);


void ReconstructGradientYPeriodicBoundary(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);


void ReconstructMassFluxOrder2YPeriodicBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2YPeriodicBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);



void ReconstructGradientNodalYPeriodicBoundary(index_t nx,
				       index_t ny,
				       RealType dy,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2YPeriodicBoundary(index_t nx, 
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

  //Rayleigh Taylor conditions

void ReconstructGradientXBoundaryRt(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);


void ReconstructMassFluxOrder2XBoundaryRt(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2XBoundaryRt(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);

void ReconstructGradientNodalXBoundaryRt(index_t nx,
				       index_t ny,
				       RealType dx,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2XBoundaryRt(index_t nx, 
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

void ReconstructGradientYBoundaryRt(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* volume_fluxes,
				  const RealType* lag_variable,
				  RealType* gradient_variable);

void ReconstructMassFluxOrder2YBoundaryRt(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2YBoundaryRt(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);


void ReconstructGradientNodalYBoundaryRt(index_t nx,
				       index_t ny,
				       RealType dy,
				       RealType dt,
				       const RealType* predicted_velocity,
				       const RealType*  lag_nodal_variable,
				       RealType*  gradient_variable);


void ProjectNodalIntensiveVariableOrder2YBoundaryRt(index_t nx, 
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


  //Wall conditions

void ReconstructMassFluxOrder2XWallBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2XWallBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);

void ProjectNodalIntensiveVariableUxXWallBoundary(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
					    const RealType* mass_flux,
					    RealType* out_vx);

void ProjectNodalIntensiveVariableUyXWallBoundary(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
					    const RealType* mass_flux,
					    RealType* out_vx);

void ProjectNodalIntensiveVariableOrder2UxXWallBoundary(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dy,
							const RealType dt,
							const RealType* lag_cell_mass,
							const RealType* out_cell_mass,
							const RealType* u_velocity_pred,
							const RealType* in_variable,
							const RealType* gradient_variable,
							const RealType* mass_flux,
							RealType* out_variable);

void ProjectNodalIntensiveVariableOrder2UyXWallBoundary(index_t nx, 
							index_t ny, 
							index_t halo_width,
							const RealType dy,
							const RealType dt,
							const RealType* lag_cell_mass,
							const RealType* out_cell_mass,
							const RealType* u_velocity_pred,
							const RealType* in_variable,
							const RealType* gradient_variable,
							const RealType* mass_flux,
							RealType* out_variable);

void ReconstructMassFluxOrder2YWallBoundary(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType*  volume_fluxes,
					const RealType*  cell_density,
					const RealType*  cell_mass_gradient,
					RealType*  mass_flux);

void ReconstructIntensiveVariableFluxOrder2YWallBoundary(index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType*  volume_fluxes,
						     const RealType*  mass_flux,
						     const RealType*  cell_variable,
						     const RealType*  cell_variable_gradient,
						     RealType*  variable_flux);

void ProjectNodalIntensiveVariableUxYWallBoundary(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,
					    const RealType* in_vy,
					    const RealType* mass_flux,
					    RealType* out_vy);

void ProjectNodalIntensiveVariableUyYWallBoundary(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,
					    const RealType* in_vy,
					    const RealType* mass_flux,
					    RealType* out_vy);

void ProjectNodalIntensiveVariableOrder2UxYWallBoundary(index_t nx, 
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

void ProjectNodalIntensiveVariableOrder2UyYWallBoundary(index_t nx, 
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

#ifdef __cplusplus
}
#endif


#endif // BOUNDARY_RAYLEIGH_TAYLOR_CONDITIONS_H



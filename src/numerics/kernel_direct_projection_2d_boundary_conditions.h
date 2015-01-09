#ifndef KERNEL_DIRECT_PROJECTION_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_DIRECT_PROJECTION_2D_BOUNDARY_CONDITIONS_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif
#include <string.h>

  // Periodic conditions


void ProjectNodalIntensiveVariablePeriodicBoundaryDirect(index_t nx, 
				                 index_t ny, 
				                 index_t halo_width,
				                 const RealType* lag_cell_mass,
				                 const RealType* out_cell_mass,
				                 const RealType* in_vx,
				                 const RealType*  mass_flux_x,
				                 const RealType*  mass_flux_y,
				                 RealType* out_vx);



void ProjectNodalIntensiveVariableOrder2PeriodicBoundaryDirect(index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable_x,
						  const RealType* gradient_variable_y,
						  const RealType* mass_flux_x,
						  const RealType* mass_flux_y,
						  RealType* out_variable);


  //Rayleigh Taylor conditions


void ProjectNodalIntensiveVariableOrder2BoundaryRtDirect(index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dy,
						  const RealType dt,
						  const RealType* lag_cell_mass,
						  const RealType* out_cell_mass,
						  const RealType* u_velocity_pred,
						  const RealType* v_velocity_pred,
						  const RealType* in_variable,
						  const RealType* gradient_variable_x,
						  const RealType* gradient_variable_y,
						  const RealType* mass_flux_x,
						  const RealType* mass_flux_y,
						  RealType* out_variable);



  //Wall conditions


void ProjectNodalIntensiveVariableUxWallBoundaryDirect(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
						       //const RealType* in_vy,
					    const RealType* mass_flux_x,
					    const RealType* mass_flux_y,
						       // RealType* out_vx,
					    RealType* out_vy);

void ProjectNodalIntensiveVariableUyWallBoundaryDirect(index_t nx, 
					    index_t ny, 
					    index_t halo_width,
					    const RealType* lag_cell_mass,
					    const RealType* out_cell_mass,  
					    const RealType* in_vx,
						       //const RealType* in_vy,
					    const RealType* mass_flux_x,
					    const RealType* mass_flux_y,
						       //RealType* out_vx,
					    RealType* out_vy);


void ProjectNodalIntensiveVariableOrder2UxWallBoundaryDirect(index_t nx, 
							     index_t ny, 
							     index_t halo_width,
							     const RealType dx,
							     const RealType dy,
							     const RealType dt,
							     const RealType* lag_cell_mass,
							     const RealType* out_cell_mass,  
							     const RealType* u_velocity_pred,
							     const RealType* v_velocity_pred,
							     const RealType* in_variable,
							     const RealType* gradient_variable_x,
							     const RealType* gradient_variable_y,
							     const RealType* mass_flux_x,
							     const RealType* mass_flux_y,
							     RealType* out_variable);



void ProjectNodalIntensiveVariableOrder2UyWallBoundaryDirect(index_t nx, 
							     index_t ny, 
							     index_t halo_width,
							     const RealType dx,
							     const RealType dy,
							     const RealType dt,
							     const RealType* lag_cell_mass,
							     const RealType* out_cell_mass,  
							     const RealType* u_velocity_pred,
							     const RealType* v_velocity_pred,
							     const RealType* in_variable,
							     const RealType* gradient_variable_x,
							     const RealType* gradient_variable_y,
							     const RealType* mass_flux_x,
							     const RealType* mass_flux_y,
							     RealType* out_variable);


#ifdef __cplusplus
}
#endif


#endif // BOUNDARY_RAYLEIGH_TAYLOR_CONDITIONS_H



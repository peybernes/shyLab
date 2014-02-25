#ifndef KERNEL_AD_PROJECTION_2D_H
#define KERNEL_AD_PROJECTION_2D_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif

void ComputeDirectionalLagrangianQuantitiesX(index_t nx, 
					     index_t ny, 
					     RealType dt,
					     RealType dx,
					     RealType dy,
					     const RealType* in_vx,
					     const RealType* cell_mass,
					     RealType* volume_fluxes,
					     RealType* directional_lagrangian_volume,
					     RealType* directional_lagrangian_density);

void ComputeDirectionalLagrangianQuantitiesY(index_t nx, 
					     index_t ny, 
					     RealType dt,
					     RealType dx,
					     RealType dy,
					     const RealType* in_vy,
					     const RealType* cell_mass,
					     RealType* volume_fluxes,
					     RealType* directional_lagrangian_volume,
					     RealType* directional_lagrangian_density);


void ReconstructMassFluxOrder1X(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType* volume_fluxes_x,
				const RealType* cell_variable,
				RealType* mass_flux);

void ReconstructMassFluxOrder1Y(index_t nx, 
				index_t ny, 
				index_t halo_width,
				const RealType* volume_fluxes_y,
				const RealType* cell_variable,
				RealType* mass_flux);

void ReconstructMassFluxOrder1BoundaryX(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* volume_fluxes_x,
					const RealType* cell_variable,
					RealType* mass_flux);

void ReconstructMassFluxOrder1BoundaryY(index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType* volume_fluxes_y,
					const RealType* cell_variable,
					RealType* mass_flux);

void ReconstructIntensiveVariableFluxOrder1X(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ReconstructIntensiveVariableFluxOrder1Y(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ReconstructIntensiveVariableFluxOrder1BoundaryX(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ReconstructIntensiveVariableFluxOrder1BoundaryY(index_t nx, 
					     index_t ny, 
					     index_t halo_width,
					     const RealType* mass_flux,
					     const RealType* cell_variable,
					     RealType* variable_flux);

void ProjectMassX(index_t nx, 
			       index_t ny, 
			       const RealType* in_cell_mass,
			       const RealType* mass_flux,
			       RealType* out_cell_mass);

void ProjectMassY(index_t nx, 
			       index_t ny, 
			       const RealType* in_cell_mass,
			       const RealType* mass_flux,
			       RealType* out_cell_mass);

void CheckFluxPeriodicalPropertyX(index_t nx,
				  index_t ny,
				  const RealType* flux);

void CheckFluxPeriodicalPropertyY(index_t nx,
				  index_t ny,
				  const RealType* flux);

void ProjectNodalIntensiveVariableX(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* in_cell_mass,
				    const RealType* in_vx,
				    const RealType*  mass_flux,
				    RealType* out_vx);

void ProjectNodalIntensiveVariableY(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* in_cell_mass,
				    const RealType* in_vx,
				    const RealType*  mass_flux,
				    RealType* out_vx);

void ProjectNodalIntensiveVariableBoundaryX(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* in_cell_mass,
				    const RealType* in_vx,
				    const RealType*  mass_flux,
				    RealType* out_vx);

void ProjectNodalIntensiveVariableBoundaryY(index_t nx, 
				    index_t ny, 
				    index_t halo_width,
				    const RealType* in_cell_mass,
				    const RealType* in_vx,
				    const RealType*  mass_flux,
				    RealType* out_vx);

void MassProjectIntensiveVariableX(index_t nx, 
				   index_t ny, 
				   const RealType* in_cell_mass,
				   //const RealType* mass_flux,
				   const RealType* in_cell_variable,
				   const RealType* in_face_variable,
				   const RealType* out_cell_mass,
				   RealType* out_cell_variable);

void MassProjectIntensiveVariableY(index_t nx, 
				   index_t ny, 
				   const RealType* in_cell_mass,
				   //const RealType* mass_flux,
				   const RealType* in_cell_variable,
				   const RealType* in_face_variable,
				   const RealType* out_cell_mass,
				   RealType* out_cell_variable);
 

#ifdef __cplusplus
}
#endif


#endif // KERNEL_AD_PROJECTION_2D_H

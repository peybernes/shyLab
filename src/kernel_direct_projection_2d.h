#ifndef KERNEL_DIRECT_PROJECTION_2D_H
#define KERNEL_DIRECT_PROJECTION_2D_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif


void ProjectMassDirect(index_t nx, 
			       index_t ny, 
			       const RealType* in_cell_mass,
			       const RealType* mass_flux_x, 
			       const RealType* mass_flux_y, 
			       RealType* out_cell_mass);


void ProjectNodalIntensiveVariableDirect(index_t nx, 
				         index_t ny, 
				         index_t halo_width,
				         const RealType* lag_cell_mass,
				         const RealType* out_cell_mass,
				         const RealType* in_vx,
				         const RealType* in_vy,
				         const RealType*  mass_flux_x,
				         const RealType*  mass_flux_y,
				         RealType* out_vx,
				         RealType* out_vy);

void ProjectNodalIntensiveVariablePeriodicBoundaryDirect(index_t nx, 
				                 index_t ny, 
				                 index_t halo_width,
				                 const RealType* lag_cell_mass,
				                 const RealType* out_cell_mass,
				                 const RealType* in_vx,
				                 const RealType* in_vy,
				                 const RealType*  mass_flux_x,
				                 const RealType*  mass_flux_y,
				                 RealType* out_vx,
				                 RealType* out_vy);

void MassProjectIntensiveVariableDirect(index_t nx, 
				        index_t ny, 
				        const RealType* in_cell_mass,
				        //const RealType* mass_flux_x,
				        //const RealType* mass_flux_y,
				        const RealType* in_cell_variable,
				        const RealType* in_face_variable_x,
				        const RealType* in_face_variable_y,
				        const RealType* out_cell_mass,
				        RealType* out_cell_variable);

 

#ifdef __cplusplus
}
#endif


#endif // KERNEL_DIRECT_PROJECTION_2D_H

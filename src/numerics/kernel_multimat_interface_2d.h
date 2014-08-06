#ifndef KERNEL_MULTIMAT_INTERFACE_2D_H
#define KERNEL_MULTIMAT_INTERFACE_2D_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif



  // ========================
  //   General functions for Geometry    
  // ========================


int sign(RealType x);


void PlaceInterface(// in
		    const RealType xmin,
		    const RealType xmax,
		    const RealType ymin,
		    const RealType ymax,
		    const RealType nX,
		    const RealType nY,
		    const RealType vol_fraction,
		    // out
		    RealType* pointA,
		    RealType* pointB);


void IntersectLineRectangle(// in
			    const RealType* point1,
			    const RealType* point2,
			    const RealType xmin,
			    const RealType xmax,
			    const RealType ymin,
			    const RealType ymax,
			    // out
			    RealType* pointA,
			    RealType* pointB);


RealType ComputePartialVolume(// in
			      const RealType* pointA,
			      const RealType* pointB,
			      const RealType xmin,
			      const RealType xmax,
			      const RealType ymin,
			      const RealType ymax,			  
			      const RealType nX,			  
			      const RealType nY);




  // ========================
  //    Functions for Projection phase   
  // ========================


void ReconstructGradientMultimatX(index_t nx,
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


void ReconstructGradientMultimatY(index_t nx,
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


void ReconstructNormalInterface(// in
				index_t nx,
				index_t ny,
				RealType dx,
				RealType dy,
				const RealType* in_variable,
				// out
				RealType* interface_normal_x,
				RealType* interface_normal_y);


void ComputeDirectionalLagrangianQuantitiesMultimatX(index_t nx, 
						     index_t ny, 
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     const RealType* in_vx,
						     const RealType* cell_mass_1,
						     const RealType* cell_mass_2,
						     const RealType* cell_volumic_fraction,
						     RealType* volume_fluxes,
						     RealType* directional_lagrangian_volume,
						     RealType* directional_lagrangian_density_1,
						     RealType* directional_lagrangian_density_2);


void ComputeDirectionalLagrangianQuantitiesMultimatY(index_t nx, 
						     index_t ny, 
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     const RealType* in_vy,
						     const RealType* cell_mass_1,
						     const RealType* cell_mass_2,
						     const RealType* cell_volumic_fraction,
						     RealType* volume_fluxes,
						     RealType* directional_lagrangian_volume,
						     RealType* directional_lagrangian_density_1,
						     RealType* directional_lagrangian_density_2);


void ComputeDirectionalLagrangianFractionalVolumeFluxesX(// in
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



void ComputeDirectionalLagrangianFractionalVolumeFluxesY(// in
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

void CheckDirectionalVolumeFluxesX(//in
				   index_t nx, 
				   index_t ny, 
				   RealType dx,
				   RealType dy,
				   const RealType* directional_lagrangian_volume,
				   const RealType* cell_volumic_fraction,
				   // out
				   RealType* volume_fluxes_1,
				   RealType* volume_fluxes_2);


void CheckDirectionalVolumeFluxesY(//in
				   index_t nx, 
				   index_t ny, 
				   RealType dx,
				   RealType dy,
				   const RealType* directional_lagrangian_volume,
				   const RealType* cell_volumic_fraction,
				   // out
				   RealType* volume_fluxes_1_y,
				   RealType* volume_fluxes_2_y);


void CheckDirectionalMassFluxX(//in
			       index_t nx, 
			       index_t ny, 
			       RealType dx,
			       RealType dy,
			       const RealType* cell_mass_1,
			       const RealType* cell_mass_2,
			       const RealType* directional_lagrangian_volume,
			       const RealType* cell_volumic_fraction,
			       const RealType* volume_fluxes_1,
			       const RealType* volume_fluxes_2,
			       // out
			       RealType* mass_flux_1,
			       RealType* mass_flux_2);


void CheckDirectionalMassFluxY(//in
			       index_t nx, 
			       index_t ny, 
			       RealType dx,
			       RealType dy,
			       const RealType* cell_mass_1,
			       const RealType* cell_mass_2,
			       const RealType* directional_lagrangian_volume,
			       const RealType* cell_volumic_fraction,
			       const RealType* volume_fluxes_1_y,
			       const RealType* volume_fluxes_2_y,
			       // out
			       RealType* mass_flux_1_y,
			       RealType* mass_flux_2_y);


void ProjectDensityMultimatX(// in
			     index_t nx, 
			     index_t ny, 
			     RealType dx,
			     RealType dy,
			     const RealType* directional_lagrangian_volume,
			     const RealType* in_cell_volumic_fraction,
			     const RealType* volume_fluxes_1_x,			     
			     const RealType* volume_fluxes_2_x,
			     const RealType* in_rho_1,
			     const RealType* in_rho_2,
			     const RealType* in_cell_mass,
			     const RealType* mass_flux_x,
			     // out
			     RealType* out_cell_mass_1,
			     RealType* out_cell_mass_2,
			     RealType* out_cell_volumic_fraction,
			     RealType* out_cell_mass,
			     RealType* out_rho_1,
			     RealType* out_rho_2);


void ProjectDensityMultimatY(// in
			     index_t nx, 
			     index_t ny, 
			     RealType dx,
			     RealType dy,
			     const RealType* directional_lagrangian_volume,
			     const RealType* in_cell_volumic_fraction,
			     const RealType* volume_fluxes_1_y,			     
			     const RealType* volume_fluxes_2_y,
			     const RealType* in_rho_1,
			     const RealType* in_rho_2,
			     const RealType* in_cell_mass,
			     const RealType* mass_flux_y,
			     // out
			     RealType* out_cell_mass_1,
			     RealType* out_cell_mass_2,
			     RealType* out_cell_volumic_fraction,
			     RealType* out_cell_mass,
			     RealType* out_rho_1,
			     RealType* out_rho_2);



  // ========================
  //    Functions for Lagrangian phase  
  // ========================


void LagrangePressurePredictedOptimisedMultimat(// in
						int nx,
						int ny,
						RealType dt,
						RealType dx,
						RealType dy,
						RealType gamma_1,
						RealType gamma_2,
						RealType pi_1,
						RealType pi_2,
						const RealType* in_mass_1,
						const RealType* in_mass_2,
						const RealType* in_rho_1,
						const RealType* in_rho_2,
						const RealType* in_energy_1,	 
						const RealType* in_energy_2,	 
						const RealType* cell_volumic_fraction,	 
						const RealType* in_velocity_x,
						const RealType* in_velocity_y,
						//out
						RealType* out_pressure,
						RealType* out_pressure_1,
						RealType* out_pressure_2,
						RealType* out_predicted_pressure,
						RealType* out_predicted_pressure_1,
						RealType* out_predicted_pressure_2,
						RealType* out_pseudo_pressure);



void LagrangeCorrectionOptimisedMultimat(// in
					 int nx,
					 int ny,
					 RealType dt,
					 RealType dx,
					 RealType dy,
					 const RealType* in_mass_1,
					 const RealType* in_mass_2,
					 const RealType* in_energy_1,
					 const RealType* in_energy_2,
					 const RealType* in_pressure_1,
					 const RealType* in_pressure_2,
					 const RealType* in_pseudo_pressure,
					 const RealType* cell_volumic_fraction,
					 const RealType* in_velocity_x,
					 const RealType* in_velocity_y,
					 // out
					 RealType* out_energy,
					 RealType* out_energy_1,
					 RealType* out_energy_2);


#ifdef __cplusplus
}
#endif


#endif // KERNEL_PROJECTION_MULTIMAT_INTERFACE_2D_H

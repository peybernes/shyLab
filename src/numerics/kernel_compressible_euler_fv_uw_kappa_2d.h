#ifndef KERNEL_COMPRESSIBLE_EULER_FV_UW_KAPPA_2D_H
#define KERNEL_COMPRESSIBLE_EULER_FV_UW_KAPPA_2D_H

#ifdef __cplusplus
extern "C" {
#endif

void CompressibleEulerPhysicalToConservative(int nx, 
					     int ny, 
					     const RealType* in_rho, 
					     const RealType* in_velocity_x,
					     const RealType* in_velocity_y,
					     const RealType* in_total_energy,
					     RealType* out_rho,
					     RealType* out_velocity_x,
					     RealType* out_velocity_y,
					     RealType* out_total_energy);

void CompressibleEulerConservativeToPhysical(int nx, 
					     int ny, 
					     const RealType* out_rho,
					     RealType* out_velocity_x,
					     RealType* out_velocity_y,
					     RealType* out_total_energy);

void CompressibleEulerFvUwKappa2dX(int nx, 
				   int ny, 
				   int halo_width,
				   RealType dt,
				   RealType dx,
				   RealType dy,
				   RealType kappa,
				   const RealType* in_rho, 
				   const RealType* in_velocity_x,
				   const RealType* in_velocity_y,
				   const RealType* in_total_energy,
				   RealType* out_rho,
				   RealType* out_velocity_x,
				   RealType* out_velocity_y,
				   RealType* out_total_energy);

void CompressibleEulerFvUwKappa2dY(int nx, 
				   int ny, 
				   int halo_width,
				   RealType dt,
				   RealType dx,
				   RealType dy,
				   RealType kappa,
				   const RealType* in_rho, 
				   const RealType* in_velocity_x,
				   const RealType* in_velocity_y,
				   const RealType* in_total_energy,
				   RealType* out_rho,
				   RealType* out_velocity_x,
				   RealType* out_velocity_y,
				   RealType* out_total_energy);

void CompressibleEulerFvUwKappa2dBoundaryConditionsX(int nx, 
						     int ny, 
						     int halo_width,
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     RealType kappa,
						     const RealType* in_rho, 
						     const RealType* in_velocity_x,
						     const RealType* in_velocity_y,
						     const RealType* in_total_energy,
						     RealType* out_rho,
						     RealType* out_velocity_x,
						     RealType* out_velocity_y,
						     RealType* out_total_energy);

void CompressibleEulerFvUwKappa2dBoundaryConditionsY(int nx, 
						     int ny, 
						     int halo_width,
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     RealType kappa,
						     const RealType* in_rho, 
						     const RealType* in_velocity_x,
						     const RealType* in_velocity_y,
						     const RealType* in_total_energy,
						     RealType* out_rho,
						     RealType* out_velocity_x,
						     RealType* out_velocity_y,
						     RealType* out_total_energy);

#ifdef __cplusplus
}
#endif


#endif // KERNEL_COMPRESSIBLE_EULER_FV_UW_KAPPA_2D_H

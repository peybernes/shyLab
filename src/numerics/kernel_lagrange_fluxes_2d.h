#ifndef KERNEL_LAGRANGE_FLUXES_2D_H
#define KERNEL_LAGRANGE_FLUXES_2D_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif

RealType TimeStepLFMix(int nx,
		       int ny,
		       const RealType dx,
		       const RealType dy,
		       const RealType CFL,
		       RealType* gamma_mix,
		       const RealType gamma_1,
		       const RealType gamma_2,
		       const RealType pi_1,
		       const RealType pi_2,
		       const RealType* density_1,
		       const RealType* density_2,
		       const RealType* pressure_1,
		       const RealType* pressure_2,
		       const RealType* in_c_1,
		       const RealType* in_c_2,
		       const RealType* in_velocity_x,
		       const RealType* in_velocity_y);

void ComputeGradientAlpha(index_t nx, 
			  index_t ny, 
			  index_t nb_mat, 
			  RealType dx,
			  RealType dy,
			  RealType** in_c_k,
			  RealType** alphak_gradx_left,
			  RealType** alphak_gradx_right,
			  RealType** alphak_grady_bot,
			  RealType** alphak_grady_top
			  );
  
#ifdef __cplusplus
}
#endif


#endif // KERNEL_LAGRANGE_FLUXES_2D_H

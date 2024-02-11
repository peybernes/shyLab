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
		       const int nb_mat,
		       const RealType CFL,
		       RealType* gamma_mix,
		       RealType* pi_prime_mix,
		       RealType* speed_of_sound_mix,
		       const RealType* gamma_k,
		       const RealType* pi_prime_k,
		       RealType* density,
		       RealType* pressure,
		       RealType** in_c_k,
		       const RealType* in_velocity_x,
		       const RealType* in_velocity_y);

void ComputeGradientAlpha(index_t nx, 
			  index_t ny, 
			  index_t nb_mat, 
			  const RealType dx,
			  const RealType dy,
			  RealType** in_c_k,
			  RealType** alphak_gradx_left,
			  RealType** alphak_gradx_right,
			  RealType** alphak_grady_bot,
			  RealType** alphak_grady_top
			  );
  
void ComputeGradientPplusPiPrime(index_t nx, 
			  index_t ny, 
			  const RealType dx,
			  const RealType dy,
			  RealType*  pressure,
			  RealType*  pi_prime_mix, 	 
			  RealType*  p_plus_pi_prime,	 
			  RealType*  p_plus_pi_prime_gradx,
			  RealType*  p_plus_pi_prime_grady	 
			  );

void ComputeGradientP(index_t nx, 
			  index_t    ny, 
			  RealType   dx,
			  RealType   dy,
			  int        nb_mat,	 
		          RealType*  pi_prime_mix,
			  RealType*  pressure,
		          RealType**  in_c_k,
		          RealType*  gamma_k,
			  RealType*  pi_prime_k,
			  RealType**  alphak_gradx_left,
			  RealType**  alphak_gradx_right,
			  RealType**  alphak_grady_bot,
			  RealType**  alphak_grady_top,
			  RealType*  p_plus_pi_prime_gradx,
			  RealType*  p_plus_pi_prime_grady,
			  RealType*  rho_e,
			  RealType*  rho_e_gradx_left,
			  RealType*  rho_e_gradx_right,
			  RealType*  rho_e_grady_bot,
			  RealType*  rho_e_grady_top,
			  RealType*  p_gradx_left,
			  RealType*  p_gradx_right,
			  RealType*  p_grady_bot,
		      RealType*  p_grady_top);

#ifdef __cplusplus
}
#endif


#endif // KERNEL_LAGRANGE_FLUXES_2D_H

#ifndef KERNEL_LAGRANGE_FLUXES_2D_BOUNDARY_CONDITIONS_H
#define KERNEL_LAGRANGE_FLUXES_2D_BOUNDARY_CONDITIONS_H

typedef int index_t;

#include"kernels.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <string.h>

void ComputeGradientAlphaBoundaryWall(index_t nx, 
			  index_t ny, 
			  index_t nb_mat, 
			  RealType dx,
			  RealType dy,
			  RealType** in_c_k,
			  RealType** alphak_gradx_left,
			  RealType** alphak_gradx_right,
			  RealType** alphak_grady_bot,
			  RealType** alphak_grady_top);

void ComputeGradientPplusPiPrimeBoundaryWall(index_t nx, 
			  index_t ny, 
			  RealType dx,
			  RealType dy,
			  RealType*  pressure,
			  RealType*  pi_prime_mix, 	 
			  RealType*  p_plus_pi_prime,	 
			  RealType*  p_plus_pi_prime_gradx,
			  RealType*  p_plus_pi_prime_grady	 
					     );

void ComputeGradientPBoundaryWall(index_t nx, 
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
  
void ComputeGradientRhoBoundaryWall(index_t nx, 
				  index_t    ny, 
				  RealType   dx,
				  RealType   dy,
				  int        nb_mat,
				  RealType** in_rho_k,
				  RealType**  in_c_k,
				  RealType**  alphak_gradx_left,
				  RealType**  alphak_gradx_right,
				  RealType**  alphak_grady_bot,
				  RealType**  alphak_grady_top,
				  RealType*  rho_gradx_left,
				  RealType*  rho_gradx_right,
				  RealType*  rho_grady_top,
				  RealType*  rho_grady_bot,
				  RealType**  rhok_gradx,
				  RealType**  rhok_grady
				    );

void ComputeHLLBoundaryPeriodic(index_t nx, 
		index_t    ny, 
		RealType   dx,
		RealType   dy,
		RealType*  rho_gradx_left,
		RealType*  rho_gradx_right,
		RealType*  rho_grady_top,
		RealType*  rho_grady_bot,
		RealType*  in_rho,
		RealType*  pressure,
		RealType*  in_u_cell,
		RealType*  in_v_cell,
		RealType*  speed_of_sound_mix,
		RealType*  p_gradx_left,
		RealType*  p_gradx_right,
		RealType*  p_grady_bot,
		RealType*  p_grady_top,
		RealType*  u_gradx_left,
		RealType*  u_gradx_right,
		RealType*  v_grady_bot,
		RealType*  v_grady_top,
		RealType*  p_xet,
		RealType*  p_yet,
		RealType*  u_et,
		RealType*  v_et
			    );

#ifdef __cplusplus
}
#endif


#endif



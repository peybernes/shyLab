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
			  const RealType* in_c_1,
			  const RealType* in_c_2,
			  RealType** alphak_gradx_left,
			  RealType** alphak_gradx_right,
			  RealType** alphak_grady_bot,
			  RealType** alphak_grady_top);

#ifdef __cplusplus
}
#endif


#endif



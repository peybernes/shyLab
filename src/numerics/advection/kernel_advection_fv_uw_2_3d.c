#include "kernels.hpp"

#define RESTRICT __restrict__

#include <assert.h>
#include <math.h>


void AdvectionFvUw23d(int nx, 
		      int ny, 
		      int nz, 
		      RealType dt,
		      RealType dx,
		      RealType dy,
		      RealType dz,
		      const RealType* in_value, 
		      const RealType* in_velocity_x,
		      const RealType* in_velocity_y,
		      const RealType* in_velocity_z,
		      RealType* out_value) {

  assert(0);

}

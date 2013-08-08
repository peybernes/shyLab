#ifndef KERNEL_ADVECTION_FV_UW_2_2D_H
#define KERNEL_ADVECTION_FV_UW_2_2D_H

void AdvectionFvUw22d(int nx, 
		      int ny, 
		      RealType dt,
		      RealType dx,
		      RealType dy,
		      const RealType* in_value, 
		      const RealType* in_velocity_x,
		      const RealType* in_velocity_y,
		      RealType* out_value);

#endif // KERNEL_ADVECTION_FV_UW_2_2D_H

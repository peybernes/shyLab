#ifndef KERNEL_ADVECTION_FV_UW_1_3D_H
#define KERNEL_ADVECTION_FV_UW_1_3D_H

void AdvectionFvUw13d(int nx, 
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
		      RealType* out_value);

#endif // KERNEL_ADVECTION_FV_UW_1_3D_H

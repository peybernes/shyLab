#include "kernel_advection_fv_uw_kappa_3d.h"

#define RESTRICT __restrict__

#include <assert.h>
#include <math.h>

void AdvectionFvUwKappa3d(int nx, 
			  int ny, 
			  int nz, 
			  RealType dt,
			  RealType dx,
			  RealType dy,
			  RealType dz,
			  const RealType* RESTRICT in_value, 
			  const RealType* RESTRICT in_velocity_x,
			  const RealType* RESTRICT in_velocity_y,
			  const RealType* RESTRICT in_velocity_z,
			  RealType* RESTRICT out_value) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);
  
  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int i = 1; i < nx - 1; ++i) {

	const int ooo = (nx * ny * k) + (nx * j) + i;
	const int moo = (nx * ny * k) + (nx * j) + i - 1;
	const int poo = (nx * ny * k) + (nx * j) + i + 1;
	const int omo = (nx * ny * k) + (nx * (j - 1)) + i;
	const int opo = (nx * ny * k) + (nx * (j + 1)) + i;
	const int oom = (nx * ny * (k - 1)) + (nx * j) + i;
	const int oop = (nx * ny * (k + 1)) + (nx * j) + i;

	const int f_moo = (nx * ny * k) + (nx * j) + i;
	const int f_poo = (nx * ny * k) + (nx * j) + i + 1;
	const int f_omo = (nx * ny * k) + (nx * j) + i;
	const int f_opo = (nx * ny * k) + (nx * (j + 1)) + i;
	const int f_oom = (nx * ny * k) + (nx * j) + i;
	const int f_oop = (nx * ny * (k + 1)) + (nx * j) + i;

	const RealType u_face_moo = in_velocity_x[f_moo];
	const RealType u_face_poo = in_velocity_x[f_poo];
	const RealType u_face_omo = in_velocity_y[f_omo];
	const RealType u_face_opo = in_velocity_y[f_opo];
	const RealType u_face_oom = in_velocity_z[f_oom];
	const RealType u_face_oop = in_velocity_z[f_oop];

	const RealType in_cell_ooo = in_value[ooo];

	const RealType in_cell_moo = in_value[moo];
	const RealType in_cell_poo = in_value[poo];
	const RealType in_cell_omo = in_value[omo];
	const RealType in_cell_opo = in_value[opo];
	const RealType in_cell_oom = in_value[oom];
	const RealType in_cell_oop = in_value[oop];

	const RealType half = 0.5;

	const RealType flux_face_moo =
	  half * ((u_face_moo + fabs(u_face_moo)) * in_cell_moo + 
		 (u_face_moo - fabs(u_face_moo)) * in_cell_ooo);

	const RealType flux_face_poo =
	  half * ((u_face_poo + fabs(u_face_poo)) * in_cell_ooo + 
		 (u_face_poo - fabs(u_face_poo)) * in_cell_poo);

	const RealType flux_face_omo =
	  half * ((u_face_omo + fabs(u_face_omo)) * in_cell_omo + 
		 (u_face_omo - fabs(u_face_omo)) * in_cell_ooo);

	const RealType flux_face_opo =
	  half * ((u_face_opo + fabs(u_face_opo)) * in_cell_ooo + 
		 (u_face_opo - fabs(u_face_opo)) * in_cell_opo);

	const RealType flux_face_oom =
	  half * ((u_face_oom + fabs(u_face_oom)) * in_cell_oom + 
		 (u_face_oom - fabs(u_face_oom)) * in_cell_ooo);

	const RealType flux_face_oop =
	  half * ((u_face_oop + fabs(u_face_oop)) * in_cell_ooo + 
		 (u_face_oop - fabs(u_face_oop)) * in_cell_oop);

	const RealType c = 0.01;

	const RealType out_cell_ooo = in_cell_ooo + 
	  c * (flux_face_moo - flux_face_poo + 
	       flux_face_omo - flux_face_opo +
	       flux_face_oom - flux_face_oop);

	out_value[ooo] = out_cell_ooo;
	
      }
    }
  }

}

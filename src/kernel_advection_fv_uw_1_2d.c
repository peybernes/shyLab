#include "kernel_advection_fv_uw_1_2d.h"

#define RESTRICT __restrict__

#include <assert.h>
#include <math.h>

void AdvectionFvUw12d(int nx, 
		      int ny, 
		      RealType dt,
		      RealType dx,
		      RealType dy,
		      const RealType* RESTRICT in_value, 
		      const RealType* RESTRICT in_velocity_x,
		      const RealType* RESTRICT in_velocity_y,
		      RealType* RESTRICT out_value) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  for (int j = 1; j < ny - 1; ++j) {
    for (int i = 1; i < nx - 1; ++i) {

      const int oo = (nx * j) + i;

      const int mo = (nx * j) + i - 1;
      const int po = (nx * j) + i + 1;
      const int om = (nx * (j - 1)) + i;
      const int op = (nx * (j + 1)) + i;

      const int f_mo = (nx * j) + i;
      const int f_po = (nx * j) + i + 1;
      const int f_om = (nx * j) + i;
      const int f_op = (nx * j) + i + 1;

      const RealType u_face_mo = in_velocity_x[f_mo];
      const RealType u_face_po = in_velocity_x[f_po];
      const RealType u_face_om = in_velocity_y[f_om];
      const RealType u_face_op = in_velocity_y[f_op];

      const RealType in_cell_oo = in_value[oo];

      const RealType in_cell_mo = in_value[mo];
      const RealType in_cell_po = in_value[po];
      const RealType in_cell_om = in_value[om];
      const RealType in_cell_op = in_value[op];

      const RealType half = 0.5;

      const RealType flux_face_mo =
	half * cy * ((u_face_mo + fabs(u_face_mo)) * in_cell_mo + 
		     (u_face_mo - fabs(u_face_mo)) * in_cell_oo);

      const RealType flux_face_po =
	half * cy * ((u_face_po + fabs(u_face_po)) * in_cell_oo + 
		     (u_face_po - fabs(u_face_po)) * in_cell_po);

      const RealType flux_face_om =
	half * cx * ((u_face_om + fabs(u_face_om)) * in_cell_om + 
		     (u_face_om - fabs(u_face_om)) * in_cell_oo);

      const RealType flux_face_op =
	half * cx * ((u_face_op + fabs(u_face_op)) * in_cell_oo + 
		     (u_face_op - fabs(u_face_op)) * in_cell_op);
      
      const RealType out_cell_oo = in_cell_oo +
	dt * (flux_face_mo - flux_face_po + flux_face_om - flux_face_op);

      out_value[oo] = out_cell_oo;
 
    }
  }

}

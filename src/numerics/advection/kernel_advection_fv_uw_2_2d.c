#include "kernel_advection_fv_uw_2_2d.h"

#include <assert.h>
#include <math.h>

#define RESTRICT __restrict

void AdvectionFvUw22d(int nx, 
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

  for (int j = 2; j < ny - 2; ++j) {
    for (int i = 2; i < nx - 2; ++i) {

      const int oo = (nx * j) + i;

      const int m1o = (nx * j) + i - 1;
      const int p1o = (nx * j) + i + 1;
      const int m2o = (nx * j) + i - 2;
      const int p2o = (nx * j) + i + 2;

      const int om1 = (nx * (j - 1)) + i;
      const int op1 = (nx * (j + 1)) + i;
      const int om2 = (nx * (j - 2)) + i;
      const int op2 = (nx * (j + 2)) + i;

      const int f_mo = (nx * j) + i;
      const int f_po = (nx * j) + i + 1;
      const int f_om = (nx * j) + i;
      const int f_op = (nx * j) + i + 1;

      const RealType u_face_mo = in_velocity_x[f_mo];
      const RealType u_face_po = in_velocity_x[f_po];
      const RealType u_face_om = in_velocity_y[f_om];
      const RealType u_face_op = in_velocity_y[f_op];

      const RealType in_cell_oo = in_value[oo];

      const RealType in_cell_mo = in_value[m1o];
      const RealType in_cell_po = in_value[p1o];
      const RealType in_cell_om = in_value[om1];
      const RealType in_cell_op = in_value[op1];

      const RealType half = 0.5;
      
      const RealType in_cell_dx_oo = half * cx * (in_value[p1o] - in_value[m1o]);
      const RealType in_cell_dx_mo = half * cx * (in_value[oo] - in_value[m2o]);
      const RealType in_cell_dx_po = half * cx * (in_value[p2o] - in_value[oo]);

      const RealType in_cell_dy_oo = half * cy * (in_value[op1] - in_value[om1]);
      const RealType in_cell_dy_om = half * cy * (in_value[oo] - in_value[om2]);
      const RealType in_cell_dy_op = half * cy * (in_value[op2] - in_value[oo]);

      const RealType flux_face_mo =
	half * cy * ((u_face_mo + fabs(u_face_mo)) * (in_cell_mo + half * dx * in_cell_dx_mo) +
		     (u_face_mo - fabs(u_face_mo)) * (in_cell_oo - half * dx * in_cell_dx_oo));

      const RealType flux_face_po =
	half * cy * ((u_face_po + fabs(u_face_po)) * (in_cell_oo + half * dx * in_cell_dx_oo) + 
		     (u_face_po - fabs(u_face_po)) * (in_cell_po - half * dx * in_cell_dx_po));

      const RealType flux_face_om =
	half * cx * ((u_face_om + fabs(u_face_om)) * (in_cell_om + half * dy * in_cell_dy_om) + 
		     (u_face_om - fabs(u_face_om)) * (in_cell_oo - half * dy * in_cell_dy_oo));

      const RealType flux_face_op =
	half * cx * ((u_face_op + fabs(u_face_op)) * (in_cell_oo + half * dy * in_cell_dy_oo) + 
		     (u_face_op - fabs(u_face_op)) * (in_cell_op - half * dy * in_cell_dy_op));
      
      const RealType out_cell_oo = in_cell_oo +
	dt * (flux_face_mo - flux_face_po + flux_face_om - flux_face_op);

      out_value[oo] = out_cell_oo;
 
    }
  }

}

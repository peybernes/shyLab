// -*- c++ -*- (for emacs users)
#include "kernel_advection_fv_uw_kappa_2d.h"

#define NDEBUG
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>

#include "kernel_tools.h"

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a < _b ? _a : _b; })



/// Implements linear (but inhomogeneous) advection, using Van Leer
/// kappa schemes.
void AdvectionFvUwKappa2d(int nx, 
			  int ny, 
			  int halo_width,
			  RealType dt,
			  RealType dx,
			  RealType dy,
			  RealType kappa,
			  const RealType* RESTRICT in_value, 
			  const RealType* RESTRICT in_velocity_x,
			  const RealType* RESTRICT in_velocity_y,
			  RealType* RESTRICT out_value) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  const RealType half = 0.5;

  for (int iy = halo_width; iy < ny - halo_width; ++iy) {
    for (int ix = halo_width; ix < nx - halo_width; ++ix) {

      const int ooo = (nx * iy) + ix;

      const int m2o = (nx * iy) + ix - 2;
      const int m1o = (nx * iy) + ix - 1;
      const int p1o = (nx * iy) + ix + 1;
      const int p2o = (nx * iy) + ix + 2;

      const int om2 = (nx * (iy - 2)) + ix;      
      const int om1 = (nx * (iy - 1)) + ix;
      const int op1 = (nx * (iy + 1)) + ix;
      const int op2 = (nx * (iy + 2)) + ix;

      const int f_mo = (nx * iy) + ix;
      const int f_po = (nx * iy) + ix + 1;
      const int f_om = (nx * iy) + ix;
      const int f_op = (nx * iy) + ix + 1;

      const RealType u_face_mo = in_velocity_x[f_mo];
      const RealType u_face_po = in_velocity_x[f_po];
      const RealType u_face_om = in_velocity_y[f_om];
      const RealType u_face_op = in_velocity_y[f_op];

      const RealType in_cell_ooo = in_value[ooo];

      const RealType in_cell_m2o = in_value[m2o];
      const RealType in_cell_m1o = in_value[m1o];
      const RealType in_cell_p1o = in_value[p1o];
      const RealType in_cell_p2o = in_value[p2o];

      const RealType in_cell_om2 = in_value[om2];
      const RealType in_cell_om1 = in_value[om1];
      const RealType in_cell_op1 = in_value[op1];
      const RealType in_cell_op2 = in_value[op2];

      // Reconstructed values at faces (left and right).

      // x direction.
      const RealType s_mo = VanAlbadaLimiter(in_cell_m1o - in_cell_m2o, in_cell_ooo - in_cell_m1o);
      assert((0.0 <= s_mo) && (s_mo <= 1.0));

      const RealType in_cell_lefft_mo = ReconstructValueLefftKappa(in_cell_m1o, s_mo, kappa, dx, in_cell_m1o - in_cell_m2o, in_cell_ooo - in_cell_m1o);
      const RealType in_cell_right_mo = ReconstructValueRightKappa(in_cell_ooo, s_mo, kappa, dx, in_cell_ooo - in_cell_m1o, in_cell_p1o - in_cell_ooo);

      //printf("in_cell_lefft_mo=%lf, in_cell_m1o=%lf, in_cell_ooo=%lf\n", in_cell_lefft_mo, in_cell_m1o, in_cell_ooo);
      assert((min(in_cell_ooo, in_cell_m1o) <= in_cell_lefft_mo) && (in_cell_lefft_mo <= max(in_cell_ooo, in_cell_m1o)));
      
      //printf("in_cell_right_mo=%lf, in_cell_m1o=%lf, in_cell_ooo=%lf\n", in_cell_right_mo, in_cell_m1o, in_cell_ooo);
      assert((min(in_cell_ooo, in_cell_m1o) <= in_cell_right_mo) && (in_cell_right_mo <= max(in_cell_ooo, in_cell_m1o)));

      const RealType s_po = VanAlbadaLimiter(in_cell_ooo - in_cell_m1o, in_cell_p1o - in_cell_ooo);
      assert((0.0 <= s_po) && (s_po <= 1.0));

      const RealType in_cell_lefft_po = ReconstructValueLefftKappa(in_cell_ooo, s_po, kappa, dx, in_cell_ooo - in_cell_m1o, in_cell_p1o - in_cell_ooo);
      const RealType in_cell_right_po = ReconstructValueRightKappa(in_cell_p1o, s_po, kappa, dx, in_cell_p1o - in_cell_ooo, in_cell_p2o - in_cell_p1o);

      //printf("in_cell_lefft_po=%lf, in_cell_p1o=%lf, in_cell_ooo=%lf\n", in_cell_lefft_po, in_cell_p1o, in_cell_ooo);
      assert((min(in_cell_ooo, in_cell_p1o) <= in_cell_lefft_po) && (in_cell_lefft_po <= max(in_cell_ooo, in_cell_p1o)));

      //printf("in_cell_right_po=%lf, in_cell_p1o=%lf, in_cell_ooo=%lf\n", in_cell_right_po, in_cell_p1o, in_cell_ooo);
      assert((min(in_cell_ooo, in_cell_p1o) <= in_cell_right_po) && (in_cell_right_po <= max(in_cell_ooo, in_cell_p1o)));

      // y direction.
      const RealType s_om = VanAlbadaLimiter(in_cell_om1 - in_cell_om2, in_cell_ooo - in_cell_om1);
      const RealType s_op = VanAlbadaLimiter(in_cell_ooo - in_cell_om1, in_cell_op1 - in_cell_ooo);

      const RealType in_cell_lefft_om = ReconstructValueLefftKappa(in_cell_om1, s_om, kappa, dx, in_cell_om1 - in_cell_om2, in_cell_ooo - in_cell_om1);
      const RealType in_cell_right_om = ReconstructValueRightKappa(in_cell_ooo, s_om, kappa, dx, in_cell_ooo - in_cell_om1, in_cell_op1 - in_cell_ooo);

      const RealType in_cell_lefft_op = ReconstructValueLefftKappa(in_cell_ooo, s_op, kappa, dx, in_cell_ooo - in_cell_om1, in_cell_op1 - in_cell_ooo);
      const RealType in_cell_right_op = ReconstructValueRightKappa(in_cell_op1, s_op, kappa, dx, in_cell_op1 - in_cell_ooo, in_cell_op2 - in_cell_op1);

      // Numerical fluxes at faces (simple upwinding for linear advection).
      const RealType flux_face_mo = half * cx * ((fabs(u_face_mo) + u_face_mo) * in_cell_lefft_mo + (u_face_mo - fabs(u_face_mo)) * in_cell_right_mo);
      const RealType flux_face_po = half * cx * ((fabs(u_face_po) + u_face_po) * in_cell_lefft_po + (u_face_po - fabs(u_face_po)) * in_cell_right_po);

      const RealType flux_face_om = half * cy * ((fabs(u_face_om) + u_face_om) * in_cell_lefft_om + (u_face_om - fabs(u_face_om)) * in_cell_right_om);
      const RealType flux_face_op = half * cy * ((fabs(u_face_op) + u_face_op) * in_cell_lefft_op + (u_face_op - fabs(u_face_op)) * in_cell_right_op);

      const RealType out_cell_ooo = in_cell_ooo + dt * (flux_face_mo - flux_face_po + flux_face_om - flux_face_op);

      out_value[ooo] = out_cell_ooo;

    }
  }
    
}

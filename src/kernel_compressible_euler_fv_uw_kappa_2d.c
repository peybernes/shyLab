// -*- c++ -*- (for emacs users)
#include "kernel_compressible_euler_fv_uw_kappa_2d.h"

#define RESTRICT __restrict__

#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#include "kernel_tools.h"

#include "iacaMarks.h"

void CompressibleEulerPhysicalToConservative(int nx, 
					     int ny, 
					     const RealType* RESTRICT in_rho, 
					     const RealType* RESTRICT in_velocity_x,
					     const RealType* RESTRICT in_velocity_y,
					     const RealType* RESTRICT in_total_energy,
					     RealType* RESTRICT out_rho,
					     RealType* RESTRICT out_velocity_x,
					     RealType* RESTRICT out_velocity_y,
					     RealType* RESTRICT out_total_energy) {


  int current_cell = 0;

  for (int iy = 0; iy < ny; ++iy) {

    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);

    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = current_cell;

      out_rho[cell_ooo] = in_rho[cell_ooo];
      out_velocity_x[cell_ooo] = in_rho[cell_ooo] * in_velocity_x[cell_ooo];
      out_velocity_y[cell_ooo] = in_rho[cell_ooo] * in_velocity_y[cell_ooo];
      out_total_energy[cell_ooo] = in_rho[cell_ooo] * in_total_energy[cell_ooo];

      ++current_cell;
    }
  }


}

void CompressibleEulerConservativeToPhysical(int nx, 
					     int ny, 
					     const RealType* RESTRICT out_rho,
					     RealType* RESTRICT out_velocity_x,
					     RealType* RESTRICT out_velocity_y,
					     RealType* RESTRICT out_total_energy) {

  int current_cell = 0;

  for (int iy = 0; iy < ny; ++iy) {

    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);

    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = current_cell;

      out_velocity_x[cell_ooo] /= out_rho[cell_ooo];
      out_velocity_y[cell_ooo] /= out_rho[cell_ooo];
      out_total_energy[cell_ooo] /= out_rho[cell_ooo];

      ++current_cell;
    }
  }

}


void CompressibleEulerFvUwKappa2dX(int nx, 
				   int ny, 
				   int halo_width,
				   RealType dt,
				   RealType dx,
				   RealType dy,
				   RealType kappa,
				   const RealType* RESTRICT in_rho, 
				   const RealType* RESTRICT in_velocity_x,
				   const RealType* RESTRICT in_velocity_y,
				   const RealType* RESTRICT in_total_energy,
				   RealType* RESTRICT out_rho,
				   RealType* RESTRICT out_velocity_x,
				   RealType* RESTRICT out_velocity_y,
				   RealType* RESTRICT out_total_energy) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  const RealType half = 0.5;

  int current_cell = 0;

  for (int iy = 0; iy < ny; ++iy) {
    
    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);

    current_cell += halo_width;

    for (int ix = halo_width; ix < nx - halo_width; ++ix) {

      //IACA_START

      const int cell_m2o = current_cell - 2;
      const int cell_m1o = current_cell - 1;
      const int cell_ooo = current_cell;
      const int cell_p1o = current_cell + 1;
      const int cell_p2o = current_cell + 2;

#include "compressible_euler_2d_data_load.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

      ++current_cell;

    }
    //IACA_END

    current_cell += halo_width;
    
  }

}

void CompressibleEulerFvUwKappa2dY(int nx, 
				   int ny, 
				   int halo_width,
				   RealType dt,
				   RealType dx,
				   RealType dy,
				   RealType kappa,
				   const RealType* RESTRICT in_rho, 
				   const RealType* RESTRICT in_velocity_x,
				   const RealType* RESTRICT in_velocity_y,
				   const RealType* RESTRICT in_total_energy,
				   RealType* RESTRICT out_rho,
				   RealType* RESTRICT out_velocity_x,
				   RealType* RESTRICT out_velocity_y,
				   RealType* RESTRICT out_total_energy) {


  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  const RealType half = 0.5;

  int current_cell = halo_width * nx;

  for (int iy = halo_width; iy < ny - halo_width; ++iy) {

    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);

    for (int ix = 0; ix < nx; ++ix) {

      // const int cell_m2o = (nx * (iy - 2)) + ix;
      // const int cell_m1o = (nx * (iy - 1)) + ix;
      // const int cell_ooo = (nx * (iy + 0)) + ix;
      // const int cell_p1o = (nx * (iy + 1)) + ix;
      // const int cell_p2o = (nx * (iy + 2)) + ix;

      const int cell_m2o = current_cell - (2 * nx);
      const int cell_m1o = current_cell - nx;
      const int cell_ooo = current_cell;
      const int cell_p1o = current_cell + nx;
      const int cell_p2o = current_cell + (2 * nx);

#include "compressible_euler_2d_data_load.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

      ++current_cell;

    }

  }
}



void CompressibleEulerFvUwKappa2dBoundaryConditionsX(int nx, 
						     int ny, 
						     int halo_width,
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     RealType kappa,
						     const RealType* RESTRICT in_rho, 
						     const RealType* RESTRICT in_velocity_x,
						     const RealType* RESTRICT in_velocity_y,
						     const RealType* RESTRICT in_total_energy,
						     RealType* RESTRICT out_rho,
						     RealType* RESTRICT out_velocity_x,
						     RealType* RESTRICT out_velocity_y,
						     RealType* RESTRICT out_total_energy) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  const RealType half = 0.5;

  // Wall boundary conditions. xmin
  for (int iy = 0; iy < ny; ++iy) {
    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);
    for (int ix = 0; ix < halo_width; ++ix) {

      const int cell_m2o = (nx * iy) + ix - 2;
      const int cell_m1o = (nx * iy) + ix - 1;
      const int cell_ooo = (nx * iy) + ix + 0;
      const int cell_p1o = (nx * iy) + ix + 1;
      const int cell_p2o = (nx * iy) + ix + 2;

      const int index_in_main_direction = ix;

#include "compressible_euler_2d_data_load_mirror_min.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

    }
  }
  
  // Wall boundary conditions. xmax
  for (int iy = 0; iy < ny; ++iy) {
    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);
    for (int ix = nx - halo_width; ix < nx; ++ix) {

      const int cell_m2o = (nx * iy) + ix - 2;
      const int cell_m1o = (nx * iy) + ix - 1;
      const int cell_ooo = (nx * iy) + ix + 0;
      const int cell_p1o = (nx * iy) + ix + 1;
      const int cell_p2o = (nx * iy) + ix + 2;

      const int index_in_main_direction = ix;
      const int size_in_main_direction = nx;

#include "compressible_euler_2d_data_load_mirror_max.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

    }
  }

//   // Periodical boundary conditions. xmin
//   for (int iy = 0; iy < ny; ++iy) {
//     for (int ix = 0; ix < halo_width; ++ix) {

//       const int cell_m2o = (nx * iy) + ix + nx - 2;
//       const int cell_m1o = (ix == 0 ? (nx * iy) + ix + nx - 1 : (nx * iy) + ix - 1);
//       const int cell_ooo = (nx * iy) + ix + 0;
//       const int cell_p1o = (nx * iy) + ix + 1;
//       const int cell_p2o = (nx * iy) + ix + 2;

// #include "compressible_euler_2d_data_load.h"
// #include "compressible_euler_2d_computation.h"

//       out_rho[cell_ooo] += out_rho_ooo;
//       out_velocity_x[cell_ooo] += out_velocity_x_ooo;
//       out_velocity_y[cell_ooo] += out_velocity_y_ooo;
//       out_total_energy[cell_ooo] += out_total_energy_ooo;

//     }
//   }

//   // Periodical boundary conditions. xmax
//   for (int iy = 0; iy < ny; ++iy) {
//     for (int ix = nx - halo_width; ix < nx; ++ix) {

//       const int cell_m2o = (nx * iy) + ix - 2;
//       const int cell_m1o = (nx * iy) + ix - 1;
//       const int cell_ooo = (nx * iy) + ix + 0;
//       const int cell_p1o = (ix == nx - 1 ? (nx * iy) + ix - nx + 1 : (nx * iy) + ix + 1);
//       const int cell_p2o = (nx * iy) + ix - nx + 2;

// #include "compressible_euler_2d_data_load.h"
// #include "compressible_euler_2d_computation.h"

//       out_rho[cell_ooo] += out_rho_ooo;
//       out_velocity_x[cell_ooo] += out_velocity_x_ooo;
//       out_velocity_y[cell_ooo] += out_velocity_y_ooo;
//       out_total_energy[cell_ooo] += out_total_energy_ooo;

//     }
//   }
//  exit(1);

}

void CompressibleEulerFvUwKappa2dBoundaryConditionsY(int nx, 
						     int ny, 
						     int halo_width,
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     RealType kappa,
						     const RealType* RESTRICT in_rho, 
						     const RealType* RESTRICT in_velocity_x,
						     const RealType* RESTRICT in_velocity_y,
						     const RealType* RESTRICT in_total_energy,
						     RealType* RESTRICT out_rho,
						     RealType* RESTRICT out_velocity_x,
						     RealType* RESTRICT out_velocity_y,
						     RealType* RESTRICT out_total_energy) {

  assert(0 < dx);
  assert(0 < dy);
  assert(0 <= dt);

  const RealType cx = 1.0 / dx;
  const RealType cy = 1.0 / dy;

  const RealType half = 0.5;

  // Wall boundary conditions. ymin
  for (int iy = 0; iy < halo_width; ++iy) {
    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_m2o = (nx * (iy - 2)) + ix;
      const int cell_m1o = (nx * (iy - 1)) + ix;
      const int cell_ooo = (nx * (iy + 0)) + ix;
      const int cell_p1o = (nx * (iy + 1)) + ix;
      const int cell_p2o = (nx * (iy + 2)) + ix;

      const int index_in_main_direction = iy;

#include "compressible_euler_2d_data_load_mirror_min.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

    }
  }

  // Wall boundary conditions. ymax
  for (int iy = ny - halo_width; iy < ny; ++iy) {
    ASSUME_ALIGNED(in_rho);
    ASSUME_ALIGNED(in_velocity_x);
    ASSUME_ALIGNED(in_velocity_y);
    ASSUME_ALIGNED(in_total_energy);
    ASSUME_ALIGNED(out_rho);
    ASSUME_ALIGNED(out_velocity_x);
    ASSUME_ALIGNED(out_velocity_y);
    ASSUME_ALIGNED(out_total_energy);
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_m2o = (nx * (iy - 2)) + ix;
      const int cell_m1o = (nx * (iy - 1)) + ix;
      const int cell_ooo = (nx * (iy + 0)) + ix;
      const int cell_p1o = (nx * (iy + 1)) + ix;
      const int cell_p2o = (nx * (iy + 2)) + ix;

      const int index_in_main_direction = iy;
      const int size_in_main_direction = ny;

#include "compressible_euler_2d_data_load_mirror_max.h"
#include "compressible_euler_2d_computation.h"

      out_rho[cell_ooo] += out_rho_ooo;
      out_velocity_x[cell_ooo] += out_velocity_x_ooo;
      out_velocity_y[cell_ooo] += out_velocity_y_ooo;
      out_total_energy[cell_ooo] += out_total_energy_ooo;

    }
  }


//  // Periodical boundary conditions. ymin
//   for (int iy = 0; iy < halo_width; ++iy) {
//     for (int ix = 0; ix < nx; ++ix) {

//       const int cell_m2o = (nx * (iy + ny - 2)) + ix;
//       const int cell_m1o = (iy == 0 ? (nx * (iy + ny - 1)) + ix : (nx * (iy - 1)) + ix);
//       const int cell_ooo = (nx * (iy + 0)) + ix;
//       const int cell_p1o = (nx * (iy + 1)) + ix;
//       const int cell_p2o = (nx * (iy + 2)) + ix;

// #include "compressible_euler_2d_data_load.h"
// #include "compressible_euler_2d_computation.h"

//       out_rho[cell_ooo] += out_rho_ooo;
//       out_velocity_x[cell_ooo] += out_velocity_x_ooo;
//       out_velocity_y[cell_ooo] += out_velocity_y_ooo;
//       out_total_energy[cell_ooo] += out_total_energy_ooo;

//     }
//   }

//   // Periodical boundary conditions. ymax
//   for (int iy = ny - halo_width; iy < ny; ++iy) {
//     for (int ix = 0; ix < nx; ++ix) {

//       const int cell_m2o = (nx * (iy - 2)) + ix;
//       const int cell_m1o = (nx * (iy - 1)) + ix;
//       const int cell_ooo = (nx * (iy + 0)) + ix;
//       const int cell_p1o = (iy == ny - 1 ? (nx * (iy - ny + 1)) + ix : (nx * (iy + 1)) + ix);
//       const int cell_p2o = (nx * (iy - ny + 2)) + ix;

// #include "compressible_euler_2d_data_load.h"
// #include "compressible_euler_2d_computation.h"

//       out_rho[cell_ooo] += out_rho_ooo;
//       out_velocity_x[cell_ooo] += out_velocity_x_ooo;
//       out_velocity_y[cell_ooo] += out_velocity_y_ooo;
//       out_total_energy[cell_ooo] += out_total_energy_ooo;

//     }
//   }

}


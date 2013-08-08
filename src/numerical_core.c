// -*- c++ -*- (for emacs users)
/// \file numerical_core.c Numerical discretization of the
/// nonlinear shallow water equation using the finite volumes method
/// framework.
///
/// Conventions specific to this file :
///
/// - All input POD values (e.g. single integer or float/double) are
/// passed by value (and, consequently, copied).
/// - All input arrays are const'd and restricted.
/// - All output arrays are restricted.
/// - All output POD values (e.g. single integer or float/double) are 
/// passed by pointer and restricted.
///
/// Moreover, comments emphasize what are the inputs/outputs, and how are
/// memory access patterns (loads and stores).
///
/// Keep in mind that, among all these functions (and neglecting I/O):
/// 
/// - NumericalFluxInterior accounts for 65% to 70% of total
/// simulation time.
/// - CellResiduals accounts for 15% to 20% of total simulation time.
/// - UpdateValuesRK1 accounts for approx 10% of simulation time.
///
/// These are typical numbers, use your favorite profiler to learn more..

#include "numerical_core.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <sys/time.h>
#include <values.h>

#define ASSERT assert

typedef enum _BcType {WALL, HEIGHTSUBC, FLOWSUBC} BcType;

double WallClock() {

  struct timeval t;
  gettimeofday(&t, NULL);

  return t.tv_sec * 1.0e3 + t.tv_usec * 1.0e-3;
}

// HMPP : creation of the grouplet
#pragma hmpp <shy> group, target=CUDA

// HMPP :  maps arguments of different codelet in the same memory space

// Mesh values. Does not change after simulation initialization.
#pragma hmpp <shy> mapbyname face_front_cell_id, face_back_cell_id
#pragma hmpp <shy> mapbyname normal_x, normal_y
#pragma hmpp <shy> mapbyname face_volume, cell_volume
#pragma hmpp <shy> mapbyname cell_signed_face1_id, cell_signed_face2_id, cell_signed_face3_id

// Cell values.
#pragma hmpp <shy> mapbyname in_H, in_U, in_V, in_Z, 
#pragma hmpp <shy> mapbyname out_H, out_U, out_V
#pragma hmpp <shy> mapbyname cell_timestep

// Boundary face values.
#pragma hmpp <shy> mapbyname HL_face, HR_face, UL_face, UR_face, VL_face, VR_face, ZL_face, ZR_face

// Face fluxes.
#pragma hmpp <shy> mapbyname flux_H, flux_HU, flux_HV, flux_bathy_left, flux_bathy_right, face_max_eigenvalue

#pragma hmpp <shy> facefromcell codelet, args[0-11].io=in, args[12-19].io=out, args[4-7].size={nb_cells}, args[8-11].size={nb_faces}, args[12-19].size={nb_boundary_faces}, args[0-3].const=true
void InterpolateCellToFace(int id_face_begin,
			   int id_face_end,
			   const int* __restrict face_back_cell_id,
			   const int* __restrict face_front_cell_id,
			   const RealType* __restrict in_rho, 
			   const RealType* __restrict in_u,
			   const RealType* __restrict in_v, 
			   const RealType* __restrict in_w,
			   const RealType* __restrict in_e,
			   const RealType* __restrict normal_x, 
			   const RealType* __restrict normal_y,
			   const RealType* __restrict normal_z,
			   RealType* __restrict rhoL_face,  
			   RealType* __restrict uL_face,
			   RealType* __restrict vL_face,  
			   RealType* __restrict wL_face,
			   RealType* __restrict eL_face,
			   RealType* __restrict rhoR_face,  
			   RealType* __restrict uR_face,
			   RealType* __restrict vR_face,  
			   RealType* __restrict wR_face,
			   RealType* __restrict eR_face) {

  for (int i = id_face_begin; i < id_face_end; ++i) {

    const int front_cell_id = face_front_cell_id[i];
    
    // Order 1 reconstruction. The face values are taken equal to the
    // corresponding cell center values.

    rhoL_face[i] = in_rho[front_cell_id];
    uL_face[i] = in_u[front_cell_id];
    vL_face[i] = in_v[front_cell_id];
    wL_face[i] = in_w[front_cell_id];
    eL_face[i] = in_e[front_cell_id];
    
    const RealType nx = normal_x[i];
    const RealType ny = normal_y[i];
    const RealType nz = normal_z[i];
      
    // Wall boundary conditions

    const RealType in_normal_velocity = 
      (in_u[front_cell_id] * nx) +
      (in_v[front_cell_id] * ny) +
      (in_w[front_cell_id] * nz);

    const RealType in_tangent_velocity =
      - (in_u[front_cell_id] * ny) +
      (in_v[front_cell_id] * nx);

    RealType out_normal_velocity = 0.0;
    RealType out_tangent_velocity = 0.0;

  }

}

void NumericalFluxCellToFace(int id_face_begin,
			     int id_face_end,
			     const int* __restrict face_back_cell_id,
			     const int* __restrict face_front_cell_id,
			     const RealType* __restrict face_areas,
			     const RealType* __restrict in_rho,
			     const RealType* __restrict in_u,
			     const RealType* __restrict in_v,
			     const RealType* __restrict in_w,
			     const RealType* __restrict in_e,
			     const RealType* __restrict normal_x, 
			     const RealType* __restrict normal_y,
			     const RealType* __restrict normal_z,
			     RealType* __restrict flux_rho, 
			     RealType* __restrict flux_rho_u, 
			     RealType* __restrict flux_rho_v,
			     RealType* __restrict flux_rho_w,
			     RealType* __restrict flux_rho_e,
			     RealType* __restrict face_max_eigen) {
  
#pragma omp parallel for
  for (int i = id_face_begin; i < id_face_end; ++i) {

    /// Beginning of memory loads.

    // Memory loads -- contiguous.

    const int back_cell_id = face_back_cell_id[i];
    const int front_cell_id = face_front_cell_id[i];

    ASSERT((0 <= back_cell_id));
    ASSERT((0 <= front_cell_id));

    const RealType face_area = face_areas[i];

    const RealType nx = normal_x[i];
    const RealType ny = normal_y[i];
    const RealType nz = normal_y[i];

    // Memory loads -- not contiguous.

    const RealType rhoL = in_rho[front_cell_id];
    const RealType uL = in_u[front_cell_id];
    const RealType vL = in_v[front_cell_id];
    const RealType wL = in_w[front_cell_id];
    const RealType eL = in_e[front_cell_id];

    const RealType rhoR = in_rho[back_cell_id];
    const RealType uR = in_u[back_cell_id];
    const RealType vR = in_v[back_cell_id];
    const RealType wR = in_w[back_cell_id];
    const RealType eR = in_e[back_cell_id];

    /// End of memory loads.

    /// Beginning of computations and memory stores (contiguous).
#include "euler_flux_hll_inl.h"
    /// End of computations and memory stores (contiguous).
  }
}


void NumericalFluxFaceToFace(int id_face_begin,
			     int id_face_end,
			     const RealType* __restrict face_areas,
			     const RealType* __restrict rhoL_face,
			     const RealType* __restrict uL_face,
			     const RealType* __restrict vL_face,
			     const RealType* __restrict wL_face,
			     const RealType* __restrict eL_face, 
			     const RealType* __restrict rhoR_face, 
			     const RealType* __restrict uR_face, 
			     const RealType* __restrict vR_face,
			     const RealType* __restrict wR_face,
			     const RealType* __restrict eR_face,
			     const RealType* __restrict normal_x, 
			     const RealType* __restrict normal_y,
			     const RealType* __restrict normal_z,
			     RealType* __restrict flux_rho, 
			     RealType* __restrict flux_rho_u, 
			     RealType* __restrict flux_rho_v,
			     RealType* __restrict flux_rho_w,
			     RealType* __restrict flux_rho_e,
			     RealType* __restrict face_max_eigen) {

#pragma omp for
  for (int i = id_face_begin; i < id_face_end; ++i) {

    /// Beginning of memory loads -- contiguous.

    const RealType face_area = face_areas[i];    
    
    const RealType nx = normal_x[i];
    const RealType ny = normal_y[i];
    const RealType nz = normal_z[i];

    const RealType rhoL = rhoL_face[i];
    const RealType uL = uL_face[i];
    const RealType vL = vL_face[i];
    const RealType wL = wL_face[i];
    const RealType eL = eL_face[i];

    const RealType rhoR = rhoR_face[i];
    const RealType uR = uR_face[i];
    const RealType vR = vR_face[i];
    const RealType wR = wR_face[i];
    const RealType eR = eR_face[i];

    /// End of memory loads -- contiguous.

    /// Beginning of computations
    #include "euler_flux_hll_inl.h"
    /// End of computations.
  }
}

void ComputeTimestep(int nb_cells,
		     RealType CFL, 
		     RealType max_timestep,
		     RealType* __restrict cell_timestep,
		     RealType* __restrict out_timestep_ptr) {

  RealType timestep = max_timestep;

#pragma hmppcg parallel, reduce (min:timestep)
#pragma omp reduce(min:timestep)
  for (int i = 0; i < nb_cells; ++i)
    timestep = fminf(timestep, cell_timestep[i]);
  
  *out_timestep_ptr = CFL * timestep;
  
}

#pragma hmpp <shy> cell_residuals codelet, args[0-14].io=in, args[15-18].io=out, args[3-10].size={nb_faces}, args[11-18].size={nb_cells}, args[0-2].const=true
void CellResiduals(int id_cell_begin, 
		   int id_cell_end,
		   int nb_faces, 
		   const int* __restrict cell_signed_face1_id,
		   const int* __restrict cell_signed_face2_id,
		   const int* __restrict cell_signed_face3_id,
		   const int* __restrict cell_signed_face4_id,
		   const RealType* __restrict flux_rho,
		   const RealType* __restrict flux_rho_u, 
		   const RealType* __restrict flux_rho_v,
		   const RealType* __restrict flux_rho_w,
		   const RealType* __restrict flux_rho_e,
		   const RealType* __restrict face_max_eigenvalue,
		   const RealType* __restrict normal_x,
		   const RealType* __restrict normal_y,
		   const RealType* __restrict normal_z,
		   const RealType* __restrict cell_volume,
		   RealType* __restrict out_rho, 
		   RealType* __restrict out_rho_u, 
		   RealType* __restrict out_rho_v,
		   RealType* __restrict out_rho_w,
		   RealType* __restrict out_rho_e,
		   RealType* __restrict cell_timestep) {

#pragma hmppcg parallel
#pragma omp parallel for
  for (int i = id_cell_begin; i < id_cell_end; ++i) {

    const int face1_id = abs(cell_signed_face1_id[i]);
    const int face2_id = abs(cell_signed_face2_id[i]);
    const int face3_id = abs(cell_signed_face3_id[i]);
    const int face4_id = abs(cell_signed_face4_id[i]);

    ASSERT((1 <= face1_id) && (face1_id <= nb_faces));
    ASSERT((1 <= face2_id) && (face2_id <= nb_faces));
    ASSERT((1 <= face3_id) && (face3_id <= nb_faces));
    ASSERT((1 <= face4_id) && (face4_id <= nb_faces));

    const int sign1 = 
      (cell_signed_face1_id[i] - 1 >= 0) - (cell_signed_face1_id[i] - 1 < 0);
    const int sign2 = 
      (cell_signed_face2_id[i] - 1 >= 0) - (cell_signed_face2_id[i] - 1 < 0);
    const int sign3 = 
      (cell_signed_face3_id[i] - 1 >= 0) - (cell_signed_face3_id[i] - 1 < 0);
    const int sign4 = 
      (cell_signed_face4_id[i] - 1 >= 0) - (cell_signed_face4_id[i] - 1 < 0);

    ASSERT((sign1 == 1) || (sign1 == -1));
    ASSERT((sign2 == 1) || (sign2 == -1));
    ASSERT((sign3 == 1) || (sign3 == -1));
    ASSERT((sign4 == 1) || (sign4 == -1));

    const RealType nx1 = normal_x[face1_id - 1];
    const RealType ny1 = normal_y[face1_id - 1];
    const RealType nz1 = normal_z[face1_id - 1];
    const RealType flux_rho_1 = flux_rho[face1_id - 1];
    const RealType flux_rho_u_1 = flux_rho_u[face1_id - 1];
    const RealType flux_rho_v_1 = flux_rho_v[face1_id - 1];
    const RealType flux_rho_w_1 = flux_rho_w[face1_id - 1];
    const RealType flux_rho_e_1 = flux_rho_e[face1_id - 1];
    const RealType max_eigen_1 = face_max_eigenvalue[face1_id - 1];

    const RealType nx2 = normal_x[face2_id - 1];
    const RealType ny2 = normal_y[face2_id - 1];
    const RealType nz2 = normal_z[face2_id - 1];
    const RealType flux_rho_2 = flux_rho[face2_id - 1];
    const RealType flux_rho_u_2 = flux_rho_u[face2_id - 1];
    const RealType flux_rho_v_2 = flux_rho_v[face2_id - 1];
    const RealType flux_rho_w_2 = flux_rho_w[face2_id - 1];
    const RealType flux_rho_e_2 = flux_rho_e[face2_id - 1];
    const RealType max_eigen_2 = face_max_eigenvalue[face2_id - 1];

    const RealType nx3 = normal_x[face3_id - 1];
    const RealType ny3 = normal_y[face3_id - 1];
    const RealType nz3 = normal_z[face3_id - 1];
    const RealType flux_rho_3 = flux_rho[face3_id - 1];
    const RealType flux_rho_u_3 = flux_rho_u[face3_id - 1];
    const RealType flux_rho_v_3 = flux_rho_v[face3_id - 1];
    const RealType flux_rho_w_3 = flux_rho_w[face3_id - 1];
    const RealType flux_rho_e_3 = flux_rho_e[face3_id - 1];
    const RealType max_eigen_3 = face_max_eigenvalue[face3_id - 1];

    const RealType nx4 = normal_x[face4_id - 1];
    const RealType ny4 = normal_y[face4_id - 1];
    const RealType nz4 = normal_z[face4_id - 1];
    const RealType flux_rho_4 = flux_rho[face4_id - 1];
    const RealType flux_rho_u_4 = flux_rho_u[face4_id - 1];
    const RealType flux_rho_v_4 = flux_rho_v[face4_id - 1];
    const RealType flux_rho_w_4 = flux_rho_v[face4_id - 1];
    const RealType flux_rho_e_4 = flux_rho_e[face4_id - 1];
    const RealType max_eigen_4 = face_max_eigenvalue[face4_id -1];

    const RealType volume = cell_volume[i];

    const RealType max_timestep = 1.0;

    cell_timestep[i] = 2 * volume /
      (max_eigen_1 + max_eigen_2 + max_eigen_3 + max_eigen_4);

    cell_timestep[i] = MIN(cell_timestep[i], max_timestep);

    out_rho[i] = volume *
      (sign1 * flux_rho_1 + sign2 * flux_rho_2 +
       sign3 * flux_rho_3 + sign4 * flux_rho_4);

    out_rho_u[i] = volume *
      (sign1 * flux_rho_u_1 + sign2 * flux_rho_u_2 +
       sign3 * flux_rho_u_3 + sign4 * flux_rho_u_4);
    
    out_rho_v[i] = volume *
      (sign1 * flux_rho_v_1 + sign2 * flux_rho_v_2 +
       sign3 * flux_rho_v_3 + sign4 * flux_rho_v_4);

    out_rho_v[i] = volume *
      (sign1 * flux_rho_w_1 + sign2 * flux_rho_w_2 +
       sign3 * flux_rho_w_3 + sign4 * flux_rho_w_4);

    out_rho_e[i] = volume *
      (sign1 * flux_rho_e_1 + sign2 * flux_rho_e_2 +
       sign3 * flux_rho_e_3 + sign4 * flux_rho_e_4);


  }
}
  
/// Update simulation values using forward Euler time integration.
///
/// \remark Memory loads are contiguous.
/// \remark Memory stores are contiguous.
///
/// \param[in] nb_cells Number of mesh cells.
/// \param[in] EPS smallness threshold for water height.
/// \param[in] CFL Courant-Friedrichs-Levy stability condition number.
/// \param[in] maximal_timestep Maximum possible timestep.
/// \param[in] in_cells_H Water height at the beginning of time step.
/// \param[in] in_cells_U Water x velocity at the beginning of time
/// step.
/// \param[in] in_cells_V Water y velocity at the beginning of time
/// step.
/// \param[in] in_cell_volume mesh cell volumes (triangle areas).
///
/// \param[out] out_cells_H Water height at the end of time step.
/// \param[out] out_cells_U Water x velocity at the end of time step.
/// \param[out] out_cells_V Water y velocity at the end of time step.
#pragma hmpp <shy> RK1 codelet, args[6-8].io=inout, args[3-8].size={nb_cells}, args[0-2].const=true
void EvolveValuesRK1(int nb_cells,
		     RealType EPS, 
		     RealType dt,
		     const RealType* __restrict in_rho,
		     const RealType* __restrict in_u,
		     const RealType* __restrict in_v,
		     const RealType* __restrict in_w,
		     const RealType* __restrict in_e,
		     RealType* __restrict out_rho,
		     RealType* __restrict out_u,
		     RealType* __restrict out_v,
		     RealType* __restrict out_w,
		     RealType* __restrict out_e) {

#pragma hmppcg grid blocksize "64x1" 
#pragma omp parallel for
  for (int i = 0; i < nb_cells; ++i) {

    //// Beginning of memory loads -- contiguous.
    const RealType in_cells_rho = in_rho[i];
    const RealType in_cells_u = in_u[i];
    const RealType in_cells_v = in_v[i];
    const RealType in_cells_w = in_w[i];
    const RealType in_cells_e = in_e[i];
    
    const RealType residual_rho = out_rho[i];
    const RealType residual_rho_u = out_u[i];
    const RealType residual_rho_v = out_v[i];
    const RealType residual_rho_w = out_w[i];
    const RealType residual_rho_e = out_e[i];
   
    const RealType out_cells_rho = in_cells_rho + (dt * residual_rho);

    const RealType out_cells_rho_u = 
      (in_cells_rho * in_cells_u) + (dt * residual_rho_u);

    const RealType out_cells_rho_v = 
      (in_cells_rho * in_cells_v) + (dt * residual_rho_v);

    const RealType out_cells_rho_w = 
      (in_cells_rho * in_cells_w) + (dt * residual_rho_w);

    const RealType out_cells_rho_e = 
      (in_cells_rho * in_cells_e) + (dt * residual_rho_e);
 
    //// Beginning of memory stores.

    out_rho[i] = out_cells_rho;

    out_u[i] = (out_cells_rho <= EPS ? 0.0 : out_cells_rho_u / out_cells_rho);
    out_v[i] = (out_cells_rho <= EPS ? 0.0 : out_cells_rho_v / out_cells_rho);
    out_w[i] = (out_cells_rho <= EPS ? 0.0 : out_cells_rho_w / out_cells_rho);

    out_e[i] = (out_cells_rho <= EPS ? 0.0 : out_cells_rho_e / out_cells_rho);

    //// End of memory stores.
  }
}

/// Copy simulation values.
/// 
/// \remark Memory loads are contiguous.
/// \remark Memory stores are contiguous.
///
/// \param[in] nb_cells Number of mesh cells ; size of values.
/// \param[in] src_H Input water height.
/// \param[in] src_U Input x velocity.
/// \param[in] src_H Input y velocity.
///
/// \param[out] dst_H Output water height.
/// \param[out] dst_U Output x velocity.
/// \param[out] dst_V Output y velocity.
void CopyValues(int nb_cells,
		const RealType* __restrict src_rho,
		const RealType* __restrict src_u,
		const RealType* __restrict src_v,
		const RealType* __restrict src_w,
		const RealType* __restrict src_e,
		RealType* __restrict dst_rho,
		RealType* __restrict dst_u,
		RealType* __restrict dst_v,
		RealType* __restrict dst_w,
		RealType* __restrict dst_e) {

#pragma omp for
  for (int i = 0; i < nb_cells; ++i) {
    
    dst_rho[i] = src_rho[i];
    dst_u[i] = src_u[i];
    dst_v[i] = src_v[i];
    dst_w[i] = src_w[i];
    dst_e[i] = src_e[i];
    
  }
}

void HMPPAllocate(int nb_boundary_faces, int nb_faces, int nb_cells)
{
#pragma hmpp <shy> allocate
}

void HMPPRelease()
{
#pragma hmpp <shy> release
}

#ifndef KERNELS_HPP_
#define KERNELS_HPP_

#ifdef __cplusplus
extern "C" {
#endif

// Finite volumes upwind methods.
#include "kernel_advection_fv_uw_1_2d.h"
#include "kernel_advection_fv_uw_1_3d.h"
#include "kernel_advection_fv_uw_2_2d.h"
#include "kernel_advection_fv_uw_2_3d.h"
#include "kernel_advection_fv_uw_kappa_2d.h"
#include "kernel_advection_fv_uw_kappa_3d.h"
#include "kernel_advection_fv_central_compact_explicit_2d.h"
#include "kernel_compressible_euler_fv_uw_kappa_2d.h"

// Lagrange Remap kernels
#include "kernel_ad_projection_2d.h"
#include "kernel_direct_projection_2d.h"
#include "kernel_lagrange_2d.h"
#include "kernel_ad_projection_2d_boundary_conditions.h"

#ifdef __cplusplus
}
#endif
#include "kernel_ad_projection_2d_boundary_interface.hpp"
#include "lagrange_2d_interface.hpp"
#include "ad_projection_2d_interface.hpp"
#include "direct_projection_2d_interface.hpp"

#endif // KERNELS_HPP_

#ifndef KERNEL_TOOLS_H
#define KERNEL_TOOLS_H

#define RESTRICT __restrict__

#define ALIGN_BYTES 64

#ifdef __GNUC__
#define ASSUME_ALIGNED(x) __builtin_assume_aligned(x, ALIGN_BYTES)
#elif __INTEL_COMPILER
#define ASSUME_ALIGNED(x) __assume_aligned(x, ALIGN_BYTES)
#else
#define ASSUME_ALIGNED(x) assert(0)
#endif


#ifdef __cplusplus
extern "C" {
#endif

static inline RealType EquationOfStatePerfectGas(RealType gamma, RealType rho, RealType ux, RealType uy, RealType total_energy) {

  const RealType one = 1.0;
  const RealType half = 0.5;

  return (gamma - one) * rho * (total_energy - half * (ux * ux + uy * uy));

}

/// Van Albada limiter (reference : Nishikawa 2008, carbuncle free solver)
static inline RealType VanAlbadaLimiter(RealType a, RealType b) {
  
  const RealType epsilon = 1.0e-6;
  const RealType two = 2.0;

  return (two * ((a * b) + epsilon) / ((a * a) + (b * b) + (two * epsilon)));
}

/// Van Leer kappa scheme
///
/// kappa = 1/3 : cubic (order 3) reconstruction (but not necessary best)
/// kappa = 1/2 : QUICK : quadratic reconstruction
/// kappa = 0 : Fromm scheme
/// kappa = 1 : central difference scheme (unstable)
/// kappa = -1 : linear upwind scheme
static inline RealType ReconstructValueLefftKappa(RealType in, RealType s, RealType kappa, RealType dx, RealType delta_minus, RealType delta_plus) {

  const RealType one_quarter = 0.25;
  const RealType one = 1.0;

  return (in + one_quarter * s * dx * (((one - kappa * s) * delta_minus) + ((one + kappa * s) * delta_plus)));

}  

static inline RealType ReconstructValueRightKappa(RealType in, RealType s, RealType kappa, RealType dx, RealType delta_minus, RealType delta_plus) {

  const RealType one_quarter = 0.25;
  const RealType one = 1.0;

  return (in - one_quarter * s * dx * (((one + kappa * s) * delta_minus + (one - kappa * s) * delta_plus)));

}  

#ifdef __cplusplus
}
#endif

#endif // KERNEL_TOOLS_H

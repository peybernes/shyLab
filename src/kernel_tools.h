#ifndef KERNEL_TOOLS_H
#define KERNEL_TOOLS_H

#define RESTRICT __restrict__

#ifdef __cplusplus
extern "C" {
#endif

static inline RealType EquationOfStatePerfectGas(RealType gamma, RealType rho, RealType ux, RealType uy, RealType total_energy) {
  
  return (1.4 - 1.0) * rho * (total_energy - 0.5 * (ux * ux + uy * uy));

}

/// Van Albada limiter (reference : Nishikawa 2008, carbuncle free solver)
static inline RealType VanAlbadaLimiter(RealType a, RealType b) {

  //return 0.0;
  return (2.0 * ((a * b) + 1.0e-6) / ((a * a) + (b * b) + (2.0 * 1.0e-6)));

}

/// Van Leer kappa scheme
///
/// kappa = 1/3 : cubic (order 3) reconstruction (but not necessary best)
/// kappa = 1/2 : QUICK : quadratic reconstruction
/// kappa = 0 : Fromm scheme
/// kappa = 1 : central difference scheme (unstable)
/// kappa = -1 : linear upwind scheme
static inline RealType ReconstructValueLefftKappa(RealType in, RealType s, RealType kappa, RealType dx, RealType delta_minus, RealType delta_plus) {

  return (in + 0.25 * s * dx * (((1.0 - kappa * s) * delta_minus) + ((1.0 + kappa * s) * delta_plus)));

}  

static inline RealType ReconstructValueRightKappa(RealType in, RealType s, RealType kappa, RealType dx, RealType delta_minus, RealType delta_plus) {

  return (in - 0.25 * s * dx * (((1.0 + kappa * s) * delta_minus + (1.0 - kappa * s) * delta_plus)));

}  

#ifdef __cplusplus
}
#endif

#endif // KERNEL_TOOLS_H

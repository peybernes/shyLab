#ifndef KERNEL_TOOLS_H
#define KERNEL_TOOLS_H

#include <cmath>
#include <xmmintrin.h>
#include <algorithm> 

#if !defined(SHY_ASM_COMMENT)
  #if (defined __GNUC__) && ( defined(__i386__) || defined(__x86_64__) )
    #define SHY_ASM_COMMENT(X)  asm("#" X)
  #else
    #define SHY_ASM_COMMENT(X)
  #endif
#endif

#if !defined(SHY_COMPILER_BARRIER)
  #if (defined __GNUC__) && ( defined(__i386__) || defined(__x86_64__) )
    #define SHY_COMPILER_BARRIER(X)  asm volatile("" ::: "memory")
  #else
    #define SHY_COMPILER_BARRIER(X)
  #endif
#endif


#define RESTRICT __restrict

#define ALIGN_BYTES 64

/* #ifdef __GNUC__ */
/* #define ASSUME_ALIGNED(x) __builtin_assume_aligned(x, ALIGN_BYTES) */
/* #elif __INTEL_COMPILER */
#define ASSUME_ALIGNED(x) __assume_aligned(x, ALIGN_BYTES)
/* #else */
/* #define ASSUME_ALIGNED(x) assert(0) */
/* #endif */


#ifdef __cplusplus
extern "C" {
#endif

  static inline RealType EquationOfState(RealType gamma, RealType rho, RealType e, RealType pi) {
    
    const RealType one = 1.0;
    
    return  (gamma - one) * rho * e - pi;
  
  }
 

  static inline RealType EquationOfStatePerfectGas(RealType gamma, RealType rho, RealType ux, RealType uy, RealType total_energy) {

    const RealType one = 1.0;
    const RealType half = 0.5;
    
    return (gamma - one) * rho * (total_energy - half * (ux * ux + uy * uy));
    
  }
 

  static inline RealType SpeedOfSound(RealType gamma, RealType rho, RealType p, RealType pi) {
        
    return  sqrt((gamma * p + pi) / rho);
  
  }


  static inline RealType MySign(RealType x) {

    return (x > 0.0) - (x < 0.0);

  }


  static inline RealType MinmodLimiter(RealType a, RealType b) {

    return (a < 0 ? -1.0 : 1.0 ) * (a * b > 0) * std::min(fabs(a),fabs(b));

  }

  /// Van Albada limiter (reference : Nishikawa 2008, carbuncle free solver)
  static inline RealType VanAlbadaLimiter(RealType a, RealType b) {
  
    const RealType epsilon = 1.0e-14;
    const RealType two = 2.0;
    
    return ( ((a * b) * (a + b )) * ((a * a + b * b) > 0.) / (a * a + b * b + epsilon) );
	//return 0.;
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

    // return (in + one_quarter * s * dx * (((one - kappa * s) * delta_minus) + ((one + kappa * s) * delta_plus)));
    return 0.0;
  }  

  static inline RealType ReconstructValueRightKappa(RealType in, RealType s, RealType kappa, RealType dx, RealType delta_minus, RealType delta_plus) {
    
    const RealType one_quarter = 0.25;
    const RealType one = 1.0;
    
    return (in - one_quarter * s * dx * (((one + kappa * s) * delta_minus + (one - kappa * s) * delta_plus)));

  }  
  
  static double diff(timespec start, timespec end) {

    timespec result;
    
    if ((end.tv_nsec-start.tv_nsec) < 0) {
      
      result.tv_sec = end.tv_sec-start.tv_sec - 1;
      result.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
      
    } else {
      
      result.tv_sec = end.tv_sec - start.tv_sec;
      result.tv_nsec = end.tv_nsec - start.tv_nsec;

    }

    const double  ms_time = 
      (1000.0 * result.tv_sec) + (result.tv_nsec / 1000000.0);

    return ms_time;
  }
  
#ifdef __cplusplus
}
#endif

#endif // KERNEL_TOOLS_H

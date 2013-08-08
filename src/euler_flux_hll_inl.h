// -*- c++ -*- (for emacs users).

// Advection equation HLL flux.
// This inline file needs :
// rhoL, rhoR, uL, uR, vL, vR, wL, wR, nx, ny, nz, face_area, gravity

const RealType EPS = 1.0e-10;

const RealType zero = 0.0;
const RealType one = 1.0;

// Speed of sound.
const RealType cL = 0.0;
const RealType cR = 0.0;

// Intermediate velocities, for use in the approximate Riemann
// problem.
const RealType uLn = (uL * nx) + (vL * ny) + (wL * nz);
const RealType uRn = (uR * nx) + (vR * ny) + (wR * nz);

const RealType max_left_eigenvalue =
  MAX(abs(uLn - cL), abs(uLn + cL));
    
const RealType max_right_eigenvalue =
  MAX(abs(uRn - cR), abs(uRn + cR));
    
const RealType max_eigenvalue = 
  MAX(max_left_eigenvalue, max_right_eigenvalue);

const RealType un_star = (0.5 * (uRn + uLn)) + (cL - cR);
const RealType c_star = (0.5 * (cL + cR)) + (.25 * (uLn - uRn));
  
RealType sL = MIN(uLn - cL, un_star - c_star);
RealType sR = MAX(uRn + cR, un_star + c_star);
    
const RealType sL_minus = MIN(sL, zero);
const RealType sR_plus = MAX(sR, zero);
  
const RealType sR_minus_sL = MAX(sR_plus - sL_minus, EPS);

const RealType t1 = sR_plus / sR_minus_sL;
const RealType t2 = one - t1;
// Alternative computation for t2.
// RealType t2 = - sL_minus / sR_minus_sL;
    
ASSERT((0.0 <= t1) && (t1 <= 1.0));
ASSERT((0.0 <= t2) && (t2 <= 1.0));
  
const RealType t3 = sR_plus * sL_minus / sR_minus_sL;

// Physical flux for left facet values.

const RealType left_physical_flux_0 = (rhoL * uLn);
const RealType left_physical_flux_1 = (rhoL * uLn * uL);
const RealType left_physical_flux_2 = (rhoL * uLn * vL);
const RealType left_physical_flux_3 = (rhoL * uLn * wL);

// Physical flux for right facet values.

const RealType right_physical_flux_0 = (rhoL * uLn);
const RealType right_physical_flux_1 = (rhoL * uLn * uL);
const RealType right_physical_flux_2 = (rhoL * uLn * vL);
const RealType right_physical_flux_3 = (rhoL * uLn * wL);

/// Beginning of memory stores -- contiguous.
flux_rho[i] = face_area *  
  (t1 * left_physical_flux_0 +
   t2 * right_physical_flux_0 +
   t3 * (rhoR - rhoL));

flux_rho_u[i] = face_area *
  (t1 * left_physical_flux_1 +
   t2 * right_physical_flux_1 +
   t3 * (uR - uL));
    
flux_rho_v[i] = face_area * 
  (t1 * left_physical_flux_2 +
   t2 * right_physical_flux_2 +
   t3 * (vR - vL));

flux_rho_w[i] = face_area *
  (t1 * left_physical_flux_3 +
   t2 * right_physical_flux_3 +
   t3 * (wR - wL));

face_max_eigen[i] = face_area * max_eigenvalue;

// End of memory stores (contiguous).

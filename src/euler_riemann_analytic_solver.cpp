#include <cassert>
#include <cmath>
#include <iostream>

// Purpose: to provide a guessed value for pressure pm in the star
// region. The choice is made according to adaptive Riemann solver
// using the PVRS, TRRS and TSRS approximate Riemann solvers. See
// Sect. 9.5 of Chapt. 9 of Ref. 1
void GuessStarPressure(RealType gamma,
		       RealType rho_left, RealType ul, RealType pl, 
		       RealType rho_right, RealType ur, RealType pr,
		       RealType *guessed_pressure_ptr) {

  // compute gamma related constants
  const RealType g1 = (gamma - 1.0) / (2.0 * gamma);
  const RealType g2 = (gamma + 1.0) / (2.0 * gamma);
  const RealType g3 = 2.0 * gamma / (gamma - 1.0);
  const RealType g4 = 2.0 / (gamma - 1.0);
  const RealType g5 = 2.0 / (gamma + 1.0);
  const RealType g6 = (gamma - 1.0) / (gamma + 1.0);
  const RealType g7 = (gamma - 1.0) / 2.0;
  const RealType g8 = gamma - 1.0;

  // compute sound speeds
  const RealType cl = sqrt(gamma * pl / rho_left);
  const RealType cr = sqrt(gamma * pr / rho_right);

  const RealType quser = 2.0;

  // compute guess pressure from PVRS Riemann solver
  const RealType cup = 0.25 * (rho_left + rho_right) * (cl + cr);
  RealType ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
  ppv = fmax(0.0, ppv);
  const RealType pmin = fmin(pl, pr);
  const RealType pmax = fmax(pl, pr);
  const RealType qmax = pmax / pmin;

  if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
    *guessed_pressure_ptr = ppv;     // select PVRS Riemann solver

  else {

    if (ppv < pmin) {
      // select Two-Rarefaction Riemann solver
      const RealType pq = pow(pl / pr, g1);
      const RealType um = (pq * ul / cl + ur / cr + g4 * (pq - 1.0)) / (pq / cl + 1.0 / cr);
      const RealType ptl = 1.0 + g7 * (ul - um) / cl;
      const RealType ptr = 1.0 + g7 * (um - ur) / cr;
      *guessed_pressure_ptr = 0.5 * (pow(pl * ptl, g3) + pow(pr * ptr, g3));

    } else {

      // select Two-Shock Riemann solver with PVRS as estimate
      const RealType gel = sqrt((g5 / rho_left)/(g6 * pl + ppv));
      const RealType ger = sqrt((g5 / rho_right)/(g6 * pr + ppv));
      *guessed_pressure_ptr = (gel * pl + ger * pr - (ur - ul)) / (gel + ger);

    }
  }
}

// purpose: to evaluate the pressure functions fl and fr in exact
// Riemann solver and their first derivatives.
static void EvalPressure(RealType gamma,
			 RealType p,
			 RealType dk,
			 RealType pk,
			 RealType ck,
			 RealType *f_ptr,
			 RealType *fd_ptr) {

  // compute gamma related constants
  const RealType g1 = (gamma - 1.0) / (2.0 * gamma);
  const RealType g2 = (gamma + 1.0) / (2.0 * gamma);
  const RealType g3 = 2.0 * gamma / (gamma - 1.0);
  const RealType g4 = 2.0 / (gamma - 1.0);
  const RealType g5 = 2.0 / (gamma + 1.0);
  const RealType g6 = (gamma - 1.0) / (gamma + 1.0);
  const RealType g7 = (gamma - 1.0) / 2.0;
  const RealType g8 = gamma - 1.0;

  if (p <= pk) {
    // rarefaction wave
    const RealType pratio = p / pk;
    *f_ptr = g4 * ck * (pow(pratio, g1) - 1.0);
    *fd_ptr = (1.0 / (dk * ck)) * pow(pratio, - g2);
  } else {
    //  shock wave
    const RealType ak = g5 / dk;
    const RealType bk = g6 * pk;
    const RealType qrt = sqrt(ak / (bk + p));
    *f_ptr = (p - pk) * qrt;
    *fd_ptr = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
  }

}

// Purpose: to compute the solution for pressure and velocity in the
// star region.
void ComputeStarValues(RealType gamma,
		       RealType rho_left, RealType ul, RealType pl, 
		       RealType rho_right, RealType ur, RealType pr, 
		       RealType pscale,
		       RealType* p_ptr, RealType* u_ptr) {

  // compute sound speeds
  RealType cl = sqrt(gamma * pl / rho_left);
  RealType cr = sqrt(gamma * pr / rho_right);

  const int nriter = 20;
  const RealType tolpre = 1.0e-6;
  RealType change, fl, fld, fr, frd, pold, pstart, udiff;

  RealType p = 0.0;
  RealType u = 0.0;

  // guessed value pstart is computed
  GuessStarPressure(gamma, rho_left, ul, pl, rho_right, ur, pr, &pstart);

  pold = pstart;
  udiff = ur - ul;

  // std::cout << "----------------------------------------\n"
  // 	    << "   Iteration number     Change\n"
  // 	    << "----------------------------------------" << std::endl;

  for (int i = 1; i <= nriter; ++i) {

    EvalPressure(gamma, pold, rho_left, pl, cl, &fl, &fld);
    EvalPressure(gamma, pold, rho_right, pr, cr, &fr, &frd);

    p = pold - (fl + fr + udiff) / (fld + frd);
    change = 2.0 * fabs((p - pold) / (p + pold));
    //std::cout << '\t' << i <<  "\t\t" << change << std::endl;
    if (change <= tolpre)
      break;
    if (p < 0.0)
      p = tolpre;

    pold = p;
  }
  // if (i > nriter)
  //   std::cout << "divergence in Newton-Raphson iteration" << std::endl;

  // compute velocity in star region
  u = 0.5*(ul + ur + fr - fl);

  // std::cout << "----------------------------------------\n"
  // 	    << "     Pressure           Velocity\n"
  // 	    << "----------------------------------------\n"
  // 	    << "     " << p / pscale << "\t\t" <<  u << '\n'
  // 	    << "----------------------------------------" << std::endl;

  *p_ptr = p;
  *u_ptr = u;
}

void RiemannAnalyticalSolver(RealType rho_left,
			     RealType u_left,
			     RealType p_left,
			     RealType rho_right,
			     RealType u_right,
			     RealType p_right,
			     RealType gamma,
			     RealType t,
			     RealType x0,
			     int n,
			     const RealType* x,
			     RealType* out_rho,
			     RealType* out_u,
			     RealType* out_p) {


  // compute gamma related constants
  const RealType g1 = (gamma - 1.0) / (2.0 * gamma);
  const RealType g2 = (gamma + 1.0) / (2.0 * gamma);
  const RealType g3 = 2.0 * gamma / (gamma - 1.0);
  const RealType g4 = 2.0 / (gamma - 1.0);
  const RealType g5 = 2.0 / (gamma + 1.0);
  const RealType g6 = (gamma - 1.0) / (gamma + 1.0);
  const RealType g7 = (gamma - 1.0) / 2.0;
  const RealType g8 = gamma - 1.0;

  RealType p_star = 0.0;
  RealType u_star = 0.0;
  
  const RealType pscale = 1.0;

  ComputeStarValues(gamma,
		    rho_left, u_left, p_left, 
		    rho_right, u_right, p_right, 
		    pscale, &p_star, &u_star);
  
  // compute sound speeds
  const RealType c_left = sqrt(gamma * p_left / rho_left);
  const RealType c_right = sqrt(gamma * p_right / rho_right);

  for (int i = 0; i < n; ++i) {

    const RealType s = (x[i] - x0) / t;

    if (s <= u_star) {

      // sampling point lies to the left of the contact discontinuity
      if (p_star <= p_left) {
	// left rarefaction
	const RealType shl = u_left - c_left;
	if (s <= shl) {
	  // sampled point is left data state
	  out_rho[i] = rho_left;
	  out_u[i] = u_left;
	  out_p[i] = p_left;
	} else {
	  const RealType cml = c_left * pow(p_star / p_left, g1);
	  const RealType stl = u_star - cml;
	  if (s > stl) {
	    // samp_lefted point is star left state
	    out_rho[i] = rho_left * pow(p_star / p_left, 1.0 / gamma);
	    out_u[i] = u_star;
	    out_p[i] = p_star;
	  } else {
	    // samp_lefted point is inside left fan
	    out_u[i] = g5 * (c_left + g7 * u_left + s);
	    const RealType c = g5 * (c_left + g7 * (u_left - s));
	    out_rho[i] = rho_left * pow(c / c_left, g4);
	    out_p[i] =  p_left * pow(c / c_left, g3);
	  }
	}
      } else {
	// left shock
	const RealType p_starl = p_star / p_left;
	const RealType sl = u_left - c_left * sqrt(g2 * p_starl + g1);
	if (s <= sl) {
	  // samp_lefted point is left data state
	  out_rho[i] = rho_left;
	  out_u[i] = u_left;
	  out_p[i] = p_left;
	} else {
	  // samp_lefted point is star left state
	  out_rho[i] = rho_left * (p_starl + g6) / (p_starl * g6 + 1.0);
	  out_u[i] = u_star;
	  out_p[i] = p_star;
	}
      }
    } else {
      // sampling point lies to the right of the contact discontinuity
      if (p_star > p_right) {

	// right shock
	const RealType p_starr = p_star / p_right;
	const RealType sr  = u_right + c_right * sqrt(g2 * p_starr + g1);
	if (s >= sr) {
	  // sampled point is right data state
	  out_rho[i] = rho_right;
	  out_u[i] = u_right;
	  out_p[i] = p_right;
	} else {
	  // sampled point is star right state
	  out_rho[i] = rho_right * (p_starr + g6) / (p_starr * g6 + 1.0);
	  out_u[i] = u_star;
	  out_p[i] = p_star;
	}
      } else {
	// right rarefaction
	const RealType shr = u_right + c_right;
	if (s >= shr) {
	  assert(0);
	  // samp_lefted point is right data state
	  out_rho[i] = rho_right;
	  out_u[i] = u_right;
	  out_p[i] = p_right;
	} else {
	  const RealType cmr = c_right * pow(p_star/p_right, g1);
	  const RealType str = u_star + cmr;
	  if (s <= str) {
	    // sampled point is star right state
	    out_rho[i] = rho_right * pow(p_star / p_right, 1.0 / gamma);
	    out_u[i] = u_star;
	    out_p[i] =  p_star;
	  } else {
	    // sampled point is inside left fan
	    out_u[i] = g5 * (- c_right + g7 * u_right + s);
	    const RealType c = g5 * (c_right - g7 * (u_right - s));
	    out_rho[i] = rho_right * pow(c / c_right, g4);
	    out_p[i] =  p_right * pow(c / c_right, g3);
	  }
	}
      }
    }
  }
}

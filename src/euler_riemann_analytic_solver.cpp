#include <cassert>
#include <cmath>

void RiemannAnalyticalSolver(RealType rho_left,
			     RealType u_left,
			     RealType p_left,
			     RealType rho_right,
			     RealType u_right,
			     RealType p_right,
			     RealType p_star, 
			     RealType u_star,
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
  
  // compute sound speeds
  const RealType c_left = sqrt(gamma * p_left / rho_left);
  const RealType c_right = sqrt(gamma * p_right / rho_right);

  for (int i = 0; i < n; ++i) {

    const RealType s = (x[i] - x0) / t;

    // Speed s is x / t.

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
	    out_rho[i] = 0.0;//rho_left * pow(p_star / p_left, 1.0 / gamma);
	    out_u[i] = u_star;
	    out_p[i] = p_star;
	  } else {
	    // samp_lefted point is inside left fan
	    out_u[i] = g5 * (c_left + g7 * u_left + s);
	    const RealType c = g5 * (c_left + g7 * (u_left - s));
	    out_rho[i] = 0.0;//rho_left * pow(c / c_left, g4);
	    out_p[i] =  p_left * pow(c / c_left, g3);
	  }
	}
      } else {
	// left shock
	const RealType p_starl = p_star / p_left;
	const RealType sl = u_left - c_left * sqrt(g2 * p_starl + g1);
	if (s <= sl) {
	  // samp_lefted point is left data state
	  out_rho[i] = 0.0;//rho_left;
	  out_u[i] = u_left;
	  out_p[i] = p_left;
	} else {
	  // samp_lefted point is star left state
	  out_rho[i] = 0.0;//rho_left * (p_starl + g6) / (p_starl * g6 + 1.0);
	  out_u[i] = u_star;
	  out_p[i] = p_star;
	}
      }
    } else {
      // samp_lefting point lies to the right of the contact discontinuity
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
	  out_rho[i] = 0.0;//rho_right * (p_starr + g6) / (p_starr * g6 + 1.0);
	  out_u[i] = u_star;
	  out_p[i] = p_star;
	}
      } else {
	// right rarefaction
	const RealType shr = u_right + c_right;
	if (s >= shr) {
	  assert(0);
	  // samp_lefted point is right data state
	  out_rho[i] = 0.0;//rho_right;
	  out_u[i] = u_right;
	  out_p[i] = p_right;
	} else {
	  const RealType cmr = c_right * pow(p_star/p_right, g1);
	  const RealType str = u_star + cmr;
	  if (s <= str) {
	    // samp_lefted point is star right state
	    out_rho[i] = 0.0;//rho_right * pow(p_star / p_right, 1.0 / gamma);
	    out_u[i] = u_star;
	    out_p[i] =  p_star;
	  } else {
	    // samp_lefted point is inside left fan
	    out_u[i] = g5 * (- c_right + g7 * u_right + s);
	    const RealType c = g5 * (c_right - g7 * (u_right - s));
	    //out_rho[i] = rho_right * pow(c / c_right, g4);
	    out_rho[i] = 0.0;
	    out_p[i] =  p_right * pow(c / c_right, g3);
	  }
	}
      }
    }
  }
}

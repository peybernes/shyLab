#ifndef EULER_RIEMANN_ANALYTIC_SOLVER_HPP
#define EULER_RIEMANN_ANALYTIC_SOLVER_HPP

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
			     RealType* out_p);


#endif // EULER_RIEMANN_ANALYTIC_SOLVER_HPP

// -*- C++-mode -*-
// Fix this ugliness :)

// Reconstructed values at faces (lefft and right).

/* // RHO  */
/* const RealType s_rho_mo = VanAlbadaLimiter(in_rho_m1o - in_rho_m2o, in_rho_ooo - in_rho_m1o); */
/* assert((0.0 <= s_rho_mo) && (s_rho_mo <= 1.0)); */

/* const RealType in_rho_lefft_mo = ReconstructValueLefftKappa(in_rho_m1o, s_rho_mo, kappa, dx, in_rho_m1o - in_rho_m2o, in_rho_ooo - in_rho_m1o); */
/* const RealType in_rho_right_mo = ReconstructValueRightKappa(in_rho_ooo, s_rho_mo, kappa, dx, in_rho_ooo - in_rho_m1o, in_rho_p1o - in_rho_ooo); */

/* const RealType s_rho_po = VanAlbadaLimiter(in_rho_ooo - in_rho_m1o, in_rho_p1o - in_rho_ooo); */
/* //printf("s_rho_po=%lf\n", s_rho_po); */
/* assert((0.0 <= s_rho_po) && (s_rho_po <= 1.0)); */

/* const RealType in_rho_lefft_po = ReconstructValueLefftKappa(in_rho_ooo, s_rho_po, kappa, dx, in_rho_ooo - in_rho_m1o, in_rho_p1o - in_rho_ooo); */
/* const RealType in_rho_right_po = ReconstructValueRightKappa(in_rho_p1o, s_rho_po, kappa, dx, in_rho_p1o - in_rho_ooo, in_rho_p2o - in_rho_p1o); */

/* // VELOCITY_X  */

/* // x direction. */
/* const RealType s_velocity_x_mo = VanAlbadaLimiter(in_velocity_x_m1o - in_velocity_x_m2o, in_velocity_x_ooo - in_velocity_x_m1o); */
/* assert((0.0 <= s_velocity_x_mo) && (s_velocity_x_mo <= 1.0)); */

/* const RealType in_velocity_x_lefft_mo = ReconstructValueLefftKappa(in_velocity_x_m1o, s_velocity_x_mo, kappa, dx, in_velocity_x_m1o - in_velocity_x_m2o, in_velocity_x_ooo - in_velocity_x_m1o); */
/* const RealType in_velocity_x_right_mo = ReconstructValueRightKappa(in_velocity_x_ooo, s_velocity_x_mo, kappa, dx, in_velocity_x_ooo - in_velocity_x_m1o, in_velocity_x_p1o - in_velocity_x_ooo); */

/* const RealType s_velocity_x_po = VanAlbadaLimiter(in_velocity_x_ooo - in_velocity_x_m1o, in_velocity_x_p1o - in_velocity_x_ooo); */
/* assert((0.0 <= s_velocity_x_po) && (s_velocity_x_po <= 1.0)); */

/* const RealType in_velocity_x_lefft_po = ReconstructValueLefftKappa(in_velocity_x_ooo, s_velocity_x_po, kappa, dx, in_velocity_x_ooo - in_velocity_x_m1o, in_velocity_x_p1o - in_velocity_x_ooo); */
/* const RealType in_velocity_x_right_po = ReconstructValueRightKappa(in_velocity_x_p1o, s_velocity_x_po, kappa, dx, in_velocity_x_p1o - in_velocity_x_ooo, in_velocity_x_p2o - in_velocity_x_p1o); */

/* // VELOCITY_Y  */

/* // x direction. */
/* const RealType s_velocity_y_mo = VanAlbadaLimiter(in_velocity_y_m1o - in_velocity_y_m2o, in_velocity_y_ooo - in_velocity_y_m1o); */
/* assert((0.0 <= s_velocity_y_mo) && (s_velocity_y_mo <= 1.0)); */

/* const RealType in_velocity_y_lefft_mo = ReconstructValueLefftKappa(in_velocity_y_m1o, s_velocity_y_mo, kappa, dx, in_velocity_y_m1o - in_velocity_y_m2o, in_velocity_y_ooo - in_velocity_y_m1o); */
/* const RealType in_velocity_y_right_mo = ReconstructValueRightKappa(in_velocity_y_ooo, s_velocity_y_mo, kappa, dx, in_velocity_y_ooo - in_velocity_y_m1o, in_velocity_y_p1o - in_velocity_y_ooo); */

/* const RealType s_velocity_y_po = VanAlbadaLimiter(in_velocity_y_ooo - in_velocity_y_m1o, in_velocity_y_p1o - in_velocity_y_ooo); */
/* assert((0.0 <= s_velocity_y_po) && (s_velocity_y_po <= 1.0)); */

/* const RealType in_velocity_y_lefft_po = ReconstructValueLefftKappa(in_velocity_y_ooo, s_velocity_y_po, kappa, dx, in_velocity_y_ooo - in_velocity_y_m1o, in_velocity_y_p1o - in_velocity_y_ooo); */
/* const RealType in_velocity_y_right_po = ReconstructValueRightKappa(in_velocity_y_p1o, s_velocity_y_po, kappa, dx, in_velocity_y_p1o - in_velocity_y_ooo, in_velocity_y_p2o - in_velocity_y_p1o); */

/* // TOTAL_ENERGY  */

/* // x direction. */
/* const RealType s_total_energy_mo = VanAlbadaLimiter(in_total_energy_m1o - in_total_energy_m2o, in_total_energy_ooo - in_total_energy_m1o); */
/* const RealType s_total_energy_po = VanAlbadaLimiter(in_total_energy_ooo - in_total_energy_m1o, in_total_energy_p1o - in_total_energy_ooo); */

/* const RealType in_total_energy_lefft_mo = ReconstructValueLefftKappa(in_total_energy_m1o, s_total_energy_mo, kappa, dx, in_total_energy_m1o - in_total_energy_m2o, in_total_energy_ooo - in_total_energy_m1o); */
/* const RealType in_total_energy_right_mo = ReconstructValueRightKappa(in_total_energy_ooo, s_total_energy_mo, kappa, dx, in_total_energy_ooo - in_total_energy_m1o, in_total_energy_p1o - in_total_energy_ooo); */

/* const RealType in_total_energy_lefft_po = ReconstructValueLefftKappa(in_total_energy_ooo, s_total_energy_po, kappa, dx, in_total_energy_ooo - in_total_energy_m1o, in_total_energy_p1o - in_total_energy_ooo); */
/* const RealType in_total_energy_right_po = ReconstructValueRightKappa(in_total_energy_p1o, s_total_energy_po, kappa, dx, in_total_energy_p1o - in_total_energy_ooo, in_total_energy_p2o - in_total_energy_p1o); */

/* // Pressure. Determined by an equation of state. Perfect gases. */
/* const RealType gamma = 1.4; */

/* //printf("in_rho_left=%lf, in_rho_right=%lf\n"); */

/* const RealType pressure_lefft_mo = //0.0; */
/*   EquationOfStatePerfectGas(gamma, in_rho_lefft_mo, in_velocity_x_lefft_mo, in_velocity_y_lefft_mo, in_total_energy_lefft_mo); */

/* const RealType pressure_right_mo = //0.0; */
/*   EquationOfStatePerfectGas(gamma, in_rho_right_mo, in_velocity_x_right_mo, in_velocity_y_right_mo, in_total_energy_right_mo); */

/* const RealType pressure_lefft_po = //0.0; */
/*   EquationOfStatePerfectGas(gamma, in_rho_lefft_po, in_velocity_x_lefft_po, in_velocity_y_lefft_po, in_total_energy_lefft_po); */

/* const RealType pressure_right_po = //0.0; */
/*   EquationOfStatePerfectGas(gamma, in_rho_right_po, in_velocity_x_right_po, in_velocity_y_right_po, in_total_energy_right_po); */

/* //printf("pressure_lefft_mo=%lf, pressure_right_mo=%lf, pressure_lefft_po=%lf, pressure_right_po=%lf\n", pressure_lefft_mo, pressure_right_mo, pressure_lefft_po, pressure_right_po); */

/* // RHO */
/* const RealType mean_flux_face_mo_rho =  */
/*   half * ((in_rho_lefft_mo * in_velocity_x_lefft_mo) + (in_rho_right_mo * in_velocity_x_right_mo)); */

/* const RealType mean_flux_face_po_rho =  */
/*   half * ((in_rho_lefft_po * in_velocity_x_lefft_po) + (in_rho_right_po * in_velocity_x_right_po)); */

/* // RHO * VELOCITY_X */
/* const RealType mean_flux_face_mo_momentum_x =  */
/*   half * (((in_rho_lefft_mo * in_velocity_x_lefft_mo * in_velocity_x_lefft_mo) + pressure_lefft_mo) +  */
/* 	  ((in_rho_right_mo * in_velocity_x_right_mo * in_velocity_x_right_mo) + pressure_right_mo)); */

/* const RealType mean_flux_face_po_momentum_x =  */
/*   half * (((in_rho_lefft_po * in_velocity_x_lefft_po * in_velocity_x_lefft_po) + pressure_lefft_po) +  */
/* 	  ((in_rho_right_po * in_velocity_x_right_po * in_velocity_x_right_po) + pressure_right_po)); */

/* //RHO_VELOCITY_Y */
/* const RealType mean_flux_face_mo_momentum_y =  */
/*   half * ((in_rho_lefft_mo * in_velocity_x_lefft_mo * in_velocity_y_lefft_mo) +  */
/* 	  ((in_rho_right_mo * in_velocity_x_right_mo * in_velocity_y_right_mo))); */

/* const RealType mean_flux_face_po_momentum_y =  */
/*   half * ((in_rho_lefft_po * in_velocity_x_lefft_po * in_velocity_y_lefft_po) +  */
/* 	  ((in_rho_right_po * in_velocity_x_right_po * in_velocity_y_right_po))); */

/* // RHO * TOTAL_ENERGY */
/* const RealType mean_flux_face_mo_rho_total_energy = */
/*   half * ((((in_rho_lefft_mo * in_total_energy_lefft_mo) + pressure_lefft_mo) * in_velocity_x_lefft_mo) +  */
/* 	  (((in_rho_right_mo * in_total_energy_right_mo) + pressure_right_mo) * in_velocity_x_right_mo)); */
      
/* const RealType mean_flux_face_po_rho_total_energy = */
/*   half * ((((in_rho_lefft_po * in_total_energy_lefft_po) + pressure_lefft_po) * in_velocity_x_lefft_po) +  */
/* 	  (((in_rho_right_po * in_total_energy_right_po) + pressure_right_po) * in_velocity_x_right_po)); */

/* // We now need average quantities for diffusing the numerical fluxes. */

/* const RealType sqrt_in_rho_lefft_mo = std::sqrt(in_rho_lefft_mo); */
/* const RealType sqrt_in_rho_right_mo = std::sqrt(in_rho_right_mo); */

/* const RealType sqrt_in_rho_lefft_po = std::sqrt(in_rho_lefft_po); */
/* const RealType sqrt_in_rho_right_po = std::sqrt(in_rho_right_po); */

/* const RealType mean_rho_mo = sqrt_in_rho_lefft_mo * sqrt_in_rho_right_mo; */
/* const RealType mean_rho_po = sqrt_in_rho_lefft_po * sqrt_in_rho_right_po; */

/* const RealType mean_velocity_mo = */
/*   ((sqrt_in_rho_lefft_mo * in_velocity_x_lefft_mo) + (sqrt_in_rho_right_mo * in_velocity_x_right_mo)) / (sqrt_in_rho_lefft_mo + sqrt_in_rho_right_mo); */

/* const RealType mean_velocity_po = */
/*   ((sqrt_in_rho_lefft_po * in_velocity_x_lefft_po) + (sqrt_in_rho_right_po * in_velocity_x_right_po)) / (sqrt_in_rho_lefft_po + sqrt_in_rho_right_po); */

/* const RealType mean_pressure_mo = */
/*   ((sqrt_in_rho_lefft_mo * pressure_lefft_mo) + (sqrt_in_rho_right_mo * pressure_right_mo)) / (sqrt_in_rho_lefft_mo + sqrt_in_rho_right_mo); */

/* const RealType mean_pressure_po = */
/*   ((sqrt_in_rho_lefft_po * pressure_lefft_po) + (sqrt_in_rho_right_po * pressure_right_po)) / (sqrt_in_rho_lefft_po + sqrt_in_rho_right_po); */

/* // Numerical fluxes at faces (Lax Friedrichs). */
/* const RealType diffusion_coefficient_mo = std::abs(mean_velocity_mo) + std::sqrt((gamma * mean_pressure_mo) / mean_rho_mo); */
/* const RealType diffusion_coefficient_po = std::abs(mean_velocity_po) + std::sqrt((gamma * mean_pressure_po) / mean_rho_po); */
								  
/* const RealType flux_face_mo_rho = cx * (mean_flux_face_mo_rho - half * diffusion_coefficient_mo * (in_rho_right_mo - in_rho_lefft_mo)); */
/* const RealType flux_face_po_rho = cx * (mean_flux_face_po_rho - half * diffusion_coefficient_po * (in_rho_right_po - in_rho_lefft_po)); */

/* const RealType flux_face_mo_momentum_x =  */
/*   cx * (mean_flux_face_mo_momentum_x -  */
/* 	half * diffusion_coefficient_mo * (in_rho_right_mo * in_velocity_x_right_mo - in_rho_lefft_mo * in_velocity_x_lefft_mo)); */

/* const RealType flux_face_po_momentum_x =  */
/*   cx * (mean_flux_face_po_momentum_x -  */
/* 	half * diffusion_coefficient_po * (in_rho_right_po * in_velocity_x_right_po - in_rho_lefft_po * in_velocity_x_lefft_po)); */

/* //printf("flux_face_momentum_mo=%lf, flux_face_momentum_po=%lf\n", flux_face_mo_momentum_x, flux_face_po_momentum_x); */

/* const RealType flux_face_mo_momentum_y =  */
/*   cx * (mean_flux_face_mo_momentum_y -  */
/* 	half * diffusion_coefficient_mo * (in_rho_right_mo * in_velocity_y_right_mo - in_rho_lefft_mo * in_velocity_y_lefft_mo)); */

/* const RealType flux_face_po_momentum_y =  */
/*   cx * (mean_flux_face_po_momentum_y -  */
/* 	half * diffusion_coefficient_po * (in_rho_right_po * in_velocity_y_right_po - in_rho_lefft_po * in_velocity_y_lefft_po)); */

/* const RealType flux_face_mo_total_energy =  */
/*   cx * (mean_flux_face_mo_rho_total_energy -  */
/* 	half * diffusion_coefficient_mo * (in_rho_right_mo * in_total_energy_right_mo - in_rho_lefft_mo * in_total_energy_lefft_mo)); */

/* const RealType flux_face_po_total_energy =  */
/*   cx * (mean_flux_face_po_rho_total_energy -  */
/* 	half * diffusion_coefficient_po * (in_rho_right_po * in_total_energy_right_po - in_rho_lefft_po * in_total_energy_lefft_po)); */

/* const RealType out_rho_ooo =  */
/*   dt * (flux_face_mo_rho - flux_face_po_rho); */


/* const RealType out_velocity_x_ooo = */
/*   dt * (flux_face_mo_momentum_x - flux_face_po_momentum_x); */

/* const RealType out_velocity_y_ooo =  */
/*   dt * (flux_face_mo_momentum_y - flux_face_po_momentum_y); */

/* const RealType out_total_energy_ooo =  */
/*   dt * (flux_face_mo_total_energy - flux_face_po_total_energy); */

const RealType out_rho_ooo = 0.0;


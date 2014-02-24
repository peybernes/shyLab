const RealType moment_flux_prev =
  (0.5 * (prev_dual_mass_flux + fabs(prev_dual_mass_flux)) * in_vx_m1) + 
  (0.5 * (prev_dual_mass_flux - fabs(prev_dual_mass_flux)) * in_vx_oo);

const RealType moment_flux_next =
  (0.5 * (next_dual_mass_flux + fabs(next_dual_mass_flux)) * in_vx_oo) + 
  (0.5 * (next_dual_mass_flux - fabs(next_dual_mass_flux)) * in_vx_p1);

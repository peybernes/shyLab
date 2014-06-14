const RealType moment_x_flux_prev_x =
  (0.5 * (prev_dual_mass_flux_x + fabs(prev_dual_mass_flux_x)) * in_vx_m1o) + 
  (0.5 * (prev_dual_mass_flux_x - fabs(prev_dual_mass_flux_x)) * in_vx_ooo);

const RealType moment_y_flux_prev_x =
  (0.5 * (prev_dual_mass_flux_x + fabs(prev_dual_mass_flux_x)) * in_vy_m1o) + 
  (0.5 * (prev_dual_mass_flux_x - fabs(prev_dual_mass_flux_x)) * in_vy_ooo);

const RealType moment_x_flux_prev_y =
  (0.5 * (prev_dual_mass_flux_y + fabs(prev_dual_mass_flux_y)) * in_vx_om1) + 
  (0.5 * (prev_dual_mass_flux_y - fabs(prev_dual_mass_flux_y)) * in_vx_ooo);

const RealType moment_y_flux_prev_y =
  (0.5 * (prev_dual_mass_flux_y + fabs(prev_dual_mass_flux_y)) * in_vy_om1) + 
  (0.5 * (prev_dual_mass_flux_y - fabs(prev_dual_mass_flux_y)) * in_vy_ooo);

const RealType moment_x_flux_next_x =
  (0.5 * (next_dual_mass_flux_x + fabs(next_dual_mass_flux_x)) * in_vx_ooo) + 
  (0.5 * (next_dual_mass_flux_x - fabs(next_dual_mass_flux_x)) * in_vx_p1o);

const RealType moment_y_flux_next_x =
  (0.5 * (next_dual_mass_flux_x + fabs(next_dual_mass_flux_x)) * in_vy_ooo) + 
  (0.5 * (next_dual_mass_flux_x - fabs(next_dual_mass_flux_x)) * in_vy_p1o);

const RealType moment_x_flux_next_y =
  (0.5 * (next_dual_mass_flux_y + fabs(next_dual_mass_flux_y)) * in_vx_ooo) + 
  (0.5 * (next_dual_mass_flux_y - fabs(next_dual_mass_flux_y)) * in_vx_op1);

const RealType moment_y_flux_next_y =
  (0.5 * (next_dual_mass_flux_y + fabs(next_dual_mass_flux_y)) * in_vy_ooo) + 
  (0.5 * (next_dual_mass_flux_y - fabs(next_dual_mass_flux_y)) * in_vy_op1);

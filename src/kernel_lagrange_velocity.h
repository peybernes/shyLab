// Kernel for Lagrangian velocity calculation 
// to be called in kernel_lagrange_2d.c  

  const RealType m_sw = in_mass[cell_SW];
  const RealType p_sw = in_pressure[cell_SW];
  const RealType q_sw = in_pseudo_pressure[cell_SW];
  const RealType m_se = in_mass[cell_SE];
  const RealType p_se = in_pressure[cell_SE];
  const RealType q_se = in_pseudo_pressure[cell_SE];
  const RealType m_nw = in_mass[cell_NW];
  const RealType p_nw = in_pressure[cell_NW];
  const RealType q_nw = in_pseudo_pressure[cell_NW];
  const RealType m_ne = in_mass[cell_NE];
  const RealType p_ne = in_pressure[cell_NE];
  const RealType q_ne = in_pseudo_pressure[cell_NE];

  const RealType one_over_mass_ooo = 4.0/(m_sw + m_se + m_nw + m_ne );

  RealType delta_p_x = 0.5 * ((p_se + q_se + p_ne + q_ne   )
				   -(p_sw + q_sw + p_nw + q_nw));
  RealType delta_p_y = 0.5 * ((p_nw + q_nw + p_ne + q_ne   )
				   -(p_sw + q_sw + p_se + q_se));

  RealType out_u_x  = in_velocity_x[node_ooo] -0.5*dt*delta_p_x *one_over_mass_ooo*dy; 
  RealType out_u_y  = in_velocity_y[node_ooo] -0.5*dt*delta_p_y *one_over_mass_ooo*dx;

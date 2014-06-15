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

SHY_ASM_COMMENT("LagrangeVelocityPredicted -- COMPUTE BEGIN");

const RealType one = 1.0;
const RealType half = 0.5;
const RealType one_quarter = 0.25;

const RealType nodal_mass_ooo = one_quarter * (m_sw + m_se + m_nw + m_ne); // 3 ADD, 1 MUL
const RealType one_over_nodal_mass_ooo = one / nodal_mass_ooo; // 1 DIV

const RealType effective_p_nw = p_nw + q_nw; // 1 ADD
const RealType effective_p_ne = p_ne + q_ne; // 1 ADD

const RealType effective_p_sw = p_sw + q_sw; // 1 ADD
const RealType effective_p_se = p_se + q_se; // 1 ADD

const RealType delta_px = half * 
  ((effective_p_se + effective_p_ne) - (effective_p_sw + effective_p_nw)); // 3 ADD, 1 MUL

const RealType delta_py = half * 
  ((effective_p_nw + effective_p_ne) - (effective_p_sw + effective_p_se)); // 3 ADD, 1 MUL

const RealType scaling_factor = half * dt * one_over_nodal_mass_ooo; // 2 MUL

const RealType out_u_x  = 
  in_velocity_x[node_ooo] - scaling_factor * delta_px * dy; // 2 MUL, 1 ADD

const RealType out_u_y  = 
  in_velocity_y[node_ooo] - scaling_factor * delta_py * dx; // 2 MUL, 1 ADD

// TOTAL : 15 ADD, 9 MUL, 1 DIV. approx 31 FLOP

SHY_ASM_COMMENT("LagrangeVelocityPredicted -- COMPUTE END");

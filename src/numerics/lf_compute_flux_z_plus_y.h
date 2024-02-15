for(int k = 0; k < nb_mat; k++) {
  
  masse_fluxes_k_y[k][face_top] = v_et[face_top] * (in_rho_k[k][cell_op1]
				    - dy * 0.5 * rhok_grady[k][cell_op1])
                                    * (in_c_k[k][cell_op1]
                                    - dy * 0.5 * alphak_grady_bot[k][cell_op1]);//*/
  
  alpha_beta_fluxes_k_y[k][face_top] = v_et[face_top]
                                        * (in_c_k[k][cell_op1] - dy * 0.5 * alphak_grady_bot[k][cell_op1])
                                        * (beta[cell_op1] - 0.5 * dy * beta_grady[cell_op1]);
  
}

rho_total_energy_fluxes_y[face_top]  = v_et[face_top] * ((in_rho[cell_op1] - dy * 0.5 * rho_grady_bot[cell_op1])
				        * 0.5 * (u2[cell_op1] - 0.5 * dy * u2_grady[cell_op1]
				        + v2[cell_op1] - 0.5 * dy * v2_grady[cell_op1])
				        + (rho_e[cell_op1] - dy * 0.5 * rho_e_grady_bot[cell_op1]))
                                        + p_yet[face_top] * v_et[face_top];

rho_U_fluxes_y[face_top]  = v_et[face_top] * (in_rho[cell_op1] - dy * 0.5 * rho_grady_bot[cell_op1])
                             * (in_u_cell[cell_op1] - 0.5 * dy * u_grady_bot[cell_op1]);// + p_yet[face_top];

rho_V_fluxes_y[face_top] = v_et[face_top]  * (in_rho[cell_op1] - dy * 0.5 * rho_grady_bot[cell_op1])
                             * (in_v_cell[cell_op1] - 0.5 * dy * v_grady_bot[cell_op1]) + p_yet[face_top];

beta_fluxes_y[face_top]   = v_et[face_top] * (beta[cell_op1] - 0.5 * dy * beta_grady[cell_op1]);


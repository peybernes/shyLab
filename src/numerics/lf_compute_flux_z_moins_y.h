for(int k = 0; k < nb_mat; k++){
  
  masse_fluxes_k_y[k][face_top] = v_et[face_top] * (in_rho_k[k][cell_ooo]
				    + dy * 0.5 * rhok_grady[k][cell_ooo])
                                    * (in_c_k[k][cell_ooo]
                                    + dy * 0.5 * alphak_grady_top[k][cell_ooo]);//*/
  
  alpha_beta_fluxes_k_y[k][face_top] = v_et[face_top]
                                        * (in_c_k[k][cell_ooo] + dy * 0.5 * alphak_grady_top[k][cell_ooo])
                                        * (beta[cell_ooo] + 0.5 * dy * beta_grady[cell_ooo]);
  
}

rho_total_energy_fluxes_y[face_top]  = v_et[face_top] * ((in_rho[cell_ooo] + dy * 0.5 * rho_grady_top[cell_ooo])
				        * 0.5 * (u2[cell_ooo] + 0.5 * dy * u2_grady[cell_ooo]
				        + v2[cell_ooo] + 0.5 * dy * v2_grady[cell_ooo])
				        + (rho_e[cell_ooo] + dy * 0.5 * rho_e_grady_top[cell_ooo]))
                                        + p_yet[face_top] * v_et[face_top];

rho_U_fluxes_y[face_top]  = v_et[face_top] * (in_rho[cell_ooo] + dy * 0.5 * rho_grady_top[cell_ooo])
                             * (in_u_cell[cell_ooo] + 0.5 * dy * u_grady_top[cell_ooo]); // + p_yet[face_top];

rho_V_fluxes_y[face_top] = v_et[face_top] *  (in_rho[cell_ooo] + dy * 0.5 * rho_grady_top[cell_ooo])
                             * (in_v_cell[cell_ooo] + 0.5 * dy * v_grady_top[cell_ooo]) + p_yet[face_top];

beta_fluxes_y[face_top]   = v_et[face_top] * (beta[cell_ooo] + 0.5 * dy * beta_grady[cell_ooo]);

 

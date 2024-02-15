for(int k = 0; k < nb_mat; k++){
  
  masse_fluxes_k_x[k][face_right] = u_et[face_right] * (in_rho_k[k][cell_ooo]
				    + dx * 0.5 * rhok_gradx[k][cell_ooo])
                                    * (in_c_k[k][cell_ooo]
                                    + dx * 0.5 * alphak_gradx_right[k][cell_ooo]);//*/
  
  alpha_beta_fluxes_k_x[k][face_right] = u_et[face_right]
                                        * (in_c_k[k][cell_ooo] + dx * 0.5 * alphak_gradx_right[k][cell_ooo])
                                        * (beta[cell_ooo] + 0.5 * dx * beta_gradx[cell_ooo]);
  
}

rho_total_energy_fluxes_x[face_right]  = u_et[face_right] * ((in_rho[cell_ooo] + dx * 0.5 * rho_gradx_right[cell_ooo])
				        * 0.5 * (u2[cell_ooo] + 0.5 * dx * u2_gradx[cell_ooo]
				        + v2[cell_ooo] + 0.5 * dx * v2_gradx[cell_ooo])
				        + (rho_e[cell_ooo] + dx * 0.5 * rho_e_gradx_right[cell_ooo]))
                                        + p_xet[face_right] * u_et[face_right];

rho_U_fluxes_x[face_right]  = u_et[face_right] * (in_rho[cell_ooo] + dx * 0.5 * rho_gradx_right[cell_ooo])
                             * (in_u_cell[cell_ooo] + 0.5 * dx * u_gradx_right[cell_ooo]) + p_xet[face_right];

rho_V_fluxes_x[face_right]  = u_et[face_right] * (in_rho[cell_ooo] + dx * 0.5 * rho_gradx_right[cell_ooo])
                             * (in_v_cell[cell_ooo] + 0.5 * dx * v_gradx_right[cell_ooo]);//+p_xet[face_right];

beta_fluxes_x[face_right]   = u_et[face_right] * (beta[cell_ooo] + 0.5 * dx * beta_gradx[cell_ooo]);

 

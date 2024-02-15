for(int k = 0; k < nb_mat; k++) {
  
  masse_fluxes_k_x[k][face_right] = u_et[face_right] * (in_rho_k[k][cell_p1o]
				    - dx * 0.5 * rhok_gradx[k][cell_p1o])
                                    * (in_c_k[k][cell_p1o]
                                    - dx * 0.5 * alphak_gradx_left[k][cell_p1o]);//*/
  
  alpha_beta_fluxes_k_x[k][face_right] = u_et[face_right]
                                        * (in_c_k[k][cell_p1o] - dx * 0.5 * alphak_gradx_left[k][cell_p1o])
                                        * (beta[cell_p1o] - 0.5 * dx * beta_gradx[cell_p1o]);
  
}

rho_total_energy_fluxes_x[face_right]  = u_et[face_right] * ((in_rho[cell_p1o] - dx * 0.5 * rho_gradx_left[cell_p1o])
				        * 0.5 * (u2[cell_p1o] - 0.5 * dx * u2_gradx[cell_p1o]
				        + v2[cell_p1o] - 0.5 * dx * v2_gradx[cell_p1o])
				        + (rho_e[cell_p1o] - dx * 0.5 * rho_e_gradx_left[cell_p1o]))
                                        + p_xet[face_right] * u_et[face_right];

rho_U_fluxes_x[face_right]  = u_et[face_right] * (in_rho[cell_p1o] - dx * 0.5 * rho_gradx_left[cell_p1o])
                             * (in_u_cell[cell_p1o] - 0.5 * dx * u_gradx_left[cell_p1o]) + p_xet[face_right];

rho_V_fluxes_x[face_right]  = u_et[face_right] * (in_rho[cell_p1o] - dx * 0.5 * rho_gradx_left[cell_p1o])
                             * (in_v_cell[cell_p1o] - 0.5 * dx * v_gradx_left[cell_p1o]);//+p_xet[face_right];

beta_fluxes_x[face_right]   = u_et[face_right] * (beta[cell_p1o] - 0.5 * dx * beta_gradx[cell_p1o]);


	p_xet[face_right] = (((in_rho[cell_ooo] + 0.5 * dx * rho_gradx_right [cell_ooo]) * (pressure[cell_p1o] - 0.5 * dx * p_gradx_left  [cell_p1o])
			   + ( in_rho[cell_p1o] - 0.5 * dx * rho_gradx_left  [cell_p1o]) * (pressure[cell_ooo] + 0.5 * dx * p_gradx_right [cell_ooo]))-
			     (( in_rho[cell_ooo] + 0.5 * dx * rho_gradx_right [cell_ooo]) *(in_rho  [cell_p1o] - 0.5 * dx * rho_gradx_left[cell_p1o]))
			    * std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_p1o]) *
			    (in_u_cell[cell_p1o] - 0.5 * dx * u_gradx_left   [cell_p1o]
			    -in_u_cell[cell_ooo] - 0.5 * dx * u_gradx_right  [cell_ooo]))
	                    /(in_rho  [cell_ooo] + 0.5 * dx * rho_gradx_right[cell_ooo]  + in_rho[cell_p1o] - 0.5 * dx * rho_gradx_left[cell_p1o]);
	
	u_et[face_right] = ((  in_rho[cell_p1o]-0.5*dx*rho_gradx_left[cell_p1o])
			   * (in_u_cell[cell_p1o] - 0.5 * dx * u_gradx_left   [cell_p1o])
			   + (in_rho   [cell_ooo] + 0.5 * dx * rho_gradx_right[cell_ooo]) * (in_u_cell[cell_ooo] + 0.5 * dx * u_gradx_right[cell_ooo])
			   - (pressure [cell_p1o] - 0.5 * dx * p_gradx_left[cell_p1o] - pressure[cell_ooo] - 0.5 * dx * p_gradx_right[cell_ooo])
			   / std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_p1o]))
	                   /(in_rho[cell_ooo] + 0.5 * dx * rho_gradx_right[cell_ooo] + in_rho[cell_p1o] - 0.5 * dx * rho_gradx_left[cell_p1o]);
	
	p_yet[face_top] = ((( in_rho [cell_ooo]  + 0.5 * dy * rho_grady_top [cell_ooo])
			  * (pressure[cell_op1]  - 0.5 * dy * p_grady_bot   [cell_op1])
			  + (in_rho  [cell_op1]  - 0.5 * dy * rho_grady_bot [cell_op1])
			  * (pressure[cell_ooo]  + 0.5 * dy * p_grady_top   [cell_ooo]))-
			 (( in_rho   [cell_ooo]  + 0.5 * dy * rho_grady_top [cell_ooo])
			  *( in_rho  [cell_op1]  - 0.5 * dy * rho_grady_bot [cell_op1]))
			  *std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_op1]) *
			  (in_v_cell[cell_op1]  - 0.5 * dy * v_grady_bot[cell_op1]
			  -in_v_cell[cell_ooo]  - 0.5 * dy * v_grady_top[cell_ooo]))
         	          /(in_rho[cell_ooo]    + 0.5 * dy * rho_grady_top[cell_ooo] + in_rho[cell_op1] - 0.5 * dy * rho_grady_bot[cell_op1]);
	
        v_et[face_top] = ((in_rho   [cell_op1] - 0.5 * dy * rho_grady_bot  [cell_op1])
			* (in_v_cell[cell_op1] - 0.5 * dy * v_grady_bot    [cell_op1])
			+ (in_rho   [cell_ooo] + 0.5 * dy * rho_grady_top  [cell_ooo]) * (in_v_cell[cell_ooo] + 0.5 * dy * v_grady_top[cell_ooo])
			- (pressure [cell_op1] - 0.5 * dy * p_grady_bot    [cell_op1]  -   pressure[cell_ooo] - 0.5 * dy * p_grady_top[cell_ooo])
			/ std::max(speed_of_sound_mix[cell_ooo],speed_of_sound_mix[cell_op1]))
	                /( in_rho   [cell_ooo] + 0.5 * dy * rho_grady_top[cell_ooo] + in_rho[cell_op1] - 0.5 * dy * rho_grady_bot[cell_op1]);


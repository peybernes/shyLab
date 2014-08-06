const RealType one_quarter = 0.25;
const RealType one_half = 0.5;

       //data load
        const RealType in_variable_ooo = in_variable[node_ooo];
        const RealType in_variable_m1o = in_variable[node_m1o];
        const RealType in_variable_p1o = in_variable[node_p1o];
        const RealType in_variable_om1 = in_variable[node_om1];
        const RealType in_variable_op1 = in_variable[node_op1];
	
	const RealType u_m1 = u_velocity_pred[node_m1o];
	const RealType u_oo = u_velocity_pred[node_ooo];
	const RealType u_p1 = u_velocity_pred[node_p1o];
	const RealType v_m1 = v_velocity_pred[node_om1];
	const RealType v_oo = v_velocity_pred[node_ooo];
	const RealType v_p1 = v_velocity_pred[node_op1];
	
	const RealType gradient_variable_x_oo = gradient_variable_x[node_ooo];
        const RealType gradient_variable_x_m1 = gradient_variable_x[node_m1o];
        const RealType gradient_variable_x_p1 = gradient_variable_x[node_p1o];
	const RealType gradient_variable_y_oo = gradient_variable_y[node_ooo];
        const RealType gradient_variable_y_m1 = gradient_variable_y[node_om1];
        const RealType gradient_variable_y_p1 = gradient_variable_y[node_op1];

        const RealType cell_mass_m1m1 = lag_cell_mass[cellm1m1];
        const RealType cell_mass_p1m1 = lag_cell_mass[cellp1m1];
        const RealType cell_mass_m1p1 = lag_cell_mass[cellm1p1];
        const RealType cell_mass_p1p1 = lag_cell_mass[cellp1p1];

        const RealType out_mass_m1m1 = out_cell_mass[cellm1m1];
        const RealType out_mass_p1m1 = out_cell_mass[cellp1m1];
        const RealType out_mass_m1p1 = out_cell_mass[cellm1p1];
        const RealType out_mass_p1p1 = out_cell_mass[cellp1p1];

        const RealType mass_flux_x_m1m1 = mass_flux_x[facexm1m1];
        const RealType mass_flux_x_oom1 = mass_flux_x[facexoom1];
        const RealType mass_flux_x_m1p1 = mass_flux_x[facexm1p1];
        const RealType mass_flux_x_oop1 = mass_flux_x[facexoop1];
        const RealType mass_flux_x_p1p1 = mass_flux_x[facexp1p1];
        const RealType mass_flux_x_p1m1 = mass_flux_x[facexp1m1];
        const RealType mass_flux_y_m1m1 = mass_flux_y[faceym1m1];
        const RealType mass_flux_y_m1oo = mass_flux_y[faceym1oo];
        const RealType mass_flux_y_m1p1 = mass_flux_y[faceym1p1];
        const RealType mass_flux_y_p1oo = mass_flux_y[faceyp1oo];
        const RealType mass_flux_y_p1p1 = mass_flux_y[faceyp1p1];
        const RealType mass_flux_y_p1m1 = mass_flux_y[faceyp1m1];
        //


        const RealType lag_node_mass_ooo = one_quarter * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
        assert(0.0 < lag_node_mass_ooo);
        const RealType out_node_mass_ooo = one_quarter * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
        assert(0.0 < out_node_mass_ooo);

        const RealType prev_dual_mass_flux_x = one_quarter * (mass_flux_x_m1m1 + mass_flux_x_oom1 + mass_flux_x_m1p1 + mass_flux_x_oop1);
	const RealType next_dual_mass_flux_x = one_quarter * (mass_flux_x_p1m1 + mass_flux_x_oom1 + mass_flux_x_p1p1 + mass_flux_x_oop1); 
        const RealType prev_dual_mass_flux_y = one_quarter * (mass_flux_y_m1m1 + mass_flux_y_m1oo + mass_flux_y_p1m1 + mass_flux_y_p1oo);
	const RealType next_dual_mass_flux_y = one_quarter * (mass_flux_y_m1p1 + mass_flux_y_m1oo + mass_flux_y_p1p1 + mass_flux_y_p1oo);

const RealType dx_lag_m1 =  dx + dt * (u_oo - u_m1);
const RealType dx_lag_p1 =  dx + dt * (u_p1 - u_oo);
const RealType dy_lag_m1 =  dy + dt * (v_oo - v_m1);
const RealType dy_lag_p1 =  dy + dt * (v_p1 - v_oo);

const RealType dt_uf_m1 = one_half * dt * (u_oo + u_m1);
const RealType dt_uf_p1 = one_half * dt * (u_oo + u_p1);
const RealType dt_vf_m1 = one_half * dt * (v_oo + v_m1);
const RealType dt_vf_p1 = one_half * dt * (v_oo + v_p1);

RealType nodal_reconstruct_variable_o2_x_m1;
RealType nodal_reconstruct_variable_o2_x_p1 ;
RealType nodal_reconstruct_variable_o2_y_m1;
RealType nodal_reconstruct_variable_o2_y_p1 ;

if (prev_dual_mass_flux_x > 0 ) {
  nodal_reconstruct_variable_o2_x_m1 = in_variable_m1o + one_half * gradient_variable_x_m1 * (dx_lag_m1 - dt_uf_m1);
 } else {
  nodal_reconstruct_variable_o2_x_m1 = in_variable_ooo + one_half * gradient_variable_x_oo * ( -dx_lag_m1 - dt_uf_m1);
 } 

if (next_dual_mass_flux_x > 0 ) {
  nodal_reconstruct_variable_o2_x_p1 = in_variable_ooo + one_half * gradient_variable_x_oo *(dx_lag_p1 - dt_uf_p1);
 } else {
  nodal_reconstruct_variable_o2_x_p1 = in_variable_p1o + one_half * gradient_variable_x_p1 *( -dx_lag_p1 - dt_uf_p1);
 }

if (prev_dual_mass_flux_y > 0 ) {
  nodal_reconstruct_variable_o2_y_m1 = in_variable_om1 + one_half * gradient_variable_y_m1 * (dy_lag_m1 - dt_vf_m1);
 } else {
  nodal_reconstruct_variable_o2_y_m1 = in_variable_ooo + one_half * gradient_variable_y_oo * ( -dy_lag_m1 - dt_vf_m1);
 } 

if (next_dual_mass_flux_y > 0 ) {
  nodal_reconstruct_variable_o2_y_p1 = in_variable_ooo + one_half * gradient_variable_y_oo *(dy_lag_p1 - dt_vf_p1);
 } else {
  nodal_reconstruct_variable_o2_y_p1 = in_variable_op1 + one_half * gradient_variable_y_p1 * ( -dy_lag_p1 - dt_vf_p1);
 }

const RealType moment_flux_prev_x = nodal_reconstruct_variable_o2_x_m1 * prev_dual_mass_flux_x ;
const RealType moment_flux_next_x = nodal_reconstruct_variable_o2_x_p1 * next_dual_mass_flux_x ;
const RealType moment_flux_prev_y = nodal_reconstruct_variable_o2_y_m1 * prev_dual_mass_flux_y ;
const RealType moment_flux_next_y = nodal_reconstruct_variable_o2_y_p1 * next_dual_mass_flux_y ; 
          
	const RealType in_moment_ooo = lag_node_mass_ooo * in_variable_ooo; 

	const RealType out_moment_ooo = 
	  in_moment_ooo + moment_flux_prev_x - moment_flux_next_x + moment_flux_prev_y - moment_flux_next_y ;

	const RealType out_nodal_variable =  out_moment_ooo / out_node_mass_ooo ;

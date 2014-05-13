  //data load
        const RealType in_variable_oo = in_variable[node_ooo];
        const RealType in_variable_m1 = in_variable[node_om1];
        const RealType in_variable_p1 = in_variable[node_op1];
	
	const RealType v_m1 = v_velocity_pred[node_om1];
	const RealType v_oo = v_velocity_pred[node_ooo];
	const RealType v_p1 = v_velocity_pred[node_op1];
	
	const RealType gradient_variable_oo = gradient_variable[node_ooo];
        const RealType gradient_variable_m1 = gradient_variable[node_om1];
        const RealType gradient_variable_p1 = gradient_variable[node_op1];

        const RealType cell_mass_m1m1 = lag_cell_mass[cellm1m1];
        const RealType cell_mass_p1m1 = lag_cell_mass[cellp1m1];
        const RealType cell_mass_m1p1 = lag_cell_mass[cellm1p1];
        const RealType cell_mass_p1p1 = lag_cell_mass[cellp1p1];

        const RealType out_mass_m1m1 = out_cell_mass[cellm1m1];
        const RealType out_mass_p1m1 = out_cell_mass[cellp1m1];
        const RealType out_mass_m1p1 = out_cell_mass[cellm1p1];
        const RealType out_mass_p1p1 = out_cell_mass[cellp1p1];

        const RealType mass_flux_m1m1 = mass_flux[faceym1m1];
        const RealType mass_flux_m1oo = mass_flux[faceym1oo];
        const RealType mass_flux_m1p1 = mass_flux[faceym1p1];
        const RealType mass_flux_p1oo = mass_flux[faceyp1oo];
        const RealType mass_flux_p1p1 = mass_flux[faceyp1p1];
        const RealType mass_flux_p1m1 = mass_flux[faceyp1m1];
        //


        const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
        assert(0.0 < lag_node_mass_ooo);
        const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
        assert(0.0 < out_node_mass_ooo);

        const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_m1oo + mass_flux_p1m1 + mass_flux_p1oo);
	const RealType next_dual_mass_flux = 0.25 * (mass_flux_m1p1 + mass_flux_m1oo + mass_flux_p1p1 + mass_flux_p1oo); 

const RealType dy_lag_m1 =  dy + dt * (v_oo - v_m1);
const RealType dy_lag_p1 =  dy + dt * (v_p1 - v_oo);

const RealType dt_vf_m1 = 0.5 * dt * (v_oo + v_m1);
const RealType dt_vf_p1 = 0.5 * dt * (v_oo + v_p1);

RealType nodal_reconstruct_variable_o2_m1;
RealType nodal_reconstruct_variable_o2_p1 ;

if (prev_dual_mass_flux > 0 ) {
  nodal_reconstruct_variable_o2_m1 = in_variable_m1 + 0.5 * gradient_variable_m1 * (dy_lag_m1 - dt_vf_m1);
 } else {
  nodal_reconstruct_variable_o2_m1 = in_variable_oo + 0.5 * gradient_variable_oo * ( -dy_lag_m1 - dt_vf_m1);
 } 

if (next_dual_mass_flux > 0 ) {
  nodal_reconstruct_variable_o2_p1 = in_variable_oo + 0.5 * gradient_variable_oo *(dy_lag_p1 - dt_vf_p1);
 } else {
  nodal_reconstruct_variable_o2_p1 = in_variable_p1 + 0.5 * gradient_variable_p1 * ( -dy_lag_p1 - dt_vf_p1);
 }

const RealType moment_flux_prev = nodal_reconstruct_variable_o2_m1 * prev_dual_mass_flux ;
const RealType moment_flux_next = nodal_reconstruct_variable_o2_p1 * next_dual_mass_flux ;
 
          
	const RealType in_moment_ooo = lag_node_mass_ooo * in_variable_oo; 
       
	const RealType out_moment_ooo = 
	  in_moment_ooo + moment_flux_prev - moment_flux_next;

	const RealType out_nodal_variable =  out_moment_ooo / out_node_mass_ooo ;

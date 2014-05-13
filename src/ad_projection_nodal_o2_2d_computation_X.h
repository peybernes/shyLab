       //data load
        const RealType in_variable_oo = in_variable[node_ooo];
        const RealType in_variable_m1 = in_variable[node_m1o];
        const RealType in_variable_p1 = in_variable[node_p1o];
	
	const RealType u_m1 = u_velocity_pred[node_m1o];
	const RealType u_oo = u_velocity_pred[node_ooo];
	const RealType u_p1 = u_velocity_pred[node_p1o];
	
	const RealType gradient_variable_oo = gradient_variable[node_ooo];
        const RealType gradient_variable_m1 = gradient_variable[node_m1o];
        const RealType gradient_variable_p1 = gradient_variable[node_p1o];

        const RealType cell_mass_m1m1 = lag_cell_mass[cellm1m1];
        const RealType cell_mass_p1m1 = lag_cell_mass[cellp1m1];
        const RealType cell_mass_m1p1 = lag_cell_mass[cellm1p1];
        const RealType cell_mass_p1p1 = lag_cell_mass[cellp1p1];

        const RealType out_mass_m1m1 = out_cell_mass[cellm1m1];
        const RealType out_mass_p1m1 = out_cell_mass[cellp1m1];
        const RealType out_mass_m1p1 = out_cell_mass[cellm1p1];
        const RealType out_mass_p1p1 = out_cell_mass[cellp1p1];

        const RealType mass_flux_m1m1 = mass_flux[facexm1m1];
        const RealType mass_flux_oom1 = mass_flux[facexoom1];
        const RealType mass_flux_m1p1 = mass_flux[facexm1p1];
        const RealType mass_flux_oop1 = mass_flux[facexoop1];
        const RealType mass_flux_p1p1 = mass_flux[facexp1p1];
        const RealType mass_flux_p1m1 = mass_flux[facexp1m1];
        //


        const RealType lag_node_mass_ooo = 0.25 * (cell_mass_m1m1 + cell_mass_p1m1 + cell_mass_m1p1 +  cell_mass_p1p1);
        assert(0.0 < lag_node_mass_ooo);
        const RealType out_node_mass_ooo = 0.25 * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
        assert(0.0 < out_node_mass_ooo);

        const RealType prev_dual_mass_flux = 0.25 * (mass_flux_m1m1 + mass_flux_oom1 + mass_flux_m1p1 + mass_flux_oop1);
	const RealType next_dual_mass_flux = 0.25 * (mass_flux_p1m1 + mass_flux_oom1 + mass_flux_p1p1 + mass_flux_oop1); 

const RealType dx_lag_m1 =  dx + dt * (u_oo - u_m1);
const RealType dx_lag_p1 =  dx + dt * (u_p1 - u_oo);

const RealType dt_uf_m1 = 0.5 * dt * (u_oo + u_m1);
const RealType dt_uf_p1 = 0.5 * dt * (u_oo + u_p1);

RealType nodal_reconstruct_variable_o2_m1;
RealType nodal_reconstruct_variable_o2_p1 ;

if (prev_dual_mass_flux > 0 ) {
  nodal_reconstruct_variable_o2_m1 = in_variable_m1 + 0.5 * gradient_variable_m1 * (dx_lag_m1 - dt_uf_m1);
 } else {
  nodal_reconstruct_variable_o2_m1 = in_variable_oo + 0.5 * gradient_variable_oo * ( -dx_lag_m1 - dt_uf_m1);
 } 

if (next_dual_mass_flux > 0 ) {
  nodal_reconstruct_variable_o2_p1 = in_variable_oo + 0.5 * gradient_variable_oo *(dx_lag_p1 - dt_uf_p1);
 } else {
  nodal_reconstruct_variable_o2_p1 = in_variable_p1 + 0.5 * gradient_variable_p1 * ( -dx_lag_p1 - dt_uf_p1);
 }

const RealType moment_flux_prev = nodal_reconstruct_variable_o2_m1 * prev_dual_mass_flux ;
const RealType moment_flux_next = nodal_reconstruct_variable_o2_p1 * next_dual_mass_flux ;
 
          
	const RealType in_moment_ooo = lag_node_mass_ooo * in_variable_oo; 
       
	const RealType out_moment_ooo = 
	  in_moment_ooo + moment_flux_prev - moment_flux_next;

	const RealType out_nodal_variable =  out_moment_ooo / out_node_mass_ooo ;

RealType moment_flux_sw, moment_flux_se, moment_flux_nw, moment_flux_ne; 

index_t node_prev_prev, node_prev_next, node_next_prev, node_next_next;
RealType delta_x_prev_prev, delta_x_prev_next, delta_x_next_prev, delta_x_next_next;
RealType delta_y_prev_prev, delta_y_prev_next, delta_y_next_prev, delta_y_next_next;
RealType delta_x_c, delta_y_c;
RealType rho_dual_corner;

RealType one_quarter = 0.25;
RealType half = 0.5;

RealType mass_flux_prev_prev = 0.0;
RealType mass_flux_prev_next = 0.0;
RealType mass_flux_next_prev = 0.0;
RealType mass_flux_next_next = 0.0;
RealType mass_flux_dual_corner = 0.0;

const RealType out_mass_m1m1 = out_cell_mass[cell_m1m1];
const RealType out_mass_p1m1 = out_cell_mass[cell_p1m1];
const RealType out_mass_m1p1 = out_cell_mass[cell_m1p1];
const RealType out_mass_p1p1 = out_cell_mass[cell_p1p1];

const RealType out_node_mass_ooo = one_quarter * (out_mass_m1m1 + out_mass_p1m1 + out_mass_m1p1 + out_mass_p1p1);
assert(0.0 < out_node_mass_ooo);


//SW corner
    { node_prev_prev = node_m1m1;
      node_prev_next = node_m1oo;
      node_next_prev = node_oom1;
      node_next_next = node_oooo;

      if ( (in_u[node_prev_prev] * in_v[node_prev_prev]) > 0.) {
	mass_flux_prev_prev = mass_corner_fluxes[node_prev_prev] * sign(in_u[node_prev_prev]);
      }
      if ( (in_u[node_prev_next] * in_v[node_prev_next]) > 0.) {
	mass_flux_prev_next = mass_corner_fluxes[node_prev_next] * sign(in_u[node_prev_next]);
      }
      if ( (in_u[node_next_prev] * in_v[node_next_prev]) > 0.) {
	mass_flux_next_prev = mass_corner_fluxes[node_next_prev] * sign(in_u[node_next_prev]);
      }
      if ( (in_u[node_next_next] * in_v[node_next_next]) > 0.) {
	mass_flux_next_next = mass_corner_fluxes[node_next_next] * sign(in_u[node_next_next]);
      }

      mass_flux_dual_corner = one_quarter * (mass_flux_prev_prev + mass_flux_prev_next + mass_flux_next_prev + mass_flux_next_next);

      moment_flux_sw = half * ( (mass_flux_dual_corner - fabs(mass_flux_dual_corner)) * in_variable[node_next_next] + (mass_flux_dual_corner + fabs(mass_flux_dual_corner)) * in_variable[node_prev_prev] );
      
    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;
mass_flux_dual_corner = 0.0;

//SE corner
    { node_prev_prev = node_oom1;
      node_prev_next = node_oooo;
      node_next_prev = node_p1m1;
      node_next_next = node_p1oo;

      if ( (in_u[node_prev_prev] * in_v[node_prev_prev]) < 0.) {
	mass_flux_prev_prev = mass_corner_fluxes[node_prev_prev] * sign(in_v[node_prev_prev]);
      }
      if ( (in_u[node_prev_next] * in_v[node_prev_next]) < 0.) {
	mass_flux_prev_next = mass_corner_fluxes[node_prev_next] * sign(in_v[node_prev_next]);
      }
      if ( (in_u[node_next_prev] * in_v[node_next_prev]) < 0.) {
	mass_flux_next_prev = mass_corner_fluxes[node_next_prev] * sign(in_v[node_next_prev]);
      }
      if ( (in_u[node_next_next] * in_v[node_next_next]) < 0.) {
	mass_flux_next_next = mass_corner_fluxes[node_next_next] * sign(in_v[node_next_next]);
      }

      mass_flux_dual_corner = one_quarter * (mass_flux_prev_prev + mass_flux_prev_next + mass_flux_next_prev + mass_flux_next_next);

      moment_flux_se = half * ( (mass_flux_dual_corner - fabs(mass_flux_dual_corner)) * in_variable[node_prev_next] + (mass_flux_dual_corner + fabs(mass_flux_dual_corner)) * in_variable[node_next_prev] );

    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;
mass_flux_dual_corner = 0.0;

//NE corner
    { node_prev_prev = node_oooo;
      node_prev_next = node_oop1;
      node_next_prev = node_p1oo;
      node_next_next = node_p1p1;

      if ( (in_u[node_prev_prev] * in_v[node_prev_prev]) > 0.) {
	mass_flux_prev_prev = - mass_corner_fluxes[node_prev_prev] * sign(in_u[node_prev_prev]);
      }
      if ( (in_u[node_prev_next] * in_v[node_prev_next]) > 0.) {
	mass_flux_prev_next = - mass_corner_fluxes[node_prev_next] * sign(in_u[node_prev_next]);
      }
      if ( (in_u[node_next_prev] * in_v[node_next_prev]) > 0.) {
	mass_flux_next_prev = - mass_corner_fluxes[node_next_prev] * sign(in_u[node_next_prev]);
      }
      if ( (in_u[node_next_next] * in_v[node_next_next]) > 0.) {
	mass_flux_next_next = - mass_corner_fluxes[node_next_next] * sign(in_u[node_next_next]);
      }

      mass_flux_dual_corner = one_quarter * (mass_flux_prev_prev + mass_flux_prev_next + mass_flux_next_prev + mass_flux_next_next);

      moment_flux_ne = half * ( (mass_flux_dual_corner - fabs(mass_flux_dual_corner)) * in_variable[node_prev_prev] + (mass_flux_dual_corner + fabs(mass_flux_dual_corner)) * in_variable[node_next_next] );

    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;
mass_flux_dual_corner = 0.0;

//NW corner
    { node_prev_prev = node_m1oo;
      node_prev_next = node_m1p1;
      node_next_prev = node_oooo;
      node_next_next = node_oop1;

      if ( (in_u[node_prev_prev] * in_v[node_prev_prev]) < 0.) {
	mass_flux_prev_prev = - mass_corner_fluxes[node_prev_prev] * sign(in_v[node_prev_prev]);
      }
      if ( (in_u[node_prev_next] * in_v[node_prev_next]) < 0.) {
	mass_flux_prev_next = - mass_corner_fluxes[node_prev_next] * sign(in_v[node_prev_next]);
      }
      if ( (in_u[node_next_prev] * in_v[node_next_prev]) < 0.) {
	mass_flux_next_prev = - mass_corner_fluxes[node_next_prev] * sign(in_v[node_next_prev]);
      }
      if ( (in_u[node_next_next] * in_v[node_next_next]) < 0.) {
	mass_flux_next_next = - mass_corner_fluxes[node_next_next] * sign(in_v[node_next_next]);
      }

      mass_flux_dual_corner = one_quarter * (mass_flux_prev_prev + mass_flux_prev_next + mass_flux_next_prev + mass_flux_next_next);

      moment_flux_nw = half * ( (mass_flux_dual_corner - fabs(mass_flux_dual_corner)) * in_variable[node_next_prev] + (mass_flux_dual_corner + fabs(mass_flux_dual_corner)) * in_variable[node_prev_next] );

    }


/*
//SW corner
    { node_prev_prev = node_m1m1;
      node_prev_next = node_m1oo;
      node_next_prev = node_oom1;
      node_next_next = node_oooo;

      delta_x_prev_prev = in_u[node_prev_prev] * dt;
      delta_x_prev_next = in_u[node_prev_next] * dt;
      delta_x_next_prev = in_u[node_next_prev] * dt;
      delta_x_next_next = in_u[node_next_next] * dt;
      delta_y_prev_prev = in_v[node_prev_prev] * dt;
      delta_y_prev_next = in_v[node_prev_next] * dt;
      delta_y_next_prev = in_v[node_next_prev] * dt;
      delta_y_next_next = in_v[node_next_next] * dt;

      if ((mass_corner_fluxes[node_prev_prev] == 0.) && (mass_corner_fluxes[node_prev_next] == 0.) && (mass_corner_fluxes[node_next_prev] == 0.) && (mass_corner_fluxes[node_next_next] == 0.)) {
	rho_dual_corner = 0.;
      } else {
	rho_dual_corner = (mass_corner_fluxes[node_prev_prev] + mass_corner_fluxes[node_prev_next] + mass_corner_fluxes[node_next_prev] + mass_corner_fluxes[node_next_next]) / (fabs(delta_x_prev_prev * delta_y_prev_prev) + fabs(delta_x_prev_next * delta_y_prev_next) + fabs(delta_x_next_prev * delta_y_next_prev) + fabs(delta_x_next_next * delta_y_next_next));
      }

      delta_x_c = one_quarter * (delta_x_prev_prev + delta_x_prev_next + delta_x_next_prev + delta_x_next_next);
      delta_y_c = one_quarter * (delta_y_prev_prev + delta_y_prev_next + delta_y_next_prev + delta_y_next_next);
     
      if ((delta_x_c > 0.) && (delta_y_c > 0.)) {
	moment_flux_sw = fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_prev_prev];
      
      } else if ((delta_x_c < 0.) && (delta_y_c < 0.)) {
	moment_flux_sw = - fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_next_next];
      
      } else {
	moment_flux_sw = 0.;
      }
      
    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;

//SE corner
    { node_prev_prev = node_oom1;
      node_prev_next = node_oooo;
      node_next_prev = node_p1m1;
      node_next_next = node_p1oo;

      delta_x_prev_prev = in_u[node_prev_prev] * dt;
      delta_x_prev_next = in_u[node_prev_next] * dt;
      delta_x_next_prev = in_u[node_next_prev] * dt;
      delta_x_next_next = in_u[node_next_next] * dt;
      delta_y_prev_prev = in_v[node_prev_prev] * dt;
      delta_y_prev_next = in_v[node_prev_next] * dt;
      delta_y_next_prev = in_v[node_next_prev] * dt;
      delta_y_next_next = in_v[node_next_next] * dt;

      if ((mass_corner_fluxes[node_prev_prev] == 0.) && (mass_corner_fluxes[node_prev_next] == 0.) && (mass_corner_fluxes[node_next_prev] == 0.) && (mass_corner_fluxes[node_next_next] == 0.)) {
	rho_dual_corner = 0.;
      } else {
	rho_dual_corner = (mass_corner_fluxes[node_prev_prev] + mass_corner_fluxes[node_prev_next] + mass_corner_fluxes[node_next_prev] + mass_corner_fluxes[node_next_next]) / (fabs(delta_x_prev_prev * delta_y_prev_prev) + fabs(delta_x_prev_next * delta_y_prev_next) + fabs(delta_x_next_prev * delta_y_next_prev) + fabs(delta_x_next_next * delta_y_next_next));
      }

      delta_x_c = one_quarter * (delta_x_prev_prev + delta_x_prev_next + delta_x_next_prev + delta_x_next_next);
      delta_y_c = one_quarter * (delta_y_prev_prev + delta_y_prev_next + delta_y_next_prev + delta_y_next_next);

      if ((delta_x_c < 0.) && (delta_y_c > 0.)) {
	moment_flux_se = fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_next_prev];
      
      } else if ((delta_x_c > 0.) && (delta_y_c < 0.)) {
	moment_flux_se = - fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_prev_next];
      
      } else {
	moment_flux_se = 0.;
      }

    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;

//NE corner
    { node_prev_prev = node_oooo;
      node_prev_next = node_oop1;
      node_next_prev = node_p1oo;
      node_next_next = node_p1p1;

      delta_x_prev_prev = in_u[node_prev_prev] * dt;
      delta_x_prev_next = in_u[node_prev_next] * dt;
      delta_x_next_prev = in_u[node_next_prev] * dt;
      delta_x_next_next = in_u[node_next_next] * dt;
      delta_y_prev_prev = in_v[node_prev_prev] * dt;
      delta_y_prev_next = in_v[node_prev_next] * dt;
      delta_y_next_prev = in_v[node_next_prev] * dt;
      delta_y_next_next = in_v[node_next_next] * dt;

      if ((mass_corner_fluxes[node_prev_prev] == 0.) && (mass_corner_fluxes[node_prev_next] == 0.) && (mass_corner_fluxes[node_next_prev] == 0.) && (mass_corner_fluxes[node_next_next] == 0.)) {
	rho_dual_corner = 0.;
      } else {
	rho_dual_corner = (mass_corner_fluxes[node_prev_prev] + mass_corner_fluxes[node_prev_next] + mass_corner_fluxes[node_next_prev] + mass_corner_fluxes[node_next_next]) / (fabs(delta_x_prev_prev * delta_y_prev_prev) + fabs(delta_x_prev_next * delta_y_prev_next) + fabs(delta_x_next_prev * delta_y_next_prev) + fabs(delta_x_next_next * delta_y_next_next));
      }

      delta_x_c = one_quarter * (delta_x_prev_prev + delta_x_prev_next + delta_x_next_prev + delta_x_next_next);
      delta_y_c = one_quarter * (delta_y_prev_prev + delta_y_prev_next + delta_y_next_prev + delta_y_next_next);
 
      if ((delta_x_c < 0.) && (delta_y_c < 0.)) {
	moment_flux_ne = fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_next_next];
      
      } else if ((delta_x_c > 0.) && (delta_y_c > 0.)) {
	moment_flux_ne = - fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_prev_prev];
      
      } else {
	moment_flux_ne = 0.;
      }

    }

mass_flux_prev_prev = 0.0;
mass_flux_prev_next = 0.0;
mass_flux_next_prev = 0.0;
mass_flux_next_next = 0.0;

//NW corner
    { node_prev_prev = node_m1oo;
      node_prev_next = node_m1p1;
      node_next_prev = node_oooo;
      node_next_next = node_oop1;

      delta_x_prev_prev = in_u[node_prev_prev] * dt;
      delta_x_prev_next = in_u[node_prev_next] * dt;
      delta_x_next_prev = in_u[node_next_prev] * dt;
      delta_x_next_next = in_u[node_next_next] * dt;
      delta_y_prev_prev = in_v[node_prev_prev] * dt;
      delta_y_prev_next = in_v[node_prev_next] * dt;
      delta_y_next_prev = in_v[node_next_prev] * dt;
      delta_y_next_next = in_v[node_next_next] * dt;

      if ((mass_corner_fluxes[node_prev_prev] == 0.) && (mass_corner_fluxes[node_prev_next] == 0.) && (mass_corner_fluxes[node_next_prev] == 0.) && (mass_corner_fluxes[node_next_next] == 0.)) {
	rho_dual_corner = 0.;
      } else {
	rho_dual_corner = (mass_corner_fluxes[node_prev_prev] + mass_corner_fluxes[node_prev_next] + mass_corner_fluxes[node_next_prev] + mass_corner_fluxes[node_next_next]) / (fabs(delta_x_prev_prev * delta_y_prev_prev) + fabs(delta_x_prev_next * delta_y_prev_next) + fabs(delta_x_next_prev * delta_y_next_prev) + fabs(delta_x_next_next * delta_y_next_next));
      }

      delta_x_c = one_quarter * (delta_x_prev_prev + delta_x_prev_next + delta_x_next_prev + delta_x_next_next);
      delta_y_c = one_quarter * (delta_y_prev_prev + delta_y_prev_next + delta_y_next_prev + delta_y_next_next);

      if ((delta_x_c > 0.) && (delta_y_c < 0.)) {
	moment_flux_nw = fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_prev_next];
      
      } else if ((delta_x_c < 0.) && (delta_y_c > 0.)) {
	moment_flux_nw = - fabs(delta_x_c * delta_y_c) * rho_dual_corner * in_variable[node_next_prev];
      
      } else {
	moment_flux_nw = 0.;
      }

    }
*/

const RealType total_variable_flux = (moment_flux_sw + moment_flux_se + moment_flux_ne + moment_flux_nw) / out_node_mass_ooo;

//printf("\n At node %d,  flux sw = %f,   flux se = %f,   flux ne = %f,   flux nw = %f,  total variable flux = %f \n",node_oooo,moment_flux_sw,moment_flux_se,moment_flux_ne,moment_flux_nw,total_variable_flux); 

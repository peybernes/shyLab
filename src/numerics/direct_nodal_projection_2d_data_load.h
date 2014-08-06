      const RealType in_vx_ooo = in_vx[node_ooo];
      const RealType in_vx_m1o = in_vx[node_m1o];
      const RealType in_vx_p1o = in_vx[node_p1o];
      const RealType in_vx_om1 = in_vx[node_om1];
      const RealType in_vx_op1 = in_vx[node_op1];

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
      const RealType mass_flux_y_p1p1 = mass_flux_y[faceyp1p1];
      const RealType mass_flux_y_p1oo = mass_flux_y[faceyp1oo];
      const RealType mass_flux_y_p1m1 = mass_flux_y[faceyp1m1];

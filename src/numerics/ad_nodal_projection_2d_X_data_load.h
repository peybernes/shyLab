      const RealType in_vx_oo = in_vx[node_ooo];
      const RealType in_vx_m1 = in_vx[node_m1o];
      const RealType in_vx_p1 = in_vx[node_p1o];

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

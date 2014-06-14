      const RealType in_vx_oo = in_vy[node_ooo];
      const RealType in_vx_m1 = in_vy[node_om1];
      const RealType in_vx_p1 = in_vy[node_op1];

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
      const RealType mass_flux_p1p1 = mass_flux[faceyp1p1];
      const RealType mass_flux_p1oo = mass_flux[faceyp1oo];
      const RealType mass_flux_p1m1 = mass_flux[faceyp1m1];      

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = (nx * iy) + ix;

      const RealType vol_fraction_11 = out_cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
      
      out_rho[cell_ooo] = vol_fraction_1 * out_rho_1[cell_ooo] + vol_fraction_2 * out_rho_2[cell_ooo];
     
      RealType out_c1 = out_cell_mass_1[cell_ooo] / out_cell_mass[cell_ooo];
      RealType out_c2 = out_cell_mass_2[cell_ooo] / out_cell_mass[cell_ooo];

      out_e[cell_ooo] = out_c1 * out_e_1[cell_ooo] + out_c2 * out_e_2[cell_ooo];
      
      
    }
  }

RealType cv_1 = 8.31 / (gamma_1 - 1.0);
RealType cv_2 = 8.31 / (gamma_2 - 1.0);

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = (nx * iy) + ix;

      out_rho[cell_ooo] =  out_cell_mass[cell_ooo] / cell_volumes[cell_ooo];

      out_c_2[cell_ooo] = 1.0 - out_c_1[cell_ooo];

      out_e_2[cell_ooo] = out_e[cell_ooo] / ( out_c_1[cell_ooo]*cv_1/cv_2 + out_c_2[cell_ooo] );

      out_e_1[cell_ooo] = cv_1 / cv_2 * out_e_2[cell_ooo];

      out_rho_2[cell_ooo] = out_rho[cell_ooo] * ( out_c_1[cell_ooo]*(gamma_1-1.0)*cv_1/((gamma_2-1.0)*cv_2) + out_c_2[cell_ooo] );

      out_rho_1[cell_ooo] = out_rho_2[cell_ooo] * (gamma_1-1.0)*cv_1/((gamma_2-1.0)*cv_2);

    }
  }

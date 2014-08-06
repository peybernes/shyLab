#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int cell_ooo = (nx * iy) + ix;
      
      out_rho[cell_ooo] = out_cell_mass[cell_ooo] / cell_volumes[cell_ooo];
  
    }
  }

#pragma omp parallel for
  for (int iy = 0; iy < ny + 1; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const int face_ooo = nx * iy + ix;      

      mass_flux_y[face_ooo] = mass_flux_1_y[face_ooo] + mass_flux_2_y[face_ooo];
      
    }
  }


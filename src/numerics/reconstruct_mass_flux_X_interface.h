#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx + 1; ++ix) {

      const int face_ooo = ((nx + 1) * iy) + ix;      

      mass_flux_x[face_ooo] = mass_flux_1_x[face_ooo] + mass_flux_2_x[face_ooo];
      
    }
  }

#pragma omp parallel for
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    for (index_t ix = 0; ix < nx + 1; ++ix) {

      const index_t node_ooo = ((nx + 1) * iy) + ix;      

      mass_corner_fluxes[node_ooo] = mass_1_corner_fluxes[node_ooo] + mass_2_corner_fluxes[node_ooo];
      
    }
  }

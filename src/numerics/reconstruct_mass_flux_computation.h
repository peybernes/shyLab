const RealType face_flux = volume_fluxes[face_ooo];

const RealType mass_flux_ooo =
  (0.5 * (face_flux + fabs(face_flux)) * prev_cell_variable) + 
  (0.5 * (face_flux - fabs(face_flux)) * next_cell_variable);


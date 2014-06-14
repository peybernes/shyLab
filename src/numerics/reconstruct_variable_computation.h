const RealType face_flux = mass_flux[face_ooo];

const RealType variable_flux_ooo =
  (0.5 * (face_flux + fabs(face_flux)) * prev_cell_variable) + 
  (0.5 * (face_flux - fabs(face_flux)) * next_cell_variable);


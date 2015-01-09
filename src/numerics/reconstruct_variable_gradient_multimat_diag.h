const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

const RealType vol_fraction_m1m1 = in_cell_volumic_fraction[cell_m1m1];
const RealType vol_fraction_m1p1 = in_cell_volumic_fraction[cell_m1p1];
const RealType vol_fraction_p1m1 = in_cell_volumic_fraction[cell_p1m1];
const RealType vol_fraction_p1p1 = in_cell_volumic_fraction[cell_p1p1];
const RealType vol_fraction_ooo = in_cell_volumic_fraction[cell_ooo];

const RealType variable_1_m1m1 = lag_variable_1[cell_m1m1]; // 1 load
const RealType variable_1_m1p1 = lag_variable_1[cell_m1p1]; // 1 load
const RealType variable_1_p1m1 = lag_variable_1[cell_p1m1]; // 1 load
const RealType variable_1_p1p1 = lag_variable_1[cell_p1p1]; // 1 load
const RealType variable_1_ooo = lag_variable_1[cell_ooo]; // 0 load

const RealType variable_2_m1m1 = lag_variable_2[cell_m1m1]; // 1 load
const RealType variable_2_m1p1 = lag_variable_2[cell_m1p1]; // 1 load
const RealType variable_2_p1m1 = lag_variable_2[cell_p1m1]; // 1 load
const RealType variable_2_p1p1 = lag_variable_2[cell_p1p1]; // 1 load
const RealType variable_2_ooo = lag_variable_2[cell_ooo]; // 0 load

const RealType d_vol_1_m1m1 = volume_fluxes_1_corner[node_m1m1];
const RealType d_vol_1_m1p1 = volume_fluxes_1_corner[node_m1p1];
const RealType d_vol_1_p1m1 = volume_fluxes_1_corner[node_p1m1];
const RealType d_vol_1_p1p1 = volume_fluxes_1_corner[node_p1p1];
const RealType d_vol_1_mm1mm1 = volume_fluxes_1_corner[node_mm1mm1];
const RealType d_vol_1_mm1pp1 = volume_fluxes_1_corner[node_mm1pp1];
const RealType d_vol_1_pp1mm1 = volume_fluxes_1_corner[node_pp1mm1];
const RealType d_vol_1_pp1pp1 = volume_fluxes_1_corner[node_pp1pp1];

const RealType d_vol_2_m1m1 = volume_fluxes_2_corner[node_m1m1];
const RealType d_vol_2_m1p1 = volume_fluxes_2_corner[node_m1p1];
const RealType d_vol_2_p1m1 = volume_fluxes_2_corner[node_p1m1];
const RealType d_vol_2_p1p1 = volume_fluxes_2_corner[node_p1p1];
const RealType d_vol_2_mm1mm1 = volume_fluxes_2_corner[node_mm1mm1];
const RealType d_vol_2_mm1pp1 = volume_fluxes_2_corner[node_mm1pp1];
const RealType d_vol_2_pp1mm1 = volume_fluxes_2_corner[node_pp1mm1];
const RealType d_vol_2_pp1pp1 = volume_fluxes_2_corner[node_pp1pp1];

const RealType diag_m1m1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_p1p1] - predicted_u[node_mm1mm1]) * dx + (predicted_v[node_p1p1] - predicted_v[node_mm1mm1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_p1m1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_pp1mm1] - predicted_u[node_m1p1]) * dx - (predicted_v[node_pp1mm1] - predicted_v[node_m1p1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_m1p1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_p1m1] - predicted_u[node_mm1pp1]) * dx - (predicted_v[node_p1m1] - predicted_v[node_mm1pp1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_p1p1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_pp1pp1] - predicted_u[node_m1m1]) * dx + (predicted_v[node_pp1pp1] - predicted_v[node_m1m1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 

RealType limited_grad_variable_1_diag, limited_grad_variable_1_antidiag, limited_grad_variable_2_diag, limited_grad_variable_2_antidiag;

// If there is a contact discontinuity around cell_ooo, we go back to order 1
if ( (d_vol_2_mm1mm1 == 0) && (d_vol_2_m1m1 == 0) && (d_vol_2_p1p1 == 0) && (d_vol_2_pp1pp1 == 0) && (vol_fraction_m1m1 == 1.0) && (vol_fraction_ooo == 1.0) && (vol_fraction_p1p1 == 1.0) ) {

  const RealType grad_1_m1m1 = (variable_1_ooo - variable_1_m1m1 ) / diag_m1m1;
  const RealType grad_1_p1p1 = (variable_1_p1p1 - variable_1_ooo ) / diag_p1p1;

  limited_grad_variable_1_diag = VanLeerLimiter(grad_1_m1m1, grad_1_p1p1);

 } else {

  limited_grad_variable_1_diag = 0.;
 }
 
if ( (d_vol_2_mm1pp1 == 0) && (d_vol_2_m1p1 == 0) && (d_vol_2_p1m1 == 0) && (d_vol_2_pp1mm1 == 0) && (vol_fraction_m1p1 == 1.0) && (vol_fraction_ooo == 1.0) && (vol_fraction_p1m1 == 1.0) ) { 

  const RealType grad_1_p1m1 = (variable_1_p1m1 - variable_1_ooo ) / diag_p1m1;
  const RealType grad_1_m1p1 = (variable_1_ooo - variable_1_m1p1 ) / diag_m1p1;

  limited_grad_variable_1_antidiag = VanLeerLimiter(grad_1_m1p1, grad_1_p1m1);

 } else {

  limited_grad_variable_1_antidiag = 0.;
 }

if ( (d_vol_1_mm1mm1 == 0) && (d_vol_1_m1m1 == 0) && (d_vol_1_p1p1 == 0) && (d_vol_1_pp1pp1 == 0) && (vol_fraction_m1m1 == 0.0) && (vol_fraction_ooo == 0.0) && (vol_fraction_p1p1 == 0.0) ) {

  const RealType grad_2_m1m1 = (variable_2_ooo - variable_2_m1m1 ) / diag_m1m1;
  const RealType grad_2_p1p1 = (variable_2_p1p1 - variable_2_ooo ) / diag_p1p1;

  limited_grad_variable_2_diag = VanLeerLimiter(grad_2_m1m1, grad_2_p1p1);

 } else {

  limited_grad_variable_2_diag = 0.;
 }
 
if ( (d_vol_1_mm1pp1 == 0) && (d_vol_1_m1p1 == 0) && (d_vol_1_p1m1 == 0) && (d_vol_1_pp1mm1 == 0) && (vol_fraction_m1p1 == 0.0) && (vol_fraction_ooo == 0.0) && (vol_fraction_p1m1 == 0.0) ) { 

  const RealType grad_2_p1m1 = (variable_2_p1m1 - variable_2_ooo ) / diag_p1m1;
  const RealType grad_2_m1p1 = (variable_2_ooo - variable_2_m1p1 ) / diag_m1p1;

  limited_grad_variable_2_antidiag = VanLeerLimiter(grad_2_m1p1, grad_2_p1m1);

 } else {

  limited_grad_variable_2_antidiag = 0.;
 }

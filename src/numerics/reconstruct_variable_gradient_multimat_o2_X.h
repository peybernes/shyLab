const RealType vol_fraction_m1o = in_cell_volumic_fraction[cell_m1o];
const RealType vol_fraction_ooo = in_cell_volumic_fraction[cell_ooo];
const RealType vol_fraction_p1o = in_cell_volumic_fraction[cell_p1o];
      
const RealType variable_1_m1o = lag_variable_1[cell_m1o];
const RealType variable_1_ooo = lag_variable_1[cell_ooo];
const RealType variable_1_p1o = lag_variable_1[cell_p1o];

const RealType variable_2_m1o = lag_variable_2[cell_m1o];
const RealType variable_2_ooo = lag_variable_2[cell_ooo];
const RealType variable_2_p1o = lag_variable_2[cell_p1o];

const RealType d_vol_1_m2o = volume_fluxes_1[face_m2o];
const RealType d_vol_1_m1o = volume_fluxes_1[face_m1o];
const RealType d_vol_1_p1o = volume_fluxes_1[face_p1o];
const RealType d_vol_1_p2o = volume_fluxes_1[face_p2o];

const RealType d_vol_2_m2o = volume_fluxes_2[face_m2o];
const RealType d_vol_2_m1o = volume_fluxes_2[face_m1o];
const RealType d_vol_2_p1o = volume_fluxes_2[face_p1o];
const RealType d_vol_2_p2o = volume_fluxes_2[face_p2o];

RealType limited_grad_variable_1;
RealType limited_grad_variable_2;

// If there is a contact discontinuity around cell_ooo, we go back to order 1
if ( (d_vol_2_m2o == 0) && (d_vol_2_m1o == 0) && (d_vol_2_p1o == 0) && (d_vol_2_p2o == 0) && (vol_fraction_m1o == 1.0) && (vol_fraction_ooo == 1.0) && (vol_fraction_p1o == 1.0) ) {

  const RealType grad_1_m1o = (variable_1_ooo - variable_1_m1o ) /
    (dx + (d_vol_1_p1o - d_vol_1_m2o) / dy);

  const RealType grad_1_p1o = (variable_1_p1o - variable_1_ooo ) /
    (dx + (d_vol_1_p2o - d_vol_1_m1o) / dy);

  limited_grad_variable_1 = VanLeerLimiter(grad_1_m1o, grad_1_p1o);

 } else {

  limited_grad_variable_1 = 0.0;      
 }

if ( (d_vol_1_m2o == 0) && (d_vol_1_m1o == 0) && (d_vol_1_p1o == 0) && (d_vol_1_p2o == 0) && (vol_fraction_m1o == 0.0) && (vol_fraction_ooo == 0.0) && (vol_fraction_p1o == 0.0) ) {

  const RealType grad_2_m1o = (variable_2_ooo - variable_2_m1o ) /
    (dx + (d_vol_2_p1o - d_vol_2_m2o) / dy);

  const RealType grad_2_p1o = (variable_2_p1o - variable_2_ooo ) /
    (dx + (d_vol_2_p2o - d_vol_2_m1o) / dy);

  limited_grad_variable_2 = VanLeerLimiter(grad_2_m1o, grad_2_p1o);

 } else {

  limited_grad_variable_2 = 0.0;     
 }

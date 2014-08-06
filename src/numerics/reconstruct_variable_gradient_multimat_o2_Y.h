const RealType vol_fraction_om1 = in_cell_volumic_fraction[cell_om1];
const RealType vol_fraction_ooo = in_cell_volumic_fraction[cell_ooo];
const RealType vol_fraction_op1 = in_cell_volumic_fraction[cell_op1];

const RealType variable_1_om1 = lag_variable_1[cell_om1];
const RealType variable_1_ooo = lag_variable_1[cell_ooo];
const RealType variable_1_op1 = lag_variable_1[cell_op1];

const RealType variable_2_om1 = lag_variable_2[cell_om1];
const RealType variable_2_ooo = lag_variable_2[cell_ooo];
const RealType variable_2_op1 = lag_variable_2[cell_op1];

const RealType d_vol_1_om2 = volume_fluxes_1[face_om2];
const RealType d_vol_1_om1 = volume_fluxes_1[face_om1];
const RealType d_vol_1_op1 = volume_fluxes_1[face_op1];
const RealType d_vol_1_op2 = volume_fluxes_1[face_op2];

const RealType d_vol_2_om2 = volume_fluxes_2[face_om2];
const RealType d_vol_2_om1 = volume_fluxes_2[face_om1];
const RealType d_vol_2_op1 = volume_fluxes_2[face_op1];
const RealType d_vol_2_op2 = volume_fluxes_2[face_op2];

RealType limited_grad_variable_1;
RealType limited_grad_variable_2;

// If there is a contact discontinuity around cell_ooo, we go back to order 1
if ( (d_vol_2_om2 == 0) && (d_vol_2_om1 == 0) && (d_vol_2_op1 == 0) && (d_vol_2_op2 == 0) && (vol_fraction_om1 == 1.0) && (vol_fraction_ooo == 1.0) && (vol_fraction_op1 == 1.0) ) {

  const RealType grad_1_om1 = (variable_1_ooo - variable_1_om1 ) /
    (dx + (d_vol_1_op1 - d_vol_1_om2) / dy);

  const RealType grad_1_op1 = (variable_1_op1 - variable_1_ooo ) /
    (dx + (d_vol_1_op2 - d_vol_1_om1) / dy);

  limited_grad_variable_1 = MinmodLimiter(grad_1_om1, grad_1_op1);

 } else {

  limited_grad_variable_1 = 0.0;      
 }

if ( (d_vol_1_om2 == 0) && (d_vol_1_om1 == 0) && (d_vol_1_op1 == 0) && (d_vol_1_op2 == 0) && (vol_fraction_om1 == 0.0) && (vol_fraction_ooo == 0.0) && (vol_fraction_op1 == 0.0) ) {

  const RealType grad_2_om1 = (variable_2_ooo - variable_2_om1 ) /
    (dx + (d_vol_2_op1 - d_vol_2_om2) / dy);

  const RealType grad_2_op1 = (variable_2_op1 - variable_2_ooo ) /
    (dx + (d_vol_2_op2 - d_vol_2_om1) / dy);

  limited_grad_variable_2 = MinmodLimiter(grad_2_om1, grad_2_op1);

 } else {

  limited_grad_variable_2 = 0.0;     
 }   

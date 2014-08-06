const RealType variable_oooo = in_variable[cell_oooo];
const RealType variable_m1oo = in_variable[cell_m1oo];
const RealType variable_p1oo = in_variable[cell_p1oo];
const RealType variable_oom1 = in_variable[cell_oom1];
const RealType variable_oop1 = in_variable[cell_oop1];
const RealType variable_m1m1 = in_variable[cell_m1m1];
const RealType variable_p1m1 = in_variable[cell_p1m1];
const RealType variable_p1p1 = in_variable[cell_p1p1];
const RealType variable_m1p1 = in_variable[cell_m1p1];

const RealType var_m1o = 0.25*(variable_m1p1 + 2*variable_m1oo + variable_m1m1);
const RealType var_p1o = 0.25*(variable_p1p1 + 2*variable_p1oo + variable_p1m1);
const RealType var_om1 = 0.25*(variable_m1m1 + 2*variable_oom1 + variable_p1m1);
const RealType var_op1 = 0.25*(variable_m1p1 + 2*variable_oop1 + variable_p1p1);

const RealType gradient_var_x = (var_p1o - var_m1o) / (2*dx);
const RealType gradient_var_y = (var_op1 - var_om1) / (2*dy);
const RealType norm_gradient_var = sqrt(gradient_var_x * gradient_var_x + gradient_var_y * gradient_var_y);

const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

const RealType variable_m1m1 = lag_variable[cell_m1m1]; // 1 load
const RealType variable_m1p1 = lag_variable[cell_m1p1]; // 1 load
const RealType variable_p1m1 = lag_variable[cell_p1m1]; // 1 load
const RealType variable_p1p1 = lag_variable[cell_p1p1]; // 1 load
const RealType variable_ooo = lag_variable[cell_ooo]; // 0 load

const RealType diag_m1m1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_p1p1] - predicted_u[node_mm1mm1]) * dx + (predicted_v[node_p1p1] - predicted_v[node_mm1mm1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_p1m1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_pp1mm1] - predicted_u[node_m1p1]) * dx - (predicted_v[node_pp1mm1] - predicted_v[node_m1p1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_m1p1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_p1m1] - predicted_u[node_mm1pp1]) * dx - (predicted_v[node_p1m1] - predicted_v[node_mm1pp1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 
const RealType diag_p1p1 = std::sqrt(dx * dx + dy * dy) + ((predicted_u[node_pp1pp1] - predicted_u[node_m1m1]) * dx + (predicted_v[node_pp1pp1] - predicted_v[node_m1m1]) * dy) * dt / std::sqrt(dx * dx + dy * dy); 

const RealType grad_m1m1 = (variable_ooo - variable_m1m1 ) / diag_m1m1;
const RealType grad_p1m1 = (variable_p1m1 - variable_ooo ) / diag_p1m1;
const RealType grad_m1p1 = (variable_ooo - variable_m1p1 ) / diag_m1p1;
const RealType grad_p1p1 = (variable_p1p1 - variable_ooo ) / diag_p1p1;


const RealType limited_grad_variable_diag = VanLeerLimiter(grad_m1m1, grad_p1p1);
const RealType limited_grad_variable_antidiag = VanLeerLimiter(grad_m1p1, grad_p1m1);

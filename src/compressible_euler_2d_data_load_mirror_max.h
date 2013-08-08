// -*- C++-mode -*-
// Fix this ugliness :)

const RealType in_rho_m2o = in_rho[cell_m2o];
const RealType in_rho_m1o = in_rho[cell_m1o];
const RealType in_rho_ooo = in_rho[cell_ooo];
const RealType in_rho_p1o = (cell_p1o >= nx ? in_rho[cell_m1o] : in_rho[cell_p1o]);
const RealType in_rho_p2o = (cell_p2o >= nx ? in_rho[cell_m2o] : in_rho[cell_p2o]);

// Note the sign change.
const RealType in_velocity_x_m2o = in_velocity_x[cell_m2o];
const RealType in_velocity_x_m1o = in_velocity_x[cell_m1o];
const RealType in_velocity_x_ooo = in_velocity_x[cell_ooo];
const RealType in_velocity_x_p2o = (cell_p2o >= nx ? -in_velocity_x[cell_m2o] : in_velocity_x[cell_p2o]);
const RealType in_velocity_x_p1o = (cell_p1o >= nx ? -in_velocity_x[cell_m1o] : in_velocity_x[cell_p1o]);

const RealType in_velocity_y_m2o = in_velocity_y[cell_m2o];
const RealType in_velocity_y_m1o = in_velocity_y[cell_m1o];
const RealType in_velocity_y_ooo = in_velocity_y[cell_ooo];
const RealType in_velocity_y_p2o = (cell_p2o >= nx ? in_velocity_y[cell_m2o] : in_velocity_y[cell_p2o]);
const RealType in_velocity_y_p1o = (cell_p1o >= nx ? in_velocity_y[cell_m1o] : in_velocity_y[cell_p1o]);

const RealType in_total_energy_m2o = in_total_energy[cell_m2o];
const RealType in_total_energy_m1o = in_total_energy[cell_m1o];
const RealType in_total_energy_ooo = in_total_energy[cell_ooo];
const RealType in_total_energy_p2o = (cell_p2o >= nx ? in_total_energy[cell_m2o] : in_total_energy[cell_p2o]);
const RealType in_total_energy_p1o = (cell_p1o >= nx ? in_total_energy[cell_m1o] : in_total_energy[cell_p1o]);

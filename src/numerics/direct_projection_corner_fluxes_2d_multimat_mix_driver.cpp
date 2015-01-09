// -*- c++ -*- (for emacs users)
#include "direct_projection_corner_fluxes_2d_multimat_mix_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"

void DirectProjectionCornerFluxes2dMultimaterialMixDriver(//in
							  const std::string BoundaryConditions, 
							  const int nx,
							  const int ny,
							  const RealType dx,
							  const RealType dy,
							  const RealType dt,
							  const RealType halo_width,
							  const RealType gamma_1,
							  const RealType gamma_2,
							  const RealType pi_1,
							  const RealType pi_2,
							  const RealType* RESTRICT predicted_u,
							  const RealType* RESTRICT predicted_v,
							  const RealType* RESTRICT cell_volumes,
							  const RealType* RESTRICT e_lag,
							  const RealType* RESTRICT u_lag,
							  const RealType* RESTRICT v_lag,
							  const RealType* RESTRICT in_cell_mass,
							  const RealType* RESTRICT in_pressure,
							  const RealType* RESTRICT in_c_1,
							  const RealType* RESTRICT in_c_2,				   
							  //out
							  RealType* RESTRICT out_u,
							  RealType* RESTRICT out_v,
							  RealType* RESTRICT out_e,
							  RealType* RESTRICT out_cell_mass,
							  RealType* RESTRICT out_e_1,
							  RealType* RESTRICT out_e_2,
							  RealType* RESTRICT out_rho,
							  RealType* RESTRICT out_rho_1,
							  RealType* RESTRICT out_rho_2,
							  RealType* RESTRICT out_c_1,
							  RealType* RESTRICT out_c_2,
							  RealType* RESTRICT directional_lagrangian_volume,
							  RealType* RESTRICT directional_lagrangian_density,
							  RealType* RESTRICT volume_fluxes_x,
							  RealType* RESTRICT volume_fluxes_y,
							  RealType* RESTRICT volume_fluxes_corner,
							  RealType* RESTRICT mass_flux_x,
							  RealType* RESTRICT mass_flux_y,
							  RealType* RESTRICT mass_corner_fluxes,
							  index_t* RESTRICT sign_x_corner_fluxes,
							  index_t* RESTRICT sign_y_corner_fluxes,
							  RealType* RESTRICT energy_flux_x,
							  RealType* RESTRICT energy_flux_y,
							  RealType* RESTRICT energy_flux_corner,
							  RealType* RESTRICT concentration_flux_x,
							  RealType* RESTRICT concentration_flux_y,
							  RealType* RESTRICT concentration_flux_corner,
							  RealType* RESTRICT density_gradient,
							  RealType* RESTRICT density_gradient_y,
							  RealType* RESTRICT density_gradient_diag,
							  RealType* RESTRICT density_gradient_antidiag,
							  RealType* RESTRICT energy_gradient,
							  RealType* RESTRICT energy_gradient_y,
							  RealType* RESTRICT energy_gradient_diag,
							  RealType* RESTRICT energy_gradient_antidiag,
							  RealType* RESTRICT concentration_gradient,
							  RealType* RESTRICT concentration_gradient_y,
							  RealType* RESTRICT concentration_gradient_diag,
							  RealType* RESTRICT concentration_gradient_antidiag,
							  RealType* RESTRICT reconstructed_density_faces_x,
							  RealType* RESTRICT reconstructed_density_faces_y,
							  RealType* RESTRICT reconstructed_energy_faces_x,
							  RealType* RESTRICT reconstructed_energy_faces_y,
							  RealType* RESTRICT reconstructed_concentration_faces_x,
							  RealType* RESTRICT reconstructed_concentration_faces_y,
							  RealType* RESTRICT gradient_u,
							  RealType* RESTRICT gradient_u_y,
							  RealType* RESTRICT gradient_u_diag,
							  RealType* RESTRICT gradient_u_antidiag,
							  RealType* RESTRICT gradient_v,
							  RealType* RESTRICT gradient_v_y,
							  RealType* RESTRICT gradient_v_diag,
							  RealType* RESTRICT gradient_v_antidiag,
							  //timing
							  std::vector<RealType>& time_compute_volume_fluxes_X,
							  std::vector<RealType>& time_gradient_X,
							  std::vector<RealType>& time_mass_reconstruct_o2_X,
							  std::vector<RealType>& time_project_mass_X,
							  std::vector<RealType>& time_reconstruct_energy_o2_X,
							  std::vector<RealType>& time_project_energy_X,
							  std::vector<RealType>& time_gradient_nodal_X,
							  std::vector<RealType>& time_project_nodal_velocity_X,
							  std::vector<RealType>& time_compute_volume_fluxes_Y,
							  std::vector<RealType>& time_gradient_Y,
							  std::vector<RealType>& time_mass_reconstruct_o2_Y,
							  std::vector<RealType>& time_project_mass_Y,
							  std::vector<RealType>& time_reconstruct_energy_o2_Y,
							  std::vector<RealType>& time_project_energy_Y,
							  std::vector<RealType>& time_gradient_nodal_Y,
							  std::vector<RealType>& time_project_nodal_velocity_Y,
							  std::vector<RealType> time_periodic_boundary) {


  //for timing
  struct timespec time_begin;
  struct timespec time_end;
  //

  //========================
  //    Projection _ mass .
  //========================
    
      
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianQuantitiesCornerFluxesXY(nx, ny, dt, dx, dy, predicted_u, predicted_v, in_cell_mass,
						       volume_fluxes_x, volume_fluxes_y, volume_fluxes_corner, sign_x_corner_fluxes, sign_y_corner_fluxes, directional_lagrangian_volume, directional_lagrangian_density);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));
  time_compute_volume_fluxes_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientX(nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_density,
		       density_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientY(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_density,
		       density_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));   

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_density,
			       density_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_density,
			       density_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2XCornerFluxes(nx,  ny,  halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_x, volume_fluxes_x, directional_lagrangian_density, density_gradient,
					 reconstructed_density_faces_x, mass_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2YCornerFluxes(nx,  ny,  halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, volume_fluxes_y, directional_lagrangian_density, density_gradient_y,
					 reconstructed_density_faces_y, mass_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2XCornerFluxesBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_x, volume_fluxes_x, directional_lagrangian_density, density_gradient,
						 reconstructed_density_faces_x, mass_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2YCornerFluxesBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, volume_fluxes_y, directional_lagrangian_density, density_gradient_y,
						 reconstructed_density_faces_y, mass_flux_y);     
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2Corner(nx, ny, dx, dy, dt, directional_lagrangian_density, reconstructed_density_faces_x, reconstructed_density_faces_y, volume_fluxes_corner, predicted_u, predicted_v, 
				  mass_corner_fluxes);
  //ReconstructMassFluxOrder2DiagCorner( nx, ny, dx, dy, dt, predicted_u, predicted_v, directional_lagrangian_density, volume_fluxes_corner, density_gradient_diag, density_gradient_antidiag,
  //				       mass_corner_fluxes);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
  time_project_mass_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2CornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, directional_lagrangian_density, reconstructed_density_faces_x, reconstructed_density_faces_y, volume_fluxes_corner, predicted_u, predicted_v, 
					  mass_corner_fluxes);     
  //ReconstructMassFluxOrder2DiagCornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, predicted_u, predicted_v, directional_lagrangian_density, volume_fluxes_corner, density_gradient_diag, density_gradient_antidiag,
  //					      mass_corner_fluxes);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassDirect(nx, ny, in_cell_mass, mass_flux_x, mass_flux_y, out_cell_mass);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
  time_project_mass_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassCornerFluxes(nx, ny, dt, mass_corner_fluxes, predicted_u, predicted_v, 
			  out_cell_mass);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
  time_project_mass_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassCornerFluxesBoundary(BoundaryConditions, nx, ny, dt, mass_corner_fluxes, predicted_u, predicted_v, 
				  out_cell_mass);     
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));


  // ========================
  //     Projection _ e. (projection of mass*energy then back to e )
  // ========================
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientX(nx, ny, dx, dy, volume_fluxes_x, e_lag,
		       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientY(nx, ny, dx, dy, volume_fluxes_y, e_lag,
		       energy_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, e_lag,
			       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, e_lag,
			       energy_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XCornerFluxes(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_x, mass_flux_x, e_lag, energy_gradient,
						      reconstructed_energy_faces_x, energy_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YCornerFluxes(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, mass_flux_y, e_lag, energy_gradient_y,
						      reconstructed_energy_faces_y, energy_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_x, mass_flux_x, e_lag, energy_gradient,
							      reconstructed_energy_faces_x, energy_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, mass_flux_y, e_lag, energy_gradient_y,
							      reconstructed_energy_faces_y, energy_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  //ReconstructIntensiveVariableFluxOrder2Corner(nx, ny, dx, dy, dt, e_lag,reconstructed_energy_faces_x, reconstructed_energy_faces_y, mass_corner_fluxes, predicted_u, predicted_v,
  //					       energy_flux_corner);
  ReconstructIntensiveVariableFluxOrder2DiagCorner(nx, ny, dx, dy, dt, predicted_u, predicted_v, e_lag, mass_corner_fluxes, energy_gradient_diag, energy_gradient_antidiag,
  						     energy_flux_corner);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  //ReconstructIntensiveVariableFluxOrder2CornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, e_lag, reconstructed_energy_faces_x, reconstructed_energy_faces_y, mass_corner_fluxes, predicted_u, predicted_v,
  //						       energy_flux_corner);
  ReconstructIntensiveVariableFluxOrder2DiagCornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, predicted_u, predicted_v, e_lag, mass_corner_fluxes, energy_gradient_diag, energy_gradient_antidiag,
  							     energy_flux_corner);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableDirect(nx, ny, in_cell_mass, e_lag, energy_flux_x, energy_flux_y, out_cell_mass,
				     out_e);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableCornerFluxes(nx, ny, dt, e_lag, energy_flux_corner, predicted_u, predicted_v, out_cell_mass,
					   out_e);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableCornerFluxesBoundary(BoundaryConditions, nx, ny, dt, e_lag, energy_flux_corner, predicted_u, predicted_v, out_cell_mass,
						   out_e);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

 
  // ========================
  //     Projection _ concentration c_1. (projection of mass*c_1 then back to c_1 )
  // ========================
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientX(nx, ny, dx, dy, volume_fluxes_x, in_c_1,
		       concentration_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientY(nx, ny, dx, dy, volume_fluxes_y, in_c_1,
		       concentration_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, in_c_1,
			       concentration_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, in_c_1,
			       concentration_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XCornerFluxes(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_x, mass_flux_x, in_c_1, concentration_gradient,
						      reconstructed_concentration_faces_x, concentration_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YCornerFluxes(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, mass_flux_y, in_c_1, concentration_gradient_y,
						      reconstructed_concentration_faces_y, concentration_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy,  dt, predicted_u, predicted_v, volume_fluxes_x, mass_flux_x, in_c_1, concentration_gradient,
							      reconstructed_concentration_faces_x, concentration_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_y, mass_flux_y, in_c_1, concentration_gradient_y,
							      reconstructed_concentration_faces_y, concentration_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  //ReconstructIntensiveVariableFluxOrder2Corner(nx, ny, dx, dy, dt, in_c_1, reconstructed_concentration_faces_x, reconstructed_concentration_faces_y, mass_corner_fluxes, predicted_u, predicted_v,
  //					       concentration_flux_corner);
  ReconstructIntensiveVariableFluxOrder2DiagCorner(nx, ny, dx, dy, dt, predicted_u, predicted_v, in_c_1, mass_corner_fluxes, concentration_gradient_diag, concentration_gradient_antidiag,
  						     concentration_flux_corner);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  //ReconstructIntensiveVariableFluxOrder2CornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, in_c_1, reconstructed_concentration_faces_x, reconstructed_concentration_faces_y, mass_corner_fluxes, predicted_u, predicted_v,
  //						       concentration_flux_corner);
  ReconstructIntensiveVariableFluxOrder2DiagCornerBoundary(BoundaryConditions, nx, ny, dx, dy, dt, predicted_u, predicted_v, in_c_1, mass_corner_fluxes, concentration_gradient_diag, concentration_gradient_antidiag,
  						     concentration_flux_corner);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableDirect(nx, ny, in_cell_mass, in_c_1, concentration_flux_x, concentration_flux_y, out_cell_mass,
				     out_c_1);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableCornerFluxes(nx, ny, dt, in_c_1, concentration_flux_corner, predicted_u, predicted_v, out_cell_mass,
					   out_c_1);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableCornerFluxesBoundary(BoundaryConditions, nx, ny, dt, in_c_1, concentration_flux_corner, predicted_u, predicted_v, out_cell_mass,
						   out_c_1);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  
  // ========================
  //     Projection _ u.   (projection of mu, then back to u)   
  // ========================
        
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalX( nx, ny, dx, dt, predicted_u,  u_lag,
			     gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalY( nx, ny, dy, dt, predicted_v,  u_lag,
			     gradient_u_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalXBoundary(BoundaryConditions, nx, ny, dx, dt, predicted_u,  u_lag, 
				    gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalYBoundary(BoundaryConditions,  nx, ny, dy, dt, predicted_v, u_lag, 
				    gradient_u_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2Direct(nx,  ny,  halo_width, dx, dy, dt, in_cell_mass, out_cell_mass, predicted_u, predicted_v, u_lag, gradient_u, gradient_u_y, mass_flux_x, mass_flux_y,
					       out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2BoundaryDirect(BoundaryConditions, "project_ux", nx,  ny,  halo_width, dx, dy, dt, in_cell_mass, out_cell_mass, predicted_u, predicted_v, u_lag, gradient_u, gradient_u_y, mass_flux_x, mass_flux_y,
					       out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableCornerFluxes(nx,  ny, dt,  halo_width, out_cell_mass, u_lag, u_lag, v_lag, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes,
					    out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableCornerFluxesBoundary(BoundaryConditions, nx,  ny, dt, halo_width, out_cell_mass, u_lag, u_lag, v_lag, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes,
  						    out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));


  // ========================
  //     Projection _ v.   (projection of mv, then back to v)   
  // ========================
        
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalX( nx, ny, dx, dt, predicted_u,  v_lag,
			     gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalY( nx, ny, dy, dt, predicted_v,  v_lag,
			     gradient_v_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_Y.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalXBoundary(BoundaryConditions, nx, ny, dx, dt, predicted_u,  v_lag, 
				    gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalYBoundary(BoundaryConditions,  nx, ny, dy, dt, predicted_v, v_lag, 
				    gradient_v_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2Direct(nx,  ny,  halo_width, dx, dy, dt, in_cell_mass, out_cell_mass, predicted_u, predicted_v, v_lag, gradient_v, gradient_v_y, mass_flux_x, mass_flux_y,
					       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2BoundaryDirect(BoundaryConditions, "project_uy", nx,  ny,  halo_width, dx, dy, dt, in_cell_mass, out_cell_mass, predicted_u, predicted_v, v_lag, gradient_v, gradient_v_y, mass_flux_x, mass_flux_y,
					       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableCornerFluxes(nx,  ny, dt,  halo_width, out_cell_mass, v_lag, u_lag, v_lag, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes,
					    out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableCornerFluxesBoundary(BoundaryConditions, nx,  ny, dt,  halo_width, out_cell_mass, v_lag, u_lag, v_lag, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes,
						    out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));


  // ========================
  //     Reconstruction of each material quantities.   
  // ========================

  #include "reconstruct_material_quantities_mix.h"

}



// -*- c++ -*- (for emacs users)
#include "direct_projection_2d_multimat_interface_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"
#include "cartesian_connectivity.h"


void DirectProjection2dMultimaterialInterfaceDriver(//in
						    const std::string BoundaryConditions,
						    const int nx,
						    const int ny,
						    const int nb_faces_x,
						    const int nb_faces_y,
						    const int nb_cells,
						    const int nb_nodes,
						    const RealType dx,
						    const RealType dy,
						    const RealType dt,
						    const RealType halo_width,
						    const RealType gamma_1,
						    const RealType gamma_2,
						    const RealType* RESTRICT predicted_u,
						    const RealType* RESTRICT predicted_v,
						    const RealType* RESTRICT e_lag,
						    const RealType* RESTRICT e_1_lag,
						    const RealType* RESTRICT e_2_lag,
						    const RealType* RESTRICT u_lag,
						    const RealType* RESTRICT v_lag,
						    const RealType* RESTRICT in_rho_1,
						    const RealType* RESTRICT in_rho_2,
						    const RealType* RESTRICT in_cell_mass,
						    const RealType* RESTRICT in_cell_mass_1,
						    const RealType* RESTRICT in_cell_mass_2,
						    const RealType* RESTRICT in_cell_volumic_fraction,
						    const RealType* RESTRICT cell_volumes,
						    // out
						    RealType* RESTRICT out_u,
						    RealType* RESTRICT out_v,
						    RealType* RESTRICT out_e,
						    RealType* RESTRICT out_cell_mass,
						    RealType* RESTRICT out_cell_mass_1,
						    RealType* RESTRICT out_cell_mass_2,
						    RealType* RESTRICT out_e_1,
						    RealType* RESTRICT out_e_2,
						    RealType* RESTRICT out_rho,
						    RealType* RESTRICT out_rho_1,
						    RealType* RESTRICT out_rho_2,
						    RealType* RESTRICT out_cell_volumic_fraction,
						    RealType* RESTRICT directional_lagrangian_volume,
						    RealType* RESTRICT directional_lagrangian_volume_y,
						    RealType* RESTRICT directional_lagrangian_density,
						    RealType* RESTRICT directional_lagrangian_density_y,
						    RealType* RESTRICT directional_lagrangian_density_1,
						    RealType* RESTRICT directional_lagrangian_density_1_y,
						    RealType* RESTRICT directional_lagrangian_density_2,
						    RealType* RESTRICT directional_lagrangian_density_2_y,
						    RealType* RESTRICT volume_fluxes_x,
						    RealType* RESTRICT volume_fluxes_y,
						    RealType* RESTRICT volume_fluxes_1_x,
						    RealType* RESTRICT volume_fluxes_1_y,
						    RealType* RESTRICT volume_fluxes_2_x,
						    RealType* RESTRICT volume_fluxes_2_y,
						    RealType* RESTRICT mass_flux_x,
						    RealType* RESTRICT mass_flux_y,
						    RealType* RESTRICT mass_flux_1_x,
						    RealType* RESTRICT mass_flux_1_y,
						    RealType* RESTRICT mass_flux_2_x,
						    RealType* RESTRICT mass_flux_2_y,
						    RealType* RESTRICT energy_flux_x,
						    RealType* RESTRICT energy_flux_y,
						    RealType* RESTRICT energy_flux_1_x,
						    RealType* RESTRICT energy_flux_1_y,
						    RealType* RESTRICT energy_flux_2_x,
						    RealType* RESTRICT energy_flux_2_y,
						    RealType* RESTRICT bool_check_fluxes_x,
						    RealType* RESTRICT bool_check_fluxes_y,
						    RealType* RESTRICT density_gradient,
						    RealType* RESTRICT density_gradient_y,
						    RealType* RESTRICT density_1_gradient,
						    RealType* RESTRICT density_1_gradient_y,
						    RealType* RESTRICT density_2_gradient,
						    RealType* RESTRICT density_2_gradient_y,
						    RealType* RESTRICT energy_gradient,
						    RealType* RESTRICT energy_gradient_y,
						    RealType* RESTRICT energy_1_gradient,
						    RealType* RESTRICT energy_1_gradient_y,
						    RealType* RESTRICT energy_2_gradient,
						    RealType* RESTRICT energy_2_gradient_y,
						    RealType* RESTRICT gradient_u,
						    RealType* RESTRICT gradient_u_y,
						    RealType* RESTRICT gradient_v,
						    RealType* RESTRICT gradient_v_y,
						    RealType* RESTRICT interface_normal_x,
						    RealType* RESTRICT interface_normal_y,
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

  //========================
  //    Projection mass .
  //========================
 
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructNormalInterface(nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructNormalInterfaceBoundary(BoundaryConditions, nx, ny, dx, dy, in_cell_volumic_fraction, interface_normal_x, interface_normal_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianQuantitiesMultimatX(nx, ny, dt, dx, dy, predicted_u, in_cell_mass_1, in_cell_mass_2, in_cell_volumic_fraction,
						  volume_fluxes_x, directional_lagrangian_volume, directional_lagrangian_density_1, directional_lagrangian_density_2);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianQuantitiesMultimatY(nx, ny, dt, dx, dy, predicted_v, in_cell_mass_1, in_cell_mass_2, in_cell_volumic_fraction,
						  volume_fluxes_y, directional_lagrangian_volume_y, directional_lagrangian_density_1_y, directional_lagrangian_density_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianFractionalVolumeFluxesX(nx, ny, dx, dy, volume_fluxes_x,directional_lagrangian_volume, in_cell_volumic_fraction, interface_normal_x, interface_normal_y,
						  volume_fluxes_1_x, volume_fluxes_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianFractionalVolumeFluxesY(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_volume_y, in_cell_volumic_fraction, interface_normal_x, interface_normal_y,
						  volume_fluxes_1_y, volume_fluxes_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianFractionalVolumeFluxesXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_volume, in_cell_volumic_fraction, interface_normal_x, interface_normal_y,
						  volume_fluxes_1_x, volume_fluxes_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
   
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianFractionalVolumeFluxesYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_volume_y, in_cell_volumic_fraction, interface_normal_x, interface_normal_y,
						  volume_fluxes_1_y, volume_fluxes_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  CheckDirectionalVolumeFluxesXY(nx, ny, dx, dy, volume_fluxes_x, volume_fluxes_y, in_cell_volumic_fraction,
				 bool_check_fluxes_x, bool_check_fluxes_y, volume_fluxes_1_x, volume_fluxes_1_y, volume_fluxes_2_x, volume_fluxes_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));
  time_compute_volume_fluxes_Y.push_back(diff(time_begin, time_end)); 

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatX(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_2_x, directional_lagrangian_density_1, directional_lagrangian_density_2,
			       density_1_gradient, density_2_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatY(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_y, volume_fluxes_2_y, directional_lagrangian_density_1_y, directional_lagrangian_density_2_y,
			       density_1_gradient_y, density_2_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatXBoundary(BoundaryConditions, nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_2_x, directional_lagrangian_density_1, directional_lagrangian_density_2,
			       density_1_gradient, density_2_gradient); 
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatYBoundary(BoundaryConditions, nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_y, volume_fluxes_2_y, directional_lagrangian_density_1_y, directional_lagrangian_density_2_y,
			       density_1_gradient_y, density_2_gradient_y); 
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2X(nx,  ny,  halo_width, dx, dy, volume_fluxes_1_x,
			     directional_lagrangian_density_1, density_1_gradient,
			     mass_flux_1_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_X.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2Y(nx,  ny,  halo_width, dx, dy, volume_fluxes_1_y,
			     directional_lagrangian_density_1_y, density_1_gradient_y,
			     mass_flux_1_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2X(nx,  ny,  halo_width, dx, dy, volume_fluxes_2_x,
			     directional_lagrangian_density_2, density_2_gradient,
			     mass_flux_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2Y(nx,  ny,  halo_width, dx, dy, volume_fluxes_2_y,
			     directional_lagrangian_density_2_y, density_2_gradient_y,
			     mass_flux_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2XBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_1_x,
				     directional_lagrangian_density_1, density_1_gradient,
				     mass_flux_1_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2YBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_1_y,
				     directional_lagrangian_density_1_y, density_1_gradient_y,
				     mass_flux_1_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2XBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_2_x,
				     directional_lagrangian_density_2, density_2_gradient,
				     mass_flux_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end)); 
   
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2YBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_2_y,
				     directional_lagrangian_density_2_y, density_2_gradient_y,
				     mass_flux_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  CheckDirectionalMassFluxXY(nx, ny, dx, dy, in_cell_mass_1, in_cell_mass_2, volume_fluxes_x, volume_fluxes_y, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_1_y, volume_fluxes_2_x, volume_fluxes_2_y,
			     mass_flux_1_x, mass_flux_1_y, mass_flux_2_x, mass_flux_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_X.push_back(diff(time_begin, time_end));
  time_mass_reconstruct_o2_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassDirect(nx, ny, in_cell_mass_1, mass_flux_1_x, mass_flux_1_y, out_cell_mass_1);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
  time_project_mass_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassDirect(nx, ny, in_cell_mass_2, mass_flux_2_x, mass_flux_2_y, out_cell_mass_2);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
  time_project_mass_Y.push_back(diff(time_begin, time_end));

#include "reconstruct_mass_flux_X_interface.h"
#include "reconstruct_mass_flux_Y_interface.h"

  ProjectDensityMultimatDirect(nx, ny, dx, dy, volume_fluxes_x, volume_fluxes_y, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_1_y, volume_fluxes_2_x, volume_fluxes_2_y, in_rho_1, in_rho_2, in_cell_mass, mass_flux_x, mass_flux_y, out_cell_mass_1, out_cell_mass_2,
			       out_cell_volumic_fraction, out_cell_mass, out_rho_1, out_rho_2);

	      
  // ========================
  //     Projection e. (projection of mass*energy then back to e )
  // ========================
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatX(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_2_x, e_1_lag, e_2_lag,
			       energy_1_gradient, energy_2_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatY(nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_y, volume_fluxes_2_y, e_1_lag, e_2_lag,
			       energy_1_gradient_y, energy_2_gradient_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatXBoundary(BoundaryConditions, nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_x, volume_fluxes_2_x, e_1_lag, e_2_lag,
				       energy_1_gradient, energy_2_gradient); 
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientMultimatYBoundary(BoundaryConditions, nx, ny, dx, dy, in_cell_volumic_fraction, volume_fluxes_1_y, volume_fluxes_2_y, e_1_lag, e_2_lag,
				       energy_1_gradient_y, energy_2_gradient_y); 
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2X(nx, ny, halo_width, dx, dy, volume_fluxes_1_x, mass_flux_1_x, e_1_lag, energy_1_gradient,
					  energy_flux_1_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2Y(nx, ny, halo_width, dx, dy, volume_fluxes_1_y, mass_flux_1_y, e_1_lag, energy_1_gradient_y,
					  energy_flux_1_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2X(nx, ny, halo_width, dx, dy, volume_fluxes_2_x, mass_flux_2_x, e_2_lag, energy_2_gradient,
					  energy_flux_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2Y(nx, ny, halo_width, dx, dy, volume_fluxes_2_y, mass_flux_2_y, e_2_lag, energy_2_gradient_y,
					  energy_flux_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_1_x, mass_flux_1_x, e_1_lag, energy_1_gradient,
						  energy_flux_1_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_1_y, mass_flux_1_y, e_1_lag, energy_1_gradient_y,
						  energy_flux_1_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_2_x, mass_flux_2_x, e_2_lag, energy_2_gradient,
						  energy_flux_2_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_2_y, mass_flux_2_y, e_2_lag, energy_2_gradient_y,
						  energy_flux_2_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableDirect(nx, ny, in_cell_mass_1, e_1_lag, energy_flux_1_x, energy_flux_1_y, out_cell_mass_1,
				     out_e_1);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableDirect(nx, ny, in_cell_mass_2, e_2_lag, energy_flux_2_x, energy_flux_2_y, out_cell_mass_2,
				     out_e_2);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  time_project_energy_Y.push_back(diff(time_begin, time_end));

    
  // ========================
  //     Projection u.   (projection of mu, then back to u)   
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


  // ========================
  //     Projection v.   (projection of mv, then back to v)   
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


  // ========================
  //   Reconstruction of thermodynamics quantities
  // ========================

#include "reconstruct_material_quantities_interface.h"
 
}


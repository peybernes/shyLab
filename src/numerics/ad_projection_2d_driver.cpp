// -*- c++ -*- (for emacs users)
#include "ad_projection_2d_driver.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"

void AdProjection2dXDriver(const std::string BoundaryConditions,
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
			   const RealType* RESTRICT predicted_u,
			   const RealType* RESTRICT e_lag,
			   const RealType* RESTRICT u_lag,
			   const RealType* RESTRICT v_lag,
			   const RealType* RESTRICT in_cell_mass,
			   RealType* RESTRICT out_u,
			   RealType* RESTRICT out_v,
			   RealType* RESTRICT out_e,
			   RealType* RESTRICT out_cell_mass,
			   RealType* RESTRICT directional_lagrangian_volume,
			   RealType* RESTRICT directional_lagrangian_density,
			   RealType* RESTRICT volume_fluxes_x,
			   RealType* RESTRICT mass_flux_x,
			   RealType* RESTRICT energy_flux_x,
			   RealType* RESTRICT density_gradient,
			   RealType* RESTRICT energy_gradient,
			   RealType* RESTRICT gradient_u,
			   RealType* RESTRICT gradient_v,
			   //timing
			   std::vector<RealType>& time_compute_volume_fluxes_X,
			   std::vector<RealType>& time_gradient_X,
			   std::vector<RealType>& time_mass_reconstruct_o2_X,
			   std::vector<RealType>& time_project_mass_X,
			   std::vector<RealType>& time_reconstruct_energy_o2_X,
			   std::vector<RealType>& time_project_energy_X,
			   std::vector<RealType>& time_gradient_nodal_X,
			   std::vector<RealType>& time_project_nodal_velocity_X,
			   std::vector<RealType> time_periodic_boundary)  {
  
  //for timing
  struct timespec time_begin;
  struct timespec time_end;

  //========================
  //    Projection X _mass .
  //========================
        
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianQuantitiesX(nx, ny, dt, dx, dy, predicted_u, in_cell_mass,
					  volume_fluxes_x, directional_lagrangian_volume, directional_lagrangian_density);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientX(nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_density,
		       density_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, directional_lagrangian_density,
			       density_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2X(nx,  ny,  halo_width, dx, dy, volume_fluxes_x,
			     directional_lagrangian_density,density_gradient,
			     mass_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2XBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_x,
				     directional_lagrangian_density,density_gradient,
				     mass_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassX(nx, ny, in_cell_mass, mass_flux_x, out_cell_mass);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_X.push_back(diff(time_begin, time_end));
    
           
  // ========================
  //     Projection X _ e. (projection of mass*energy then back to e )
  // ========================
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientX(nx, ny, dx, dy, volume_fluxes_x, e_lag,
		       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientXBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_x, e_lag,
			       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2X(nx, ny, halo_width, dx, dy, volume_fluxes_x, mass_flux_x, e_lag, energy_gradient,
					  energy_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2XBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_x, mass_flux_x, e_lag, energy_gradient,
						  energy_flux_x);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableX(nx, ny , in_cell_mass, e_lag, energy_flux_x, out_cell_mass,
				out_e);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_X.push_back(diff(time_begin, time_end));
  
  
     
  // ========================
  //     Projection X _ u.   (projection of mu, then back to u)   
  // ========================
  

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalX( nx, ny, dx, dt, predicted_u,  u_lag,
			     gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_X.push_back(diff(time_begin, time_end));
  
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalXBoundary(BoundaryConditions, nx, ny, dx, dt, predicted_u,  u_lag, 
				    gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2X(nx, ny,  halo_width, dx, dt, in_cell_mass, out_cell_mass, predicted_u, u_lag, gradient_u, mass_flux_x,
				       out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2XBoundary(BoundaryConditions, "project_ux", nx,  ny,  halo_width, dx, dt, in_cell_mass, out_cell_mass, predicted_u, u_lag, gradient_u, mass_flux_x,
					       out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

  // ========================
  //     Projection X _ v.   (projection of mv, then back to v)   
  // ========================

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalX( nx, ny, dx, dt, predicted_u,  v_lag,
			     gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_X.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalXBoundary(BoundaryConditions, nx, ny, dx, dt, predicted_u,  v_lag,
				    gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
    
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2X(nx, ny,  halo_width, dx, dt, in_cell_mass, out_cell_mass, predicted_u, v_lag, gradient_v, mass_flux_x,
				       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_X.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2XBoundary(BoundaryConditions, "project_uy", nx,  ny,  halo_width, dx, dt, in_cell_mass, out_cell_mass, predicted_u, v_lag, gradient_v, mass_flux_x,
					       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));

          
}

void AdProjection2dYDriver(const std::string BoundaryConditions,
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
			   const RealType* RESTRICT predicted_v,
			   const RealType* RESTRICT e_lag,
			   const RealType* RESTRICT u_lag,
			   const RealType* RESTRICT v_lag,
			   const RealType* RESTRICT in_cell_mass,
			   RealType* RESTRICT out_u,
			   RealType* RESTRICT out_v,
			   RealType* RESTRICT out_e,
			   RealType* RESTRICT out_cell_mass,
			   RealType* RESTRICT directional_lagrangian_volume_y,
			   RealType* RESTRICT directional_lagrangian_density_y,
			   RealType* RESTRICT volume_fluxes_y,
			   RealType* RESTRICT mass_flux_y,
			   RealType* RESTRICT energy_flux_y,
			   RealType* RESTRICT density_gradient,
			   RealType* RESTRICT energy_gradient,
			   RealType* RESTRICT gradient_u,
			   RealType* RESTRICT gradient_v,
			   //timing
			   std::vector<RealType>& time_compute_volume_fluxes_Y,
			   std::vector<RealType>& time_gradient_Y,
			   std::vector<RealType>& time_mass_reconstruct_o2_Y,
			   std::vector<RealType>& time_project_mass_Y,
			   std::vector<RealType>& time_reconstruct_energy_o2_Y,
			   std::vector<RealType>& time_project_energy_Y,
			   std::vector<RealType>& time_gradient_nodal_Y,
			   std::vector<RealType>& time_project_nodal_velocity_Y,		     
			   std::vector<RealType> time_periodic_boundary)  {
  
  //for timing
  struct timespec time_begin;
  struct timespec time_end;
  //

     
  //===================
  // //Projection Y.
  //===================
     

  //==============
  // mass Y
  //==============

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ComputeDirectionalLagrangianQuantitiesY(nx, ny, dt, dx, dy, predicted_v, in_cell_mass,
					  volume_fluxes_y, directional_lagrangian_volume_y, directional_lagrangian_density_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_compute_volume_fluxes_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientY(nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_density_y,
		       density_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, directional_lagrangian_density_y,
			       density_gradient);

  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2Y(nx,  ny,  halo_width, dx, dy, volume_fluxes_y, directional_lagrangian_density_y, density_gradient,
			     mass_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_mass_reconstruct_o2_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructMassFluxOrder2YBoundary(BoundaryConditions, nx,  ny,  halo_width, dx, dy, volume_fluxes_y, directional_lagrangian_density_y,density_gradient,
				     mass_flux_y);
     
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectMassY(nx, ny, in_cell_mass, mass_flux_y, out_cell_mass);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_mass_Y.push_back(diff(time_begin, time_end));
     
     
  //==============
  //energy Y
  //==============
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientY(nx, ny, dx, dy, volume_fluxes_y, e_lag,
		       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientYBoundary(BoundaryConditions, nx, ny, dx, dy, volume_fluxes_y, e_lag,
			       energy_gradient);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2Y(nx, ny, halo_width, dx, dy, volume_fluxes_y, mass_flux_y, e_lag, energy_gradient,
					  energy_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_reconstruct_energy_o2_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructIntensiveVariableFluxOrder2YBoundary(BoundaryConditions, nx, ny, halo_width, dx, dy, volume_fluxes_y, mass_flux_y, e_lag, energy_gradient,
						  energy_flux_y);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  MassProjectIntensiveVariableY(nx, ny , in_cell_mass,e_lag, energy_flux_y, out_cell_mass,
				out_e);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_energy_Y.push_back(diff(time_begin, time_end));
        
  //==============
  // velocity Y   (u)
  //==============
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalY( nx, ny, dy, dt, predicted_v,  u_lag,
			     gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalYBoundary(BoundaryConditions,  nx, ny, dy, dt, predicted_v, u_lag, 
				    gradient_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2Y(nx, ny, halo_width, dy, dt, in_cell_mass, out_cell_mass, predicted_v, u_lag, gradient_u, mass_flux_y,
				       out_u);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2YBoundary(BoundaryConditions, "project_ux", nx, ny, halo_width, dy, dt, in_cell_mass, out_cell_mass, predicted_v, u_lag, gradient_u, mass_flux_y,
					       out_u);
     
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
     
     
  //==============
  // velocity Y   (v)
  //==============
     
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalY( nx, ny, dy, dt, predicted_v,  v_lag,
			     gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_gradient_nodal_Y.push_back(diff(time_begin, time_end));

  clock_gettime(CLOCK_REALTIME, &time_begin);
  ReconstructGradientNodalYBoundary(BoundaryConditions,  nx, ny, dy, dt, predicted_v,  v_lag, 
				    gradient_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2Y(nx,  ny,  halo_width, dy, dt, in_cell_mass, out_cell_mass, predicted_v, v_lag, gradient_v, mass_flux_y,
				       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_project_nodal_velocity_Y.push_back(diff(time_begin, time_end));
     
  clock_gettime(CLOCK_REALTIME, &time_begin);
  ProjectNodalIntensiveVariableOrder2YBoundary(BoundaryConditions, "project_uy", nx,  ny,  halo_width, dx, dt, in_cell_mass, out_cell_mass, predicted_v, v_lag, gradient_v, mass_flux_y,
					       out_v);
  clock_gettime(CLOCK_REALTIME, &time_end);
  time_periodic_boundary.push_back(diff(time_begin, time_end));
     
}

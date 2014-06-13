// -*- c++ -*- (for emacs users)
#include "direct_projection_2d_interface.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <malloc.h>

#include "kernel_tools.h"

void DirectProjection_2d(const std::string BoundaryConditions, 
			 const int nx,
		     const int ny,
		     const RealType dx,
		     const RealType dy,
		     const RealType dt,
		     const RealType halo_width,
		     const RealType* RESTRICT predicted_u,
		     const RealType* RESTRICT predicted_v,
		     RealType* RESTRICT e_lag,
		     RealType* RESTRICT u_lag,
		     RealType* RESTRICT v_lag,
		     RealType* RESTRICT in_u,
		     RealType* RESTRICT in_v,
		     RealType* RESTRICT in_e,
		     RealType* RESTRICT in_cell_mass,
		     RealType* RESTRICT out_u,
		     RealType* RESTRICT out_v,
		     RealType* RESTRICT out_e,
		     RealType* RESTRICT out_cell_mass,
		     RealType* RESTRICT directional_lagrangian_volume,
		     RealType* RESTRICT directional_lagrangian_density,
		     RealType* RESTRICT directional_lagrangian_volume_y,
		     RealType* RESTRICT directional_lagrangian_density_y,
		     RealType* RESTRICT volume_fluxes_x,
		     RealType* RESTRICT volume_fluxes_y,
		     RealType* RESTRICT mass_flux_x,
		     RealType* RESTRICT mass_flux_y,
		     RealType* RESTRICT energy_flux_x,
		     RealType* RESTRICT energy_flux_y,
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
		     std::vector<RealType>& time_periodic_boundary
		     ) {


    //for timing
      struct timespec time1;
      struct timespec time2;
    //



    //========================
    //    Projection _ mass .
    //========================
    
      
     clock_gettime(CLOCK_REALTIME, &time1);
     ComputeDirectionalLagrangianQuantitiesX(nx, ny, dt, dx, dy, predicted_u,
					     in_cell_mass,
					     volume_fluxes_x, directional_lagrangian_volume, directional_lagrangian_density);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_compute_volume_fluxes_0.push_back(diff(time1, time2));
    
     //     CheckFluxPeriodicalPropertyX(nx, ny, volume_fluxes_x);

     clock_gettime(CLOCK_REALTIME, &time1);
     ComputeDirectionalLagrangianQuantitiesY(nx, ny, dt, dx, dy, predicted_v,
					     in_cell_mass,
					     volume_fluxes_y, directional_lagrangian_volume_y, directional_lagrangian_density_y);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_compute_volume_fluxes_0.push_back(diff(time1, time2));
    
     //     CheckFluxPeriodicalPropertyY(nx, ny, volume_fluxes_y);

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructMassFluxOrder1X(nx, ny, halo_width, volume_fluxes_x, directional_lagrangian_density, mass_flux_x);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_reconstruct_0.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructMassFluxOrder1PeriodicBoundaryX(nx, ny, halo_width, volume_fluxes_x, directional_lagrangian_density, mass_flux_x);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_PeriodicBoundary.push_back(diff(time1, time2));

     // CheckFluxPeriodicalPropertyX(nx, ny, mass_flux_x);

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructMassFluxOrder1Y(nx, ny, halo_width, volume_fluxes_y, directional_lagrangian_density_y, mass_flux_y);
     clock_gettime(CLOCK_REALTIME, &time2);
     //time_reconstruct_0.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructMassFluxOrder1PeriodicBoundaryY(nx, ny, halo_width, volume_fluxes_y, directional_lagrangian_density_y, mass_flux_y);
     clock_gettime(CLOCK_REALTIME, &time2);
     //time_PeriodicBoundary.push_back(diff(time1, time2));

     // CheckFluxPeriodicalPropertyY(nx, ny, mass_flux_y);

     clock_gettime(CLOCK_REALTIME, &time1);
     ProjectMassDirect(nx, ny, in_cell_mass, mass_flux_x, mass_flux_y, out_cell_mass);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_project_mass_0.push_back(diff(time1, time2));


     // ========================
     //     Projection _ e. (projection of mass*energy then back to e )
     // ========================
     
     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructIntensiveVariableFluxOrder1X(nx, ny, halo_width, mass_flux_x, e_lag,
					    energy_flux_x);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_reconstruct_0.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructIntensiveVariableFluxOrder1PeriodicBoundaryX(nx, ny, halo_width, mass_flux_x, e_lag,
						     energy_flux_x );
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_PeriodicBoundary.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructIntensiveVariableFluxOrder1Y(nx, ny, halo_width, mass_flux_y, e_lag,
					    energy_flux_y);
     clock_gettime(CLOCK_REALTIME, &time2);
     //time_reconstruct_0.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ReconstructIntensiveVariableFluxOrder1PeriodicBoundaryY(nx, ny, halo_width, mass_flux_y, e_lag,
						    energy_flux_y);
     clock_gettime(CLOCK_REALTIME, &time2);
     //time_PeriodicBoundary.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     MassProjectIntensiveVariableDirect(nx, ny, in_cell_mass, e_lag, energy_flux_x, energy_flux_y, out_cell_mass,
				   out_e);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_project_intensive_variable_0.push_back(diff(time1, time2));

     
     // ========================
     //     Projection _ u,v.   (projection of mu,mv then back to u,v)   
     // ========================
        
     clock_gettime(CLOCK_REALTIME, &time1);
     ProjectNodalIntensiveVariableDirect(nx, ny, halo_width, in_cell_mass, out_cell_mass, u_lag, v_lag, mass_flux_x, mass_flux_y,
					 out_u, out_v);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_project_nodal_intensive_variable_0.push_back(diff(time1, time2));

     clock_gettime(CLOCK_REALTIME, &time1);
     ProjectNodalIntensiveVariablePeriodicBoundaryDirect(nx, ny, halo_width, in_cell_mass, out_cell_mass, u_lag, v_lag, mass_flux_x, mass_flux_y,
						 out_u, out_v);
     clock_gettime(CLOCK_REALTIME, &time2);
     // time_PeriodicBoundary.push_back(diff(time1, time2));




}

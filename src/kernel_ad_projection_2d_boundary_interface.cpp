// -*- c++ -*- (for emacs users)
#include "kernel_ad_projection_2d_boundary_interface.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"

void BoundaryVelocityPrediction(std::string BoundaryConditions,
					int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					const RealType* RESTRICT in_mass,
					const RealType* RESTRICT in_pressure,
					const RealType* RESTRICT in_pseudo_pressure,
					const RealType* RESTRICT in_velocity_x,
					const RealType* RESTRICT in_velocity_y,
					RealType* RESTRICT out_velocity_x,
					RealType* RESTRICT out_velocity_y)
{
  if (BoundaryConditions == "periodic") {
    PeriodicBoundaryVelocityPrediction(nx, ny, dt, dx, dy, in_mass, in_pressure, in_pseudo_pressure, in_velocity_x, in_velocity_y, out_velocity_x, out_velocity_y);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    RtBoundaryVelocityPrediction(nx, ny, dt, dx, dy, in_mass, in_pressure, in_pseudo_pressure, in_velocity_x, in_velocity_y, out_velocity_x, out_velocity_y);
  }
  else{
    PeriodicBoundaryVelocityPrediction(nx, ny, dt, dx, dy, in_mass, in_pressure, in_pseudo_pressure, in_velocity_x, in_velocity_y, out_velocity_x, out_velocity_y);
  }
  
}

void BoundaryCopy(std::string BoundaryConditions,
			  int nx,
			  int ny,
			  RealType* RESTRICT in_velocity_x,
			  RealType* RESTRICT in_velocity_y)
{
  if (BoundaryConditions == "periodic") {
    PeriodicBoundaryCopy(nx, ny, in_velocity_x, in_velocity_y);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    RtBoundaryCopy(nx, ny, in_velocity_x, in_velocity_y);
  }
  else{
    PeriodicBoundaryCopy(nx, ny, in_velocity_x, in_velocity_y);
  }

}

void ReconstructGradientXBoundary(std::string BoundaryConditions,
				  index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT volume_fluxes,
				  const RealType* RESTRICT lag_variable,
				  RealType* RESTRICT gradient_variable) {

  if (BoundaryConditions == "periodic") {
    ReconstructGradientXPeriodicBoundary(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructGradientXBoundaryRt(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable);
  }
  else{
    ReconstructGradientXPeriodicBoundary(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable);
  }
}

  void ReconstructMassFluxOrder2XBoundary(std::string BoundaryConditions,
					index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_density,
					const RealType* RESTRICT cell_density_gradient,
					RealType* RESTRICT mass_flux) {
    
    if (BoundaryConditions == "periodic") {
      ReconstructMassFluxOrder2XPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
    }
    else if (BoundaryConditions == "RayleighTaylor") {
      ReconstructMassFluxOrder2XBoundaryRt(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
    }
    else{
      ReconstructMassFluxOrder2XPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
    }   
  }

void ReconstructIntensiveVariableFluxOrder2XBoundary(std::string BoundaryConditions,
						     index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType* RESTRICT volume_fluxes,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     const RealType* RESTRICT cell_variable_gradient,
						     RealType* RESTRICT variable_flux) {
  
  if (BoundaryConditions == "periodic") {
    ReconstructIntensiveVariableFluxOrder2XPeriodicBoundary(nx, ny, halo_width, dx,dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructIntensiveVariableFluxOrder2XBoundaryRt(nx, ny, halo_width, dx,dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }
  else{
    ReconstructIntensiveVariableFluxOrder2XPeriodicBoundary(nx, ny, halo_width, dx,dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }     
}


void ReconstructGradientNodalXBoundary(std::string BoundaryConditions,
				       index_t nx,
				       index_t ny,
				       RealType dx,
				       RealType dt,
				       const RealType* RESTRICT predicted_velocity,
				       const RealType* RESTRICT lag_nodal_variable,
					       RealType* RESTRICT gradient_variable){
  
  if (BoundaryConditions == "periodic") {
    ReconstructGradientNodalXPeriodicBoundary(nx, ny, dx, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructGradientNodalXBoundaryRt(nx, ny, dx, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  else{
    ReconstructGradientNodalXPeriodicBoundary(nx, ny, dx, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  
}

void ProjectNodalIntensiveVariableOrder2XBoundary(std::string BoundaryConditions,
						  index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dx,
						  const RealType dt,
						  const RealType* RESTRICT lag_cell_mass,
						  const RealType* RESTRICT out_cell_mass,
						  const RealType* RESTRICT u_velocity_pred,
						  const RealType* RESTRICT in_variable,
						  const RealType* RESTRICT gradient_variable,
						  const RealType* RESTRICT mass_flux,
						  RealType* RESTRICT out_variable) {
  
  if (BoundaryConditions == "periodic") {
    ProjectNodalIntensiveVariableOrder2XPeriodicBoundary(nx, ny, halo_width, dx, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ProjectNodalIntensiveVariableOrder2XBoundaryRt(nx, ny, halo_width, dx, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  else{
    ProjectNodalIntensiveVariableOrder2XPeriodicBoundary(nx, ny, halo_width, dx, dt, lag_cell_mass, out_cell_mass, u_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  
}

void ReconstructGradientYBoundary(std::string BoundaryConditions,
				  index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT volume_fluxes,
				  const RealType* RESTRICT lag_variable,
				  RealType* RESTRICT gradient_variable) {
  
  if (BoundaryConditions == "periodic") {
    ReconstructGradientYPeriodicBoundary(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructGradientYBoundaryRt(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable); 
  }
  else{
    ReconstructGradientYPeriodicBoundary(nx, ny, dx, dy, volume_fluxes, lag_variable, gradient_variable);
  }
  
}

void ReconstructMassFluxOrder2YBoundary(std::string BoundaryConditions,
					index_t nx, 
					index_t ny, 
					index_t halo_width,
					const RealType dx,
					const RealType dy,
					const RealType* RESTRICT volume_fluxes,
					const RealType* RESTRICT cell_density,
					const RealType* RESTRICT cell_density_gradient,
					RealType* RESTRICT mass_flux) {
  
  if (BoundaryConditions == "periodic") {
    ReconstructMassFluxOrder2YPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructMassFluxOrder2YBoundaryRt(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
  }
  else{
    ReconstructMassFluxOrder2YPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, cell_density, cell_density_gradient, mass_flux);
  }
  
}

void ReconstructIntensiveVariableFluxOrder2YBoundary(std::string BoundaryConditions,
						     index_t nx, 
						     index_t ny, 
						     index_t halo_width,
						     const RealType dx,
						     const RealType dy,
						     const RealType* RESTRICT volume_fluxes,
						     const RealType* RESTRICT mass_flux,
						     const RealType* RESTRICT cell_variable,
						     const RealType* RESTRICT cell_variable_gradient,
						     RealType* RESTRICT variable_flux) {
  
  if (BoundaryConditions == "periodic") {
    ReconstructIntensiveVariableFluxOrder2YPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructIntensiveVariableFluxOrder2YBoundaryRt(nx, ny, halo_width, dx, dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }
  else{
    ReconstructIntensiveVariableFluxOrder2YPeriodicBoundary(nx, ny, halo_width, dx, dy, volume_fluxes, mass_flux, cell_variable, cell_variable_gradient, variable_flux);
  }
  
}

void ReconstructGradientNodalYBoundary(std::string BoundaryConditions,
				       index_t nx,
				       index_t ny,
				       RealType dy,
				       RealType dt,
				       const RealType* RESTRICT predicted_velocity,
				       const RealType* RESTRICT lag_nodal_variable,
				       RealType* RESTRICT gradient_variable) {
  
  if (BoundaryConditions == "periodic") {
    ReconstructGradientNodalYPeriodicBoundary(nx, ny, dy, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ReconstructGradientNodalYBoundaryRt(nx, ny, dy, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  else{
    ReconstructGradientNodalYPeriodicBoundary(nx, ny, dy, dt, predicted_velocity, lag_nodal_variable, gradient_variable);
  }
  
}

void ProjectNodalIntensiveVariableOrder2YBoundary(std::string BoundaryConditions,
						  index_t nx, 
						  index_t ny, 
						  index_t halo_width,
						  const RealType dy,
						  const RealType dt,
						  const RealType* RESTRICT lag_cell_mass,
						  const RealType* RESTRICT out_cell_mass,
						  const RealType* RESTRICT v_velocity_pred,
						  const RealType* RESTRICT in_variable,
						  const RealType* RESTRICT gradient_variable,
						  const RealType* RESTRICT mass_flux,
						  RealType* RESTRICT out_variable) {
  
  if (BoundaryConditions == "periodic") {
    ProjectNodalIntensiveVariableOrder2YPeriodicBoundary(nx, ny, halo_width, dy, dt, lag_cell_mass, out_cell_mass, v_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  else if (BoundaryConditions == "RayleighTaylor") {
    ProjectNodalIntensiveVariableOrder2YBoundaryRt(nx, ny, halo_width, dy, dt, lag_cell_mass, out_cell_mass, v_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  else{
    ProjectNodalIntensiveVariableOrder2YPeriodicBoundary(nx, ny, halo_width, dy, dt, lag_cell_mass, out_cell_mass, v_velocity_pred, in_variable, gradient_variable, mass_flux, out_variable);
  }
  
}

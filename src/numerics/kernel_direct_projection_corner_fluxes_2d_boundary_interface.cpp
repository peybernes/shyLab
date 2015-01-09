// -*- c++ -*- (for emacs users)
#include "kernel_direct_projection_corner_fluxes_2d_boundary_interface.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"



void ReconstructGradientDiagBoundary(std::string BoundaryConditions,
				     index_t nx,
				     index_t ny,
				     RealType dx,
				     RealType dy,
				     RealType dt,
				     const RealType* RESTRICT predicted_u,
				     const RealType* RESTRICT predicted_v,
				     const RealType* RESTRICT lag_variable,
				     RealType* RESTRICT gradient_variable_diag,
				     RealType* RESTRICT gradient_variable_antidiag) {
  
  if (BoundaryConditions == "periodic") {

    ReconstructGradientDiagPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, lag_variable, gradient_variable_diag, gradient_variable_antidiag);

  } else if (BoundaryConditions == "Wall") {

    ReconstructGradientDiagBoundaryWall(nx, ny, dx, dy, dt, predicted_u, predicted_v, lag_variable, gradient_variable_diag, gradient_variable_antidiag);

  } else {

    ReconstructGradientDiagPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, lag_variable, gradient_variable_diag, gradient_variable_antidiag);

  }

} //endReconstructGradientDiagBoundary 


void ReconstructMassFluxOrder2XCornerFluxesBoundary(std::string BoundaryConditions,
						    index_t nx, 
						    index_t ny, 
						    index_t halo_width,
						    const RealType dx,
						    const RealType dy,
						    const RealType dt,
						    const RealType* RESTRICT predicted_u,
						    const RealType* RESTRICT predicted_v,
						    const RealType* RESTRICT volume_fluxes_tot,
						    const RealType* RESTRICT volume_fluxes,
						    const RealType* RESTRICT cell_density,
						    const RealType* RESTRICT cell_density_gradient,
						    RealType* RESTRICT reconstructed_density_faces_x,
						    RealType* RESTRICT mass_flux) {

  if (BoundaryConditions == "periodic") {

    ReconstructMassFluxOrder2XCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_x, mass_flux);

  } else if (BoundaryConditions == "Wall") {

    ReconstructMassFluxOrder2XCornerFluxesBoundaryWall(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_x, mass_flux);

  } else {

    ReconstructMassFluxOrder2XCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_x, mass_flux);

  }

} //end ReconstructMassFluxOrder2XCornerFluxesBoundary



void ReconstructMassFluxOrder2YCornerFluxesBoundary(std::string BoundaryConditions,
						    index_t nx, 
						    index_t ny, 
						    index_t halo_width,
						    const RealType dx,
						    const RealType dy,
						    const RealType dt,
						    const RealType* RESTRICT predicted_u,
						    const RealType* RESTRICT predicted_v,
						    const RealType* RESTRICT volume_fluxes_tot,
						    const RealType* RESTRICT volume_fluxes,
						    const RealType* RESTRICT cell_density,
						    const RealType* RESTRICT cell_density_gradient,
						    RealType* RESTRICT reconstructed_density_faces_y,
						    RealType* RESTRICT mass_flux) {

  if (BoundaryConditions == "periodic") {

    ReconstructMassFluxOrder2YCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_y, mass_flux);

  } else if (BoundaryConditions == "Wall") {

    ReconstructMassFluxOrder2YCornerFluxesBoundaryWall(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_y, mass_flux);

  } else {

    ReconstructMassFluxOrder2YCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v, volume_fluxes_tot, volume_fluxes, cell_density, cell_density_gradient,reconstructed_density_faces_y, mass_flux);

  }

} //end ReconstructMassFluxOrder2YCornerFluxesBoundary




void ReconstructMassFluxOrder2CornerBoundary(std::string BoundaryConditions,
					     index_t nx, 
					     index_t ny,
					     RealType dx,
					     RealType dy,
					     RealType dt,
					     const RealType* RESTRICT directional_lagrangian_density,
					     const RealType* RESTRICT reconstructed_density_faces_x,
					     const RealType* RESTRICT reconstructed_density_faces_y,
					     const RealType* RESTRICT volume_fluxes_corner,
					     const RealType* RESTRICT in_vx, 
					     const RealType* RESTRICT in_vy, 
					     RealType* RESTRICT mass_corner_fluxes) {

  if (BoundaryConditions == "periodic") {
    
    ReconstructMassFluxOrder2CornerPeriodicBoundary(nx, ny, dx, dy, dt, directional_lagrangian_density, reconstructed_density_faces_x, reconstructed_density_faces_y, volume_fluxes_corner, in_vx, in_vy, mass_corner_fluxes);

  } else if (BoundaryConditions == "Wall") {
    
    ReconstructMassFluxOrder2CornerBoundaryWall(nx, ny, dx, dy, dt, directional_lagrangian_density, reconstructed_density_faces_x, reconstructed_density_faces_y, volume_fluxes_corner, in_vx, in_vy, mass_corner_fluxes);

  } else {

    ReconstructMassFluxOrder2CornerPeriodicBoundary(nx, ny, dx, dy, dt, directional_lagrangian_density, reconstructed_density_faces_x, reconstructed_density_faces_y, volume_fluxes_corner, in_vx, in_vy, mass_corner_fluxes);

  }

} //end ReconstructMassFluxorder2CornerBoundary


void ReconstructMassFluxOrder2DiagCornerBoundary(std::string BoundaryConditions,
						 index_t nx, 
						 index_t ny,
						 RealType dx,
						 RealType dy,
						 RealType dt,
						 const RealType* RESTRICT predicted_u,
						 const RealType* RESTRICT predicted_v,
						 const RealType* RESTRICT directional_lagrangian_density,
						 const RealType* RESTRICT volume_fluxes_corner, 
						 const RealType* RESTRICT density_gradient_diag, 
						 const RealType* RESTRICT density_gradient_antidiag, 
						 //out
						 RealType* RESTRICT mass_corner_fluxes) {

  if (BoundaryConditions == "periodic") {

    ReconstructMassFluxOrder2DiagCornerPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, directional_lagrangian_density, volume_fluxes_corner, density_gradient_diag, density_gradient_antidiag, mass_corner_fluxes);

  } else if (BoundaryConditions == "Wall") {

    ReconstructMassFluxOrder2DiagCornerBoundaryWall(nx, ny, dx, dy, dt, predicted_u, predicted_v, directional_lagrangian_density, volume_fluxes_corner, density_gradient_diag, density_gradient_antidiag, mass_corner_fluxes);

  } else {

    ReconstructMassFluxOrder2DiagCornerPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, directional_lagrangian_density, volume_fluxes_corner, density_gradient_diag, density_gradient_antidiag, mass_corner_fluxes);

  }

} //end ReconstructMassFluxorder2DiagCornerBoundary



void ProjectMassCornerFluxesBoundary(std::string BoundaryConditions,
				     index_t nx, 
				     index_t ny,
				     RealType dt,
				     const RealType* RESTRICT mass_corner_fluxes,
				     const RealType* RESTRICT in_vx, 
				     const RealType* RESTRICT in_vy, 
				     RealType* RESTRICT out_cell_mass) {

  if (BoundaryConditions == "periodic") {
    
    ProjectMassCornerFluxesPeriodicBoundary(nx, ny, dt, mass_corner_fluxes, in_vx, in_vy, out_cell_mass);

  } else {

  }

} //end ProjectMassCornerFluxesBoundary



void ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 index_t halo_width,
								 const RealType dx,
								 const RealType dy,
								 const RealType dt,
								 const RealType* RESTRICT predicted_u,
								 const RealType* RESTRICT predicted_v,
								 const RealType* RESTRICT volume_fluxes,
								 const RealType* RESTRICT mass_flux,
								 const RealType* RESTRICT cell_variable,
								 const RealType* RESTRICT cell_variable_gradient,
								 RealType* RESTRICT reconstructed_variable_faces_x,
								 RealType* RESTRICT variable_flux) {

  if (BoundaryConditions == "periodic") {

   ReconstructIntensiveVariableFluxOrder2XCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_x,variable_flux); 

  } else if (BoundaryConditions == "Wall") {

   ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundaryWall(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_x,variable_flux); 

  } else {

   ReconstructIntensiveVariableFluxOrder2XCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_x,variable_flux); 

  } 

}//end ReconstructIntensiveVariableFluxOrder2XCornerFluxesBoundary



void ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary(std::string BoundaryConditions,
								 index_t nx, 
								 index_t ny, 
								 index_t halo_width,
								 const RealType dx,
								 const RealType dy,
								 const RealType dt,
								 const RealType* RESTRICT predicted_u,
								 const RealType* RESTRICT predicted_v,
								 const RealType* RESTRICT volume_fluxes,
								 const RealType* RESTRICT mass_flux,
								 const RealType* RESTRICT cell_variable,
								 const RealType* RESTRICT cell_variable_gradient,
								 RealType* RESTRICT reconstructed_variable_faces_y,
								 RealType* RESTRICT variable_flux) {

  if (BoundaryConditions == "periodic") {

   ReconstructIntensiveVariableFluxOrder2YCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_y,variable_flux); 

  } else if (BoundaryConditions == "Wall") {

   ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundaryWall(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_y,variable_flux); 

  } else {

   ReconstructIntensiveVariableFluxOrder2YCornerFluxesPeriodicBoundary(nx, ny, halo_width, dx, dy, dt, predicted_u, predicted_v,volume_fluxes, mass_flux, cell_variable, cell_variable_gradient,reconstructed_variable_faces_y,variable_flux); 

  }

} //end ReconstructIntensiveVariableFluxOrder2YCornerFluxesBoundary



void ReconstructIntensiveVariableFluxOrder2CornerBoundary(std::string BoundaryConditions,
							  index_t nx, 
							  index_t ny,
							  RealType dx,
							  RealType dy,
							  RealType dt,
							  const RealType* RESTRICT cell_variable,
							  const RealType* RESTRICT reconstructed_variable_faces_x,
							  const RealType* RESTRICT reconstructed_variable_faces_y,
							  const RealType* RESTRICT mass_corner_fluxes,
							  const RealType* RESTRICT in_vx, 
							  const RealType* RESTRICT in_vy, 
							  RealType* RESTRICT variable_flux_corner) {

  if (BoundaryConditions == "periodic") {
    
    ReconstructIntensiveVariableFluxOrder2CornerPeriodicBoundary(nx, ny, dx, dy, dt, cell_variable, reconstructed_variable_faces_x, reconstructed_variable_faces_y, mass_corner_fluxes, in_vx, in_vy, variable_flux_corner);

  } else  if (BoundaryConditions == "Wall") {
    
    ReconstructIntensiveVariableFluxOrder2CornerBoundaryWall(nx, ny, dx, dy, dt, cell_variable, reconstructed_variable_faces_x, reconstructed_variable_faces_y, mass_corner_fluxes, in_vx, in_vy, variable_flux_corner);

  } else  {

    ReconstructIntensiveVariableFluxOrder2CornerPeriodicBoundary(nx, ny, dx, dy, dt, cell_variable, reconstructed_variable_faces_x, reconstructed_variable_faces_y, mass_corner_fluxes, in_vx, in_vy, variable_flux_corner);
  }

} //end ReconstructIntensiveVariableFluxOrder2CornerBoundary


void ReconstructIntensiveVariableFluxOrder2DiagCornerBoundary(std::string BoundaryConditions,
							      index_t nx, 
							      index_t ny,
							      RealType dx,
							      RealType dy,
							      RealType dt,
							      const RealType* RESTRICT predicted_u,
							      const RealType* RESTRICT predicted_v,
							      const RealType* RESTRICT cell_variable,
							      const RealType* RESTRICT mass_corner_fluxes, 
							      const RealType* RESTRICT variable_gradient_diag, 
							      const RealType* RESTRICT variable_gradient_antidiag,
							      RealType* RESTRICT variable_flux_corner) {

  if (BoundaryConditions == "periodic") {
    
    ReconstructIntensiveVariableFluxOrder2DiagCornerPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, cell_variable, mass_corner_fluxes, variable_gradient_diag, variable_gradient_antidiag, variable_flux_corner);

  } else  if (BoundaryConditions == "Wall") {
    
    ReconstructIntensiveVariableFluxOrder2DiagCornerBoundaryWall(nx, ny, dx, dy, dt, predicted_u, predicted_v, cell_variable, mass_corner_fluxes, variable_gradient_diag, variable_gradient_antidiag, variable_flux_corner);

  } else  {

    ReconstructIntensiveVariableFluxOrder2DiagCornerPeriodicBoundary(nx, ny, dx, dy, dt, predicted_u, predicted_v, cell_variable, mass_corner_fluxes, variable_gradient_diag, variable_gradient_antidiag, variable_flux_corner);
  }

} //end ReconstructIntensiveVariableFluxOrder2DiagCornerBoundary



void MassProjectIntensiveVariableCornerFluxesBoundary(std::string BoundaryConditions,
						      index_t nx, 
						      index_t ny,
						      RealType dt,
						      const RealType* RESTRICT in_cell_variable,
						      const RealType* RESTRICT variable_flux_corner, 
						      const RealType* RESTRICT in_vx, 
						      const RealType* RESTRICT in_vy, 
						      const RealType* RESTRICT out_cell_mass,
						      RealType* RESTRICT out_cell_variable) {

  if (BoundaryConditions == "periodic") {
    
    MassProjectIntensiveVariableCornerFluxesPeriodicBoundary(nx, ny, dt, in_cell_variable, variable_flux_corner, in_vx, in_vy, out_cell_mass, out_cell_variable);

  } else {

  }

} //end MassProjectIntensiveVariableCornerFluxesBoundary  




void ProjectNodalIntensiveVariableCornerFluxesBoundary(std::string BoundaryConditions,
						       index_t nx, 
						       index_t ny, 
						       RealType dt,
						       index_t halo_width,
						       const RealType* RESTRICT out_cell_mass,
						       const RealType* RESTRICT in_variable,
						       const RealType* RESTRICT in_u,
						       const RealType* RESTRICT in_v,
						       const RealType* RESTRICT mass_corner_fluxes,
						       const index_t* RESTRICT sign_x_corner_fluxes,
						       const index_t* RESTRICT sign_y_corner_fluxes,
						       RealType* RESTRICT out_variable) {

  if (BoundaryConditions == "periodic") {

    ProjectNodalIntensiveVariableCornerFluxesPeriodicBoundary(nx, ny, dt, halo_width, out_cell_mass, in_variable, in_u, in_v, mass_corner_fluxes, sign_x_corner_fluxes, sign_y_corner_fluxes, out_variable);

  } else {

  }

} //end ProjectNodalIntensiveVariableCornerFluxesBoundary

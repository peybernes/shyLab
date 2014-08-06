// -*- c++ -*- (for emacs users)
#include "kernel_multimat_interface_2d.h"

//#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <math.h>

#include "cartesian_connectivity.h" 
#include "kernel_tools.h"




  // ========================
  //   General functions for Geometry    
  // ========================



int sign(RealType x) {

  return (x > 0.0) ? 1 : ((x < 0.0) ? -1 : 0);

}

void PlaceInterface(// in
		    const RealType xmin,
		    const RealType xmax,
		    const RealType ymin,
		    const RealType ymax,
		    const RealType nX,
		    const RealType nY,
		    const RealType vol_fraction,
		    // out
		    RealType* pointA,
		    RealType* pointB) {

  RealType xA,xB,yA,yB,Nx,Ny,Vol_fraction,tmp,Vol_inf,Vol_sup,a,b;
  int NX_NY,NXNY;

  const RealType nnX = nX;
  const RealType nnY = nY;

  RealType Test_1 = nnX * (xmax - xmin);
  RealType test_1 = fabs(Test_1);
  RealType test_2 = fabs(nnY * (ymax - ymin));

  // If vertical or horizontal normal
  if (test_1 < 1.0e-8) {
 
    xA = 0.0;
    yA = 0.5*( (1+sign(nnY))*vol_fraction + (1-sign(nnY))*(1-vol_fraction) );
    xB = 1.0;
    yB = yA ;

  } else if (test_2 < 1.0e-8) {

    yA = 0.0;
    xA = 0.5*( (1+sign(nnX))*vol_fraction + (1-sign(nnX))*(1-vol_fraction) );
    yB = 1.0;
    xB = xA;        

  } else {
        
    // Normalization, we work on square [0,1]x[0,1]
        Nx = fabs(nnX*(xmax-xmin)) ;
        Ny = fabs(nnY*(ymax-ymin)) ;
        NX_NY = sign(1-Nx/Ny) ;
        NXNY = sign(nnX*nnY) ;
        Vol_fraction = vol_fraction;

        if (nnX < 0.0) {
            Vol_fraction = 1-vol_fraction ;
        }
        
        if (Nx < Ny) {
            tmp = Nx ;
            Nx = Ny ;
            Ny = tmp ;
        }
        
        Vol_inf = 0.5*Ny/Nx ; 
        Vol_sup = 1.0-Vol_inf ;
        a = sqrt(2*Vol_fraction*Ny/Nx) ;
        b = sqrt(2*Vol_fraction*Nx/Ny) ;

	// If smaller than Vol_inf or greater than Vol_sup, we have a triangle        
        if ( (Vol_fraction < Vol_inf) || (Vol_fraction == Vol_inf) ) {

            xA = 0.0 ;
            yA = 0.5*( (1+NX_NY)*0.5*( (1-NXNY)*(1-a) + (1+NXNY)*a ) + (1-NX_NY)*0.5*( (1-NXNY)*(1-b) + (1+NXNY)*b ) ) ;
            xB = 0.5*( (1+NX_NY)*b + (1-NX_NY)*a ) ;
            yB = 0.5*(1-NXNY) ;
            
	} else if (Vol_fraction > Vol_sup) {

            a = sqrt(2*(1-Vol_fraction)*Ny/Nx) ;
            b = sqrt(2*(1-Vol_fraction)*Nx/Ny) ;
            xA = 0.5*( (1+NX_NY)*(1-b) + (1-NX_NY)*(1-a) ) ;
            yA = 0.5*(1+NXNY) ; 
            xB = 1.0 ;
            yB = 0.5*( (1+NX_NY)*0.5*( (1-NXNY)*a + (1+NXNY)*(1-a) ) + (1-NX_NY)*0.5*( (1-NXNY)*b + (1+NXNY)*(1-b) ) );   
            
        } else {

            xA = 0.5*( (1-NX_NY)*(Vol_fraction-Vol_inf) ) ;
            yA = 0.5*( (1-NX_NY)*0.5*(1+NXNY) + (1+NX_NY)*0.5*( (1-NXNY)*(1.0-Vol_fraction-Vol_inf) + (1+NXNY)*(Vol_fraction+Vol_inf) ) ) ;
            xB = 0.5*( (1-NX_NY)*(Vol_fraction+Vol_inf) + (1+NX_NY) ) ; 
            yB = 0.5*( (1-NX_NY)*0.5*(1-NXNY) + (1+NX_NY)*0.5*( (1-NXNY)*(1.0-Vol_fraction+Vol_inf) + (1+NXNY)*(Vol_fraction-Vol_inf) ) ) ;
           
	}            
        
  }
  
  // Inverse transformation of the normalization  
  pointA[0] = xmin + xA*(xmax-xmin) ;
  pointA[1] = ymin + yA*(ymax-ymin) ;
  pointB[0] = xmin + xB*(xmax-xmin) ;
  pointB[1] = ymin + yB*(ymax-ymin) ;
          
} // end PlaceInterface



void IntersectLineRectangle(// in
			    const RealType* point1,
			    const RealType* point2,
			    const RealType xmin,
			    const RealType xmax,
			    const RealType ymin,
			    const RealType ymax,
			    // out
			    RealType* pointA,
			    RealType* pointB) {


  RealType a,b,c,x12,y23,x34,y41,x1,x2,y1,y2;

  x1 = point1[0];
  x2 = point2[0];
  y1 = point1[1];
  y2 = point2[1];

  // In the whole function, if no intersection, pointA[0] = 10^20
  // Default values for x12, x34, y23, y41
  x12 = 1.0e20;
  x34 = 1.0e20;
  y23 = 1.0e20;
  y41 = 1.0e20;

  // Compute the coefficients of the line
  a = y1-y2;
  b = -(x1-x2);
  c = -a*x1 - b*y1;

  // Intersection between the line and the sides of the rectangle
  if (a != 0.0) {
    x12 = -(b*ymin + c)/a ;
    x34 = -(b*ymax + c)/a ;
  }

  if (b != 0.0) {
    y23 = -(a*xmax + c)/b ;
    y41 = -(a*xmin + c)/b ;
  }
 
  // Disjunction of the cases
  if ( ((xmin < x12)||(xmin == x12)) && ((x12 < xmax)||(x12 == xmax)) ) {
    pointA[0] = x12 ;
    pointA[1] = ymin ;
    
    if ( ((ymin < y23)||(ymin == y23)) && ((y23 < ymax)||(y23 == ymax)) ) {
      pointB[0] = xmax ;
      pointB[1] = y23 ;
        
    } else if ( ((xmin < x34)||(xmin == x34)) && ((x34 < xmax)||(x34 == xmax)) ) {
      pointB[0] = x34 ;
      pointB[1] = ymax ;
        
    } else if ( ((ymin < y41)||(ymin == y41)) && ((y41 < ymax)||(y41 == ymax)) ) {
      pointB[0] = xmin ;
      pointB[1] = y41 ;
    
    } else {
      pointA[0] = 1.0e20;
    }
                
  } else if ( ((ymin < y23)||(ymin == y23)) && ((y23 < ymax)||(y23 == ymax)) ) {
    pointA[0] = xmax ;
    pointA[1] = y23 ;
        
    if ( ((xmin < x34)||(xmin == x34)) && ((x34 < xmax)||(x34 == xmax)) ) {
      pointB[0] = x34 ;
      pointB[1] = ymax ;
            
    } else if ( ((ymin < y41)||(ymin == y41)) && ((y41 < ymax)||(y41 == ymax)) ) {
      pointB[0] = xmin ;
      pointB[1] = y41 ;
            
    } else {
      pointA[0] = 1.0e20;
    }
        
  } else if ( ((xmin < x34)||(xmin == x34)) && ((x34 < xmax)||(x34 == xmax)) ) {

    if ( ((ymin < y41)||(ymin == y41)) && ((y41 < ymax)||(y41 == ymax)) ) {
      pointA[0] = x34 ;
      pointA[1] = ymax ;
      pointB[0] = xmin ;
      pointB[1] = y41 ;    
    }
      
  } else {
    pointA[0] = 1.0e20;
  }   

} // end IntersectLineRectangle
    
    

RealType ComputePartialVolume(// in
			      const RealType* pointA,
			      const RealType* pointB,
			      const RealType xmin,
			      const RealType xmax,
			      const RealType ymin,
			      const RealType ymax,			  
			      const RealType nX,			  
			      const RealType nY) {


  RealType XA,XB,YA,YB,Nx,Ny,XI,YI,XC,YC,Vol,ps,volfraction;

  // Normalization
 XA = (pointA[0]-xmin)/(xmax-xmin);   YA = (pointA[1]-ymin)/(ymax-ymin);
 XB = (pointB[0]-xmin)/(xmax-xmin);   YB = (pointB[1]-ymin)/(ymax-ymin);
 Nx = nX*(xmax-xmin);                 Ny = nY*(ymax-ymin);

  // We define I the middle of [AB]
 XI = 0.5*(XA+XB);             YI = 0.5*(YA+YB);

  // If the interface crosses two opposite sides
 if ( (XA+XB == 1.0) && (fabs(XA-XB) == 1.0) ) {

    volfraction = 0.5*( YI*(1+sign(nY)) + (1-YI)*(1-sign(nY)) );

 } else if (YA+YB == 1.0) {

    volfraction = 0.5*( XI*(1+sign(nX)) + (1-XI)*(1-sign(nX)) );
    
    // Else, we have a triangle, and we can find its right angle (vertex C)
 } else {
        
        XC = floor(XA) + floor(XB);
        YC = floor(YA) + floor(YB);
        Vol = 0.5*( fabs(YB-YC)*fabs(XA-XC) + fabs(YA-YC)*fabs(XB-XC) ); 
        
	// The sign of the scalar product helps us choose between the volume of the triangle and its complementary
        ps = (XI-XC)*Nx + (YI-YC)*Ny;
        
        volfraction = 0.5*( (1-sign(ps))*(1-Vol) + (1+sign(ps))*Vol );
        
 }

 return volfraction;

} // end ComputePartialVolume





  // ========================
  //    Functions for Projection phase   
  // ========================


void ReconstructGradientMultimatX(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT in_cell_volumic_fraction,
				  const RealType* RESTRICT volume_fluxes_1,
				  const RealType* RESTRICT volume_fluxes_2,
				  const RealType* RESTRICT lag_variable_1,
				  const RealType* RESTRICT lag_variable_2,
				  RealType* RESTRICT gradient_variable_1,
				  RealType* RESTRICT gradient_variable_2) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
    
      const int cell_ooo = iy * nx + ix;
      const int cell_m1o = CellCellM1O(cell_ooo, nx);
      const int cell_p1o = CellCellP1O(cell_ooo, nx);
      
      const int face_m2o = CellFaceM1O(cell_m1o, iy, nx);
      const int face_m1o = CellFaceM1O(cell_ooo, iy, nx);
      const int face_p1o = CellFaceP1O(cell_ooo, iy, nx);
      const int face_p2o = CellFaceP1O(cell_p1o, iy, nx);
      
#include "reconstruct_variable_gradient_multimat_o2_X.h"

      gradient_variable_1[cell_ooo] = limited_grad_variable_1;
      gradient_variable_2[cell_ooo] = limited_grad_variable_2;
      
    }
  } 

} // end ReconstructGradientMultimatX


void ReconstructGradientMultimatY(index_t nx,
				  index_t ny,
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT in_cell_volumic_fraction,
				  const RealType* RESTRICT volume_fluxes_1,
				  const RealType* RESTRICT volume_fluxes_2,
				  const RealType* RESTRICT lag_variable_1,
				  const RealType* RESTRICT lag_variable_2,
				  RealType* RESTRICT gradient_variable_1,
				  RealType* RESTRICT gradient_variable_2) {
  
#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const int cell_ooo = iy * nx + ix;
      const int cell_om1 = CellCellOM1(cell_ooo, nx);
      const int cell_op1 = CellCellOP1(cell_ooo, nx);
      
      const int face_om2 = CellFaceOM1(cell_om1, iy - 1, nx);
      const int face_om1 = CellFaceOM1(cell_ooo, iy, nx);
      const int face_op1 = CellFaceOP1(cell_ooo, iy, nx);
      const int face_op2 = CellFaceOP1(cell_op1, iy + 1, nx);

#include "reconstruct_variable_gradient_multimat_o2_Y.h"
   
      gradient_variable_1[cell_ooo] = limited_grad_variable_1;
      gradient_variable_2[cell_ooo] = limited_grad_variable_2;

    }
  }

} // end ReconstructGradientMultimatY



void ReconstructNormalInterface(// in
				index_t nx,
				index_t ny,
				RealType dx,
				RealType dy,
				const RealType* RESTRICT in_variable,
				// out
				RealType* RESTRICT interface_normal_x,
				RealType* RESTRICT interface_normal_y) {


#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {
    
      const int cell_oooo = iy * nx + ix;

      if ( (in_variable[cell_oooo] == 0.0) || (in_variable[cell_oooo] == 1.0) ) {

	interface_normal_x[cell_oooo] = 0.0;
	interface_normal_y[cell_oooo] = 0.0;

      } else {
	
	const int cell_m1oo = CellCellM1O(cell_oooo, nx);
	const int cell_p1oo = CellCellP1O(cell_oooo, nx);
	const int cell_oom1 = CellCellOM1(cell_oooo, nx);
	const int cell_oop1 = CellCellOP1(cell_oooo, nx);
	const int cell_m1m1 = CellCellM1M1(cell_oooo, nx);
	const int cell_p1m1 = CellCellP1M1(cell_oooo, nx);
	const int cell_p1p1 = CellCellP1P1(cell_oooo, nx);
	const int cell_m1p1 = CellCellM1P1(cell_oooo, nx);
	
#include "reconstruct_variable_gradient_Youngs.h"
	
	interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
      
      }
    
    }
  }

} // end ReconstructNormalInterface




void ComputeDirectionalLagrangianQuantitiesMultimatX(index_t nx, 
						     index_t ny, 
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     const RealType* RESTRICT in_vx,
						     const RealType* RESTRICT cell_mass_1,
						     const RealType* RESTRICT cell_mass_2,
						     const RealType* RESTRICT cell_volumic_fraction,
						     RealType* RESTRICT volume_fluxes,
						     RealType* RESTRICT directional_lagrangian_volume,
						     RealType* RESTRICT directional_lagrangian_density_1,
						     RealType* RESTRICT directional_lagrangian_density_2) {
#pragma omp parallel 
  {
#pragma omp  for
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx + 1; ++ix) {
      
	const index_t face_ooo = ((nx + 1) * iy) + ix;
      
	const index_t prev_node = FaceXNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceXNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];
	const RealType vx_next = in_vx[next_node];

	const RealType half = 0.5;

	const RealType face_velocity = half * (vx_prev + vx_next);
      
	RealType volume_flux = dt * face_velocity * dy;

	if (fabs(volume_flux) / (dx * dy) < 1.0e-6) {
	  volume_flux = 0.0;
	}
      
	volume_fluxes[face_ooo] = volume_flux;

      }
    }  
    
#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	const RealType volume_flux_prev = volume_fluxes[prev_face];
	const RealType volume_flux_next = volume_fluxes[next_face];
	const RealType current_mass_1 = cell_mass_1[cell_ooo];
	const RealType current_mass_2 = cell_mass_2[cell_ooo];

	const RealType current_directional_lagrangian_volume =
	  (dx * dy) - volume_flux_prev + volume_flux_next;
	const RealType current_directional_lagrangian_volume_1 = vol_fraction_1 * current_directional_lagrangian_volume;
	const RealType current_directional_lagrangian_volume_2 = vol_fraction_2 * current_directional_lagrangian_volume;

	RealType current_directional_lagrangian_density_1, current_directional_lagrangian_density_2;

	if (vol_fraction_1 == 0.0) {
	  current_directional_lagrangian_density_1 = 0.0;
	  current_directional_lagrangian_density_2 = current_mass_2 / current_directional_lagrangian_volume;

	} else if (vol_fraction_2 == 0.0) {
	  current_directional_lagrangian_density_1 = current_mass_1 / current_directional_lagrangian_volume;
	  current_directional_lagrangian_density_2 = 0.0;

	} else {
	  current_directional_lagrangian_density_1 =
	    current_mass_1 / current_directional_lagrangian_volume_1;
	  current_directional_lagrangian_density_2 =
	    current_mass_2 / current_directional_lagrangian_volume_2;

	}

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume;
      
	directional_lagrangian_density_1[cell_ooo] = current_directional_lagrangian_density_1;
	directional_lagrangian_density_2[cell_ooo] = current_directional_lagrangian_density_2;

	assert(0.0 < current_directional_lagrangian_volume);
	assert( (0.0 < current_directional_lagrangian_density_1) || (0.0 == current_directional_lagrangian_density_1) );
	assert( (0.0 < current_directional_lagrangian_density_2) || (0.0 == current_directional_lagrangian_density_2) );
      
      }
    }  
  }
} // end ComputeDirectionalLagrangianQuantitiesMultimatX




void ComputeDirectionalLagrangianQuantitiesMultimatY(index_t nx, 
						     index_t ny, 
						     RealType dt,
						     RealType dx,
						     RealType dy,
						     const RealType* RESTRICT in_vy,
						     const RealType* RESTRICT cell_mass_1,
						     const RealType* RESTRICT cell_mass_2,
						     const RealType* RESTRICT cell_volumic_fraction,
						     RealType* RESTRICT volume_fluxes,
						     RealType* RESTRICT directional_lagrangian_volume,
						     RealType* RESTRICT directional_lagrangian_density_1,
						     RealType* RESTRICT directional_lagrangian_density_2) {
#pragma omp parallel
  {
#pragma omp for
    for (index_t iy = 0; iy < ny + 1; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = (nx * iy) + ix;
      
	const index_t prev_node = FaceYNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceYNodeP1(face_ooo, iy, nx);

	const RealType vy_prev = in_vy[prev_node];
	const RealType vy_next = in_vy[next_node];

	const RealType half = 0.5;
	
	const RealType face_velocity = half * (vy_prev + vy_next);
      
	RealType volume_flux = dt * face_velocity * dx;
      
	if (fabs(volume_flux) / (dx * dy) < 1.0e-6) {
	  volume_flux = 0.0;
	}

	volume_fluxes[face_ooo] = volume_flux;

      }
    }  

#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	const RealType volume_flux_prev = volume_fluxes[prev_face];
	const RealType volume_flux_next = volume_fluxes[next_face];
	const RealType current_mass_1 = cell_mass_1[cell_ooo];
	const RealType current_mass_2 = cell_mass_2[cell_ooo];

	const RealType current_directional_lagrangian_volume =
	  (dx * dy) - volume_flux_prev + volume_flux_next;
	const RealType current_directional_lagrangian_volume_1 = vol_fraction_1 * current_directional_lagrangian_volume;
	const RealType current_directional_lagrangian_volume_2 = vol_fraction_2 * current_directional_lagrangian_volume;

	RealType current_directional_lagrangian_density_1, current_directional_lagrangian_density_2;

	if (vol_fraction_1 == 0.0) {
	  current_directional_lagrangian_density_1 = 0.0;
	  current_directional_lagrangian_density_2 = current_mass_2 / current_directional_lagrangian_volume;

	} else if (vol_fraction_2 == 0.0) {
	  current_directional_lagrangian_density_1 = current_mass_1 / current_directional_lagrangian_volume;
	  current_directional_lagrangian_density_2 = 0.0;

	} else {
	  current_directional_lagrangian_density_1 =
	    current_mass_1 / current_directional_lagrangian_volume_1;
	  current_directional_lagrangian_density_2 =
	    current_mass_2 / current_directional_lagrangian_volume_2;

	}

	directional_lagrangian_volume[cell_ooo] = current_directional_lagrangian_volume;
      
	directional_lagrangian_density_1[cell_ooo] = current_directional_lagrangian_density_1;
	directional_lagrangian_density_2[cell_ooo] = current_directional_lagrangian_density_2;

	assert(0.0 < current_directional_lagrangian_volume);
	assert( (0.0 < current_directional_lagrangian_density_1) || (0.0 == current_directional_lagrangian_density_1) );
	assert( (0.0 < current_directional_lagrangian_density_2) || (0.0 == current_directional_lagrangian_density_2) );
         
      }
    }  
  }
} // end ComputeDirectionalLagrangianQuantitiesMultimatY




void ComputeDirectionalLagrangianFractionalVolumeFluxesX(// in
							 index_t nx, 
							 index_t ny, 
							 RealType dx,
							 RealType dy,
							 const RealType* RESTRICT volume_fluxes,
							 const RealType* RESTRICT directional_lagrangian_volume,
							 const RealType* RESTRICT cell_volumic_fraction,
							 const RealType* RESTRICT interface_normal_x,
							 const RealType* RESTRICT interface_normal_y,
							 // out
							 RealType* RESTRICT volume_fluxes_1,
							 RealType* RESTRICT volume_fluxes_2) {

  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 1; ix < nx; ++ix) {

	const index_t face_ooo = ((nx + 1) * iy) + ix;

	if ( volume_fluxes[face_ooo] == 0.0 ) {

	  volume_fluxes_1[face_ooo] = 0.0;
	  volume_fluxes_2[face_ooo] = 0.0;

	} else {

	 index_t cell_ooo = FaceXCellM1(face_ooo, iy, nx); 

	  if ( volume_fluxes[face_ooo] < 0.0 ) {
	     cell_ooo = FaceXCellP1(face_ooo, iy, nx);
	  }

	  RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	  RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	  RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	  RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);
	 
	  RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_flux = volume_fluxes[face_ooo];

	  // If only one material in cell_ooo
	  if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	    volume_fluxes_1[face_ooo] = vol_fraction_1 * vol_flux;
	    volume_fluxes_2[face_ooo] = vol_fraction_2 * vol_flux;

	  } else { 

	    index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	    index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
	    
	    // Algebraic values of left and right displacement
	    RealType delta_x_prev = volume_fluxes[prev_face] / dy;
	    RealType delta_x_next = volume_fluxes[next_face] / dy;
	    
	    // Vertices of the rectangular approximation of lagrangian cell_ooo
	    RealType xmin = (cell_ooo - cell_ooo / nx * nx ) * dx + delta_x_prev;
	    RealType xmax = (cell_ooo - cell_ooo / nx * nx + 1) * dx + delta_x_next;
	    RealType ymin = cell_ooo / nx * dy;
	    RealType ymax = ymin + dy;

	    // Init of variables for placement of the interface on rectangular approximation of lagrangian cell_ooo
	    RealType normal_x = interface_normal_x[cell_ooo];
	    RealType normal_y = interface_normal_y[cell_ooo];
	    // A and B are the two points that define the interface
	    RealType pointA[2];
	    RealType pointB[2];

	    PlaceInterface(xmin, xmax, ymin, ymax, normal_x, normal_y, vol_fraction_1, pointA, pointB);
	  
	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = ix * dx;
	    RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * x + (1 - sign(vol_flux)) * (x + delta_x_prev) );
	    RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * (x + delta_x_next) + (1 - sign(vol_flux)) * x );

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin, ymax, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    if (vol_flux > 0.0) {
	      if (normal_x < 0.0) {
		volume_fluxes_1[face_ooo] = vol_flux;
		volume_fluxes_2[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1[face_ooo] = 0.0;
		volume_fluxes_2[face_ooo] = vol_flux;
	      }
	    } else {
	      if (normal_x > 0.0) {
		volume_fluxes_1[face_ooo] = vol_flux;
		volume_fluxes_2[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1[face_ooo] = 0.0;
		volume_fluxes_2[face_ooo] = vol_flux;
	      }
	    }
	  	    
	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

	      // If the whole quantity of one material initially in cell_ooo is in dVol

	      if ( (vol_flux > 0.0) && (pointA_d[0] > x) && (pointB_d[0] > x) ) {
		if (normal_x < 0.0) {
		  volume_fluxes_1[face_ooo] = vol_1_lag;
		  volume_fluxes_2[face_ooo] = vol_flux - vol_1_lag;
		} else { 
		  volume_fluxes_1[face_ooo] = vol_flux - vol_2_lag;
		  volume_fluxes_2[face_ooo] = vol_2_lag;
		}
	   
	      } else if ( (vol_flux < 0.0) && (pointA_d[0] < x) && (pointB_d[0] < x) ) {		  
		if (normal_x > 0.0) {
		  volume_fluxes_1[face_ooo] = - vol_1_lag;
		  volume_fluxes_2[face_ooo] = vol_flux + vol_1_lag;
		} else { 
		  volume_fluxes_1[face_ooo] = vol_flux + vol_2_lag;
		  volume_fluxes_2[face_ooo] = - vol_2_lag;
		}
	      
	      } else {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin, ymax, normal_x, normal_y);

		if (vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));

		}

		volume_fluxes_1[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2[face_ooo] = vol_fraction_2_d * vol_flux;
	    	       
	      }
	      
	    }

	  }
	}  
    }
  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesX



void ComputeDirectionalLagrangianFractionalVolumeFluxesY(// in
							 index_t nx, 
							 index_t ny, 
							 RealType dx,
							 RealType dy,
							 const RealType* RESTRICT volume_fluxes_y,
							 const RealType* RESTRICT directional_lagrangian_volume,
							 const RealType* RESTRICT cell_volumic_fraction,
							 const RealType* RESTRICT interface_normal_x,
							 const RealType* RESTRICT interface_normal_y,
							 // out
							 RealType* RESTRICT volume_fluxes_1_y,
							 RealType* RESTRICT volume_fluxes_2_y) {

#pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = nx * iy + ix;

	if ( volume_fluxes_y[face_ooo] == 0.0 ) {

	  volume_fluxes_1_y[face_ooo] = 0.0;
	  volume_fluxes_2_y[face_ooo] = 0.0;

	} else {

	  index_t cell_ooo = FaceYCellM1(face_ooo, iy, nx); 

	  if ( volume_fluxes_y[face_ooo] < 0.0 ) {
	    cell_ooo = FaceYCellP1(face_ooo, iy, nx);
	  }

	  const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	  const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	  const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	  const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	  RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];
	  RealType vol_flux = volume_fluxes_y[face_ooo];

	  // If only one material in cell_ooo
	  if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	    volume_fluxes_1_y[face_ooo] = vol_fraction_1 * volume_fluxes_y[face_ooo];
	    volume_fluxes_2_y[face_ooo] = vol_fraction_2 * volume_fluxes_y[face_ooo];

	  } else { 

	    index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	    index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	    // Algebraic values of up and down displacement
	    RealType delta_y_prev = volume_fluxes_y[prev_face] / dx;
	    RealType delta_y_next = volume_fluxes_y[next_face] / dx;

	    // Vertices of the rectangular approximation of lagrangian cell_ooo
	    RealType xmin = (cell_ooo - cell_ooo / nx * nx ) * dx;
	    RealType xmax = xmin + dx;
	    RealType ymin = cell_ooo / nx * dy + delta_y_prev;
	    RealType ymax = (cell_ooo / nx + 1) * dy + delta_y_next;

	    // Init of variables for placement of the interface on rectangular approximation of lagrangian cell_ooo
	    RealType normal_x = interface_normal_x[cell_ooo];
	    RealType normal_y = interface_normal_y[cell_ooo];
	    RealType pointA[2];
	    RealType pointB[2];
	    
	    PlaceInterface(xmin, xmax, ymin, ymax, normal_x, normal_y, vol_fraction_1, pointA, pointB);

	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType y = iy * dy;
	    RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * y + (1 - sign(vol_flux)) * (y + delta_y_prev) );
	    RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * (y + delta_y_next) + (1 - sign(vol_flux)) * y );

	    IntersectLineRectangle(pointA, pointB, xmin, xmax, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    if (vol_flux > 0.0) {
	      if (normal_y < 0.0) {
		volume_fluxes_1_y[face_ooo] = vol_flux;
		volume_fluxes_2_y[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_y[face_ooo] = 0.0;
		volume_fluxes_2_y[face_ooo] = vol_flux;
	      }
	    } else {
	      if (normal_y > 0.0) {
		volume_fluxes_1_y[face_ooo] = vol_flux;
		volume_fluxes_2_y[face_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_y[face_ooo] = 0.0;
		volume_fluxes_2_y[face_ooo] = vol_flux;
	      }
	    }
	  	    
	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

	      // If the whole quantity of one material initially in cell_ooo is in dVol

	      if ( (vol_flux > 0.0) && (pointA_d[1] > y) && (pointB_d[1] > y) ) {
		if (normal_y < 0.0) {
		  volume_fluxes_1_y[face_ooo] = vol_1_lag;
		  volume_fluxes_2_y[face_ooo] = vol_flux - vol_1_lag;
		} else { 
		  volume_fluxes_1_y[face_ooo] = vol_flux - vol_2_lag;
		  volume_fluxes_2_y[face_ooo] = vol_2_lag;
		}
	   
	      } else if ( (vol_flux < 0.0) && (pointA_d[1] < y) && (pointB_d[1] < y) ) {		  
		if (normal_y > 0.0) {
		  volume_fluxes_1_y[face_ooo] = - vol_1_lag;
		  volume_fluxes_2_y[face_ooo] = vol_flux + vol_1_lag;
		} else { 
		  volume_fluxes_1_y[face_ooo] = vol_flux + vol_2_lag;
		  volume_fluxes_2_y[face_ooo] = - vol_2_lag;
		}
	      
	      } else {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin, xmax, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < 1.0e-4) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		
		}

		volume_fluxes_1_y[face_ooo] = vol_fraction_1_d * vol_flux;
		volume_fluxes_2_y[face_ooo] = vol_fraction_2_d * vol_flux;
	    	       
	      }
	      
	    }

	  }
	}  
    }
  }

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesY




void CheckDirectionalVolumeFluxesX(//in
				   index_t nx, 
				   index_t ny, 
				   RealType dx,
				   RealType dy,
				   const RealType* RESTRICT directional_lagrangian_volume,
				   const RealType* RESTRICT cell_volumic_fraction,
				   // out
				   RealType* RESTRICT volume_fluxes_1,
				   RealType* RESTRICT volume_fluxes_2) {
				   
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      const RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      const RealType prev_vol_flux_1 = volume_fluxes_1[prev_face];
      const RealType next_vol_flux_1 = volume_fluxes_1[next_face];
      const RealType prev_vol_flux_2 = volume_fluxes_2[prev_face];
      const RealType next_vol_flux_2 = volume_fluxes_2[next_face];  

      const RealType vol_proj_1 = vol_1_lag + prev_vol_flux_1 - next_vol_flux_1;
      const RealType vol_proj_2 = vol_2_lag + prev_vol_flux_2 - next_vol_flux_2;

      if ( (vol_proj_1 != 0) && (vol_proj_1 < 1.0e-12)) {

	RealType total_vol_flux_1 = prev_vol_flux_1 - next_vol_flux_1;
	RealType prev_fraction_1 = prev_vol_flux_1 / total_vol_flux_1;
	RealType next_fraction_1 = next_vol_flux_1 / total_vol_flux_1;
	RealType new_total_vol_flux_1 = - vol_1_lag;

	volume_fluxes_1[prev_face] = prev_fraction_1 * new_total_vol_flux_1;
	volume_fluxes_1[next_face] = next_fraction_1 * new_total_vol_flux_1;
      }

      if ( (vol_proj_2 != 0) && (vol_proj_2 < 1.0e-12) ) {

	RealType total_vol_flux_2 = prev_vol_flux_2 - next_vol_flux_2;
	RealType prev_fraction_2 = prev_vol_flux_2 / total_vol_flux_2;
	RealType next_fraction_2 = next_vol_flux_2 / total_vol_flux_2;
	RealType new_total_vol_flux_2 = - vol_2_lag;

	volume_fluxes_2[prev_face] = prev_fraction_2 * new_total_vol_flux_2;
	volume_fluxes_2[next_face] = next_fraction_2 * new_total_vol_flux_2;
      }      

    }
  }
} // end CheckDirectionalVolumeFluxesX




void CheckDirectionalVolumeFluxesY(//in
				   index_t nx, 
				   index_t ny, 
				   RealType dx,
				   RealType dy,
				   const RealType* RESTRICT directional_lagrangian_volume,
				   const RealType* RESTRICT cell_volumic_fraction,
				   // out
				   RealType* RESTRICT volume_fluxes_1_y,
				   RealType* RESTRICT volume_fluxes_2_y) {
				   
#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      const RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      const RealType prev_vol_flux_1 = volume_fluxes_1_y[prev_face];
      const RealType next_vol_flux_1 = volume_fluxes_1_y[next_face];
      const RealType prev_vol_flux_2 = volume_fluxes_2_y[prev_face];
      const RealType next_vol_flux_2 = volume_fluxes_2_y[next_face];      

      const RealType vol_proj_1 = vol_1_lag + prev_vol_flux_1 - next_vol_flux_1;
      const RealType vol_proj_2 = vol_2_lag + prev_vol_flux_2 - next_vol_flux_2;     

      if ( (vol_proj_1 != 0) && (vol_proj_1 < 1.0e-12)) {

	RealType total_vol_flux_1 = prev_vol_flux_1 - next_vol_flux_1;
	RealType prev_fraction_1 = prev_vol_flux_1 / total_vol_flux_1;
	RealType next_fraction_1 = next_vol_flux_1 / total_vol_flux_1;
	RealType new_total_vol_flux_1 = - vol_1_lag;

	volume_fluxes_1_y[prev_face] = prev_fraction_1 * new_total_vol_flux_1;
	volume_fluxes_1_y[next_face] = next_fraction_1 * new_total_vol_flux_1;
      }

      if ( (vol_proj_2 != 0) && (vol_proj_2 < 1.0e-12) ) {

	RealType total_vol_flux_2 = prev_vol_flux_2 - next_vol_flux_2;
	RealType prev_fraction_2 = prev_vol_flux_2 / total_vol_flux_2;
	RealType next_fraction_2 = next_vol_flux_2 / total_vol_flux_2;
	RealType new_total_vol_flux_2 = - vol_2_lag;

	volume_fluxes_2_y[prev_face] = prev_fraction_2 * new_total_vol_flux_2;
	volume_fluxes_2_y[next_face] = next_fraction_2 * new_total_vol_flux_2;
      }      

    }
  }
} // end CheckDirectionalVolumeFluxesY



void CheckDirectionalMassFluxX(//in
			       index_t nx, 
			       index_t ny, 
			       RealType dx,
			       RealType dy,
			       const RealType* RESTRICT cell_mass_1,
			       const RealType* RESTRICT cell_mass_2,
			       const RealType* RESTRICT directional_lagrangian_volume,
			       const RealType* RESTRICT cell_volumic_fraction,
			       const RealType* RESTRICT volume_fluxes_1,
			       const RealType* RESTRICT volume_fluxes_2,
			       // out
			       RealType* RESTRICT mass_flux_1,
			       RealType* RESTRICT mass_flux_2) {
				   
  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      const RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      const RealType prev_vol_flux_1 = volume_fluxes_1[prev_face];
      const RealType next_vol_flux_1 = volume_fluxes_1[next_face];
      const RealType prev_vol_flux_2 = volume_fluxes_2[prev_face];
      const RealType next_vol_flux_2 = volume_fluxes_2[next_face];  

      if (prev_vol_flux_1 == - vol_1_lag) {
	mass_flux_1[prev_face] = - cell_mass_1[cell_ooo];
      }

      if (next_vol_flux_1 == vol_1_lag) {
	mass_flux_1[next_face] = cell_mass_1[cell_ooo];
      }
     
      if (prev_vol_flux_2 == - vol_2_lag) {
	mass_flux_2[prev_face] = - cell_mass_2[cell_ooo];
      }

      if (next_vol_flux_2 == vol_2_lag) {
	mass_flux_2[next_face] = cell_mass_2[cell_ooo];
      }

    }
  }
} // end CheckDirectionalMassFluxX


void CheckDirectionalMassFluxY(//in
			       index_t nx, 
			       index_t ny, 
			       RealType dx,
			       RealType dy,
			       const RealType* RESTRICT cell_mass_1,
			       const RealType* RESTRICT cell_mass_2,
			       const RealType* RESTRICT directional_lagrangian_volume,
			       const RealType* RESTRICT cell_volumic_fraction,
			       const RealType* RESTRICT volume_fluxes_1_y,
			       const RealType* RESTRICT volume_fluxes_2_y,
			       // out
			       RealType* RESTRICT mass_flux_1_y,
			       RealType* RESTRICT mass_flux_2_y) {
				   
  //#pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      const RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      const RealType prev_vol_flux_1 = volume_fluxes_1_y[prev_face];
      const RealType next_vol_flux_1 = volume_fluxes_1_y[next_face];
      const RealType prev_vol_flux_2 = volume_fluxes_2_y[prev_face];
      const RealType next_vol_flux_2 = volume_fluxes_2_y[next_face];  

      if (prev_vol_flux_1 == - vol_1_lag) {
	mass_flux_1_y[prev_face] = - cell_mass_1[cell_ooo];
      }

      if (next_vol_flux_1 == vol_1_lag) {
	mass_flux_1_y[next_face] = cell_mass_1[cell_ooo];
      }
     
      if (prev_vol_flux_2 == - vol_2_lag) {
	mass_flux_2_y[prev_face] = - cell_mass_2[cell_ooo];
      }

      if (next_vol_flux_2 == vol_2_lag) {
	mass_flux_2_y[next_face] = cell_mass_2[cell_ooo];
      }

    }
  }
} // end CheckDirectionalMassFluxY



void ProjectDensityMultimatX(// in
			     index_t nx, 
			     index_t ny, 
			     RealType dx,
			     RealType dy,
			     const RealType* RESTRICT directional_lagrangian_volume,
			     const RealType* RESTRICT in_cell_volumic_fraction,
			     const RealType* RESTRICT volume_fluxes_1_x,			     
			     const RealType* RESTRICT volume_fluxes_2_x,
			     const RealType* RESTRICT in_rho_1,
			     const RealType* RESTRICT in_rho_2,
			     const RealType* RESTRICT in_cell_mass,
			     const RealType* RESTRICT mass_flux_x,
			     // out
			     RealType* RESTRICT out_cell_mass_1,
			     RealType* RESTRICT out_cell_mass_2,
			     RealType* RESTRICT out_cell_volumic_fraction,
			     RealType* RESTRICT out_cell_mass,
			     RealType* RESTRICT out_rho_1,
			     RealType* RESTRICT out_rho_2) {

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;

      const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);

      const RealType out_mass_tot = in_cell_mass[cell_ooo] + mass_flux_x[prev_face] - mass_flux_x[next_face];
      const RealType out_mass_total = out_cell_mass_1[cell_ooo] + out_cell_mass_2[cell_ooo];
      const RealType out_c11 = out_cell_mass_1[cell_ooo] / out_mass_total;
      const RealType out_c22 = out_cell_mass_2[cell_ooo] / out_mass_total;
      const RealType out_c1 = out_c11 / (out_c11 + out_c22);
      const RealType out_c2 = out_c22 / (out_c11 + out_c22);
      const RealType out_mass_1 = out_c1 * out_mass_tot;
      const RealType out_mass_2 = out_c2 * out_mass_tot;

      out_cell_mass_1[cell_ooo] = out_mass_1;
      out_cell_mass_2[cell_ooo] = out_mass_2;
      out_cell_mass[cell_ooo] = out_mass_tot;

      const RealType vol_fraction_11 = in_cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);      

      const RealType lagrangian_volume = directional_lagrangian_volume[cell_ooo];
      const RealType volume_flux_1_prev = volume_fluxes_1_x[prev_face];
      const RealType volume_flux_1_next = volume_fluxes_1_x[next_face];
      const RealType volume_flux_2_prev = volume_fluxes_2_x[prev_face];
      const RealType volume_flux_2_next = volume_fluxes_2_x[next_face];

      RealType out_volume_11 = vol_fraction_1 * lagrangian_volume + volume_flux_1_prev - volume_flux_1_next;
      RealType out_volume_22 = vol_fraction_2 * lagrangian_volume + volume_flux_2_prev - volume_flux_2_next;
      RealType out_vol_fraction_1 =  out_volume_11 / (out_volume_11 + out_volume_22);
      RealType out_vol_fraction_2 =  out_volume_22 / (out_volume_11 + out_volume_22);      

      RealType out_volume_1 = out_vol_fraction_1 * dx * dy;
      RealType out_volume_2 = out_vol_fraction_2 * dx * dy;

      RealType out_density_1, out_density_2;
 
      if (out_vol_fraction_1 < 1.0e-15) {
	
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_cell_mass[cell_ooo];
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if (out_vol_fraction_2 < 1.0e-15) {
	
	out_cell_mass_1[cell_ooo] = out_cell_mass[cell_ooo];
	out_cell_mass_2[cell_ooo] = 0.0;
	out_density_1 =  out_cell_mass_1[cell_ooo] / (dx * dy);
	out_density_2 = 0.0;
	
      } else {

	out_density_1 = out_cell_mass_1[cell_ooo] / out_volume_1;
	out_density_2 = out_cell_mass_2[cell_ooo] / out_volume_2;
	
      }

      out_cell_volumic_fraction[cell_ooo] = out_vol_fraction_1;
      out_rho_1[cell_ooo] = out_density_1;
      out_rho_2[cell_ooo] = out_density_2;

    }
  }

} // end ProjectDensityMultimatX



void ProjectDensityMultimatY(// in
			     index_t nx, 
			     index_t ny, 
			     RealType dx,
			     RealType dy,
			     const RealType* RESTRICT directional_lagrangian_volume,
			     const RealType* RESTRICT in_cell_volumic_fraction,
			     const RealType* RESTRICT volume_fluxes_1_y,			     
			     const RealType* RESTRICT volume_fluxes_2_y,
			     const RealType* RESTRICT in_rho_1,
			     const RealType* RESTRICT in_rho_2,
			     const RealType* RESTRICT in_cell_mass,
			     const RealType* RESTRICT mass_flux_y,
			     // out
			     RealType* RESTRICT out_cell_mass_1,
			     RealType* RESTRICT out_cell_mass_2,
			     RealType* RESTRICT out_cell_volumic_fraction,
			     RealType* RESTRICT out_cell_mass,
			     RealType* RESTRICT out_rho_1,
			     RealType* RESTRICT out_rho_2) {

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = (nx * iy) + ix;

      const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

      const RealType out_mass_tot = in_cell_mass[cell_ooo] + mass_flux_y[prev_face] - mass_flux_y[next_face];
      const RealType out_mass_total = out_cell_mass_1[cell_ooo] + out_cell_mass_2[cell_ooo];
      const RealType out_c11 = out_cell_mass_1[cell_ooo] / out_mass_total;
      const RealType out_c22 = out_cell_mass_2[cell_ooo] / out_mass_total;
      const RealType out_c1 = out_c11 / (out_c11 + out_c22);
      const RealType out_c2 = out_c22 / (out_c11 + out_c22);
      const RealType out_mass_1 = out_c1 * out_mass_tot;
      const RealType out_mass_2 = out_c2 * out_mass_tot;

      out_cell_mass_1[cell_ooo] = out_mass_1;
      out_cell_mass_2[cell_ooo] = out_mass_2;
      out_cell_mass[cell_ooo] = out_mass_tot;

      const RealType vol_fraction_11 = in_cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);      

      const RealType lagrangian_volume = directional_lagrangian_volume[cell_ooo];
      const RealType volume_flux_1_prev = volume_fluxes_1_y[prev_face];
      const RealType volume_flux_1_next = volume_fluxes_1_y[next_face];
      const RealType volume_flux_2_prev = volume_fluxes_2_y[prev_face];
      const RealType volume_flux_2_next = volume_fluxes_2_y[next_face];

      RealType out_vol_fraction_11 = (vol_fraction_1 * lagrangian_volume + volume_flux_1_prev - volume_flux_1_next) / (dx * dy);
      RealType out_vol_fraction_22 = (vol_fraction_2 * lagrangian_volume + volume_flux_2_prev - volume_flux_2_next) / (dx * dy);
      RealType out_vol_fraction_1 = out_vol_fraction_11 / (out_vol_fraction_11 + out_vol_fraction_22);
      RealType out_vol_fraction_2 = out_vol_fraction_22 / (out_vol_fraction_11 + out_vol_fraction_22);      

      RealType out_volume_11 = out_vol_fraction_1 * dx * dy;
      RealType out_volume_22 = out_vol_fraction_2 * dx * dy;
      RealType out_volume_1 = dx * dy * out_volume_11 / (out_volume_11 + out_volume_22);
      RealType out_volume_2 = dx * dy * out_volume_22 / (out_volume_11 + out_volume_22);

      RealType out_density_1, out_density_2;
 
      if (out_vol_fraction_1 < 1.0e-15) {
	
	out_vol_fraction_1 = 0.0;
	out_vol_fraction_2 = 1.0;
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_cell_mass[cell_ooo];
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if (out_vol_fraction_2 < 1.0e-15) {
	
	out_vol_fraction_1 = 1.0;
	out_vol_fraction_2 = 0.0;
	out_cell_mass_1[cell_ooo] = out_cell_mass[cell_ooo];
	out_cell_mass_2[cell_ooo] = 0.0;
	out_density_1 =  out_cell_mass_1[cell_ooo] / (dx * dy);
	out_density_2 = 0.0;
	
      } else {

	out_density_1 = out_cell_mass_1[cell_ooo] / out_volume_1;
	out_density_2 = out_cell_mass_2[cell_ooo] / out_volume_2;
	
      }
      
      out_cell_volumic_fraction[cell_ooo] = out_vol_fraction_1;
      out_rho_1[cell_ooo] = out_density_1;
      out_rho_2[cell_ooo] = out_density_2;

    }
  }

} // end ProjectDensityMultimatY



  // ========================
  //    Functions for Lagrangian phase  
  // ========================


void LagrangePressurePredictedOptimisedMultimat(// in
						int nx,
						int ny,
						RealType dt,
						RealType dx,
						RealType dy,
						RealType gamma_1,
						RealType gamma_2,
						RealType pi_1,
						RealType pi_2,
						const RealType* RESTRICT in_mass_1,
						const RealType* RESTRICT in_mass_2,
						const RealType* RESTRICT in_rho_1,
						const RealType* RESTRICT in_rho_2,
						const RealType* RESTRICT in_energy_1,	 
						const RealType* RESTRICT in_energy_2,	 
						const RealType* RESTRICT cell_volumic_fraction,	 
						const RealType* RESTRICT in_velocity_x,
						const RealType* RESTRICT in_velocity_y,
						// out
						RealType* RESTRICT out_pressure,
						RealType* RESTRICT out_pressure_1,
						RealType* RESTRICT out_pressure_2,
						RealType* RESTRICT out_predicted_pressure,
						RealType* RESTRICT out_predicted_pressure_1,
						RealType* RESTRICT out_predicted_pressure_2,
						RealType* RESTRICT out_pseudo_pressure) {
  
  const RealType half = 0.5;
  const RealType one = 1.0;

  const RealType one_over_dx = one / dx;
  const RealType one_over_dy = one / dy;
  
#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {

      SHY_ASM_COMMENT("LagrangePressurePredicted -- INNER LOOP BEGIN");

      // BEGIN DATA LOAD.
      const int cell_ooo = nx * iy + ix; 
   
      const int node_SW = CellNodeM1M1(cell_ooo, iy, nx);
      const int node_SE = CellNodeP1M1(cell_ooo, iy, nx);
      const int node_NW = CellNodeM1P1(cell_ooo, iy, nx);
      const int node_NE = CellNodeP1P1(cell_ooo, iy, nx);

      const RealType ux_sw = in_velocity_x[node_SW];
      const RealType ux_se = in_velocity_x[node_SE];
      const RealType ux_nw = in_velocity_x[node_NW];
      const RealType ux_ne = in_velocity_x[node_NE];

      const RealType uy_sw = in_velocity_y[node_SW];
      const RealType uy_se = in_velocity_y[node_SE];
      const RealType uy_nw = in_velocity_y[node_NW];
      const RealType uy_ne = in_velocity_y[node_NE];

      const RealType mass_1_ooo = in_mass_1[cell_ooo];
      const RealType mass_2_ooo = in_mass_2[cell_ooo];
      const RealType mass_ooo = mass_1_ooo + mass_2_ooo;

      const RealType e_1_ooo = in_energy_1[cell_ooo];
      const RealType e_2_ooo = in_energy_2[cell_ooo];

      const RealType rho_1_ooo = in_rho_1[cell_ooo];     
      const RealType rho_2_ooo = in_rho_2[cell_ooo]; 
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      // END DATA LOAD.

      // BEGIN COMPUTE.
            
      const RealType rho_ooo = one_over_dx * one_over_dy * mass_ooo; // 2 MUL            
      const RealType one_over_rho_1_ooo = one / rho_1_ooo; // 1 DIV
      const RealType one_over_rho_2_ooo = one / rho_2_ooo; // 1 DIV
      const RealType one_over_rho_ooo = one / rho_ooo; // 1 DIV

      const RealType p_1_ooo = EquationOfState(gamma_1, rho_1_ooo, e_1_ooo, pi_1); // 1 MUL, 1 FMA
      const RealType p_2_ooo = EquationOfState(gamma_2, rho_2_ooo, e_2_ooo, pi_2); // 1 MUL, 1 FMA
      const RealType p_ooo = vol_fraction_1 * p_1_ooo + vol_fraction_2 * p_2_ooo;
      
      const RealType delta_ux = half * ((ux_se + ux_ne) - (ux_sw + ux_nw)); // 2 ADD, 1 FMA
      const RealType delta_uy = half * ((uy_nw + uy_ne) - (uy_sw + uy_se)); // 2 ADD, 1 FMA
      
      const RealType delta_volume = half * dt * (delta_ux * dy + delta_uy * dx); // 3 MUL, 1 FMA
      
      const RealType delta_velocity = delta_ux + delta_uy; // 1 ADD

      // Formulas below valid for perfect gas law and stiffened gas.
      const RealType cs_1 = SpeedOfSound(gamma_1, rho_ooo, p_ooo, pi_1);
      const RealType cs_2 = SpeedOfSound(gamma_2, rho_ooo, p_ooo, pi_2);
      const RealType cs  = std::max(cs_1, cs_2); 

      // Pseudo viscosity (same for both materials).
      const RealType linear_pseudo_coeff = one;
      const RealType quadratic_pseudo_coeff = one;

      const RealType delta_velocity_minus = - half * (delta_velocity - fabs(delta_velocity)); // 1 MUL, 1 ADD, 1 ABS
      
      const RealType q_ooo = rho_ooo * delta_velocity_minus * 
	((linear_pseudo_coeff * cs) + (quadratic_pseudo_coeff * delta_velocity_minus)); // 3 MUL, 1 FMA
      
      const RealType div_u_ooo =
	(one_over_dx * delta_ux + one_over_dy * delta_uy); // 1 MUL, 1 FMA

      RealType e_1_lag_ooo, e_2_lag_ooo, predicted_rho_1_ooo, predicted_rho_2_ooo;
      
      if (vol_fraction_1 == 0.0) {

	e_1_lag_ooo = 0.0; 
	e_2_lag_ooo = e_2_ooo 
	  - half * dt *  vol_fraction_2 * (p_2_ooo + q_ooo) * div_u_ooo * one_over_rho_2_ooo; // 1 FMA, 3 MUL, 1 ADD
	predicted_rho_1_ooo = 0.0;
	predicted_rho_2_ooo = mass_2_ooo / (dx * dy + delta_volume) ; // 1 DIV, 1 FMA

      } else if (vol_fraction_2 == 0.0) {

	e_1_lag_ooo = e_1_ooo 
	  - half * dt * vol_fraction_1 * (p_1_ooo + q_ooo) * div_u_ooo * one_over_rho_1_ooo; // 1 FMA, 3 MUL, 1 ADD     
	e_2_lag_ooo = 0.0; 
	predicted_rho_1_ooo = mass_1_ooo / (dx * dy + delta_volume); // 1 DIV, 1 FMA
	predicted_rho_2_ooo = 0.0;

      } else {

	e_1_lag_ooo = e_1_ooo 
	  - half * dt * vol_fraction_1 * (p_1_ooo + q_ooo) * div_u_ooo * one_over_rho_1_ooo; // 1 FMA, 3 MUL, 1 ADD
	e_2_lag_ooo = e_2_ooo 
	  - half * dt *  vol_fraction_2 * (p_2_ooo + q_ooo) * div_u_ooo * one_over_rho_2_ooo; // 1 FMA, 3 MUL, 1 ADD     
	predicted_rho_1_ooo = mass_1_ooo / ( (dx * dy + delta_volume) * vol_fraction_1 ); // 1 DIV, 1 FMA
	predicted_rho_2_ooo = mass_2_ooo / ( (dx * dy + delta_volume) * vol_fraction_2 ); // 1 DIV, 1 FMA

      }

      const RealType out_predicted_p_1_ooo = EquationOfState(gamma_1, predicted_rho_1_ooo, e_1_lag_ooo, pi_1); // 1 MUL, 1 FMA
      const RealType out_predicted_p_2_ooo = EquationOfState(gamma_2, predicted_rho_2_ooo, e_2_lag_ooo, pi_2); // 1 MUL, 1 FMA
      const RealType out_predicted_p_ooo = vol_fraction_1 * out_predicted_p_1_ooo + vol_fraction_2 * out_predicted_p_2_ooo;

      // END COMPUTE.
      // Summary : 7 ADD, 17 MUL, 9 FMA, 2 DIV, 1 ABS, 1 SQRT
      // Flop = 7 + 17 + (9 x 2) + (2 x 7) + 7 = 63 (indicative only...).
      // 3 loads, 3 stores = 24 bytes loads, 24 bytes stores.

      // BEGIN DATA STORE.
      out_pressure[cell_ooo] = p_ooo ;
      out_pressure_1[cell_ooo] = p_1_ooo ;
      out_pressure_2[cell_ooo] = p_2_ooo ;
      out_predicted_pressure[cell_ooo] = out_predicted_p_ooo;
      out_predicted_pressure_1[cell_ooo] = out_predicted_p_1_ooo;
      out_predicted_pressure_2[cell_ooo] = out_predicted_p_2_ooo;
      out_pseudo_pressure[cell_ooo] = q_ooo;
      // END DATA STORE.

      SHY_ASM_COMMENT("LagrangePressurePredicted - INNER LOOP END");

    }
  }
} // end LagrangePressurePredictedOptimisedMultimat



void LagrangeCorrectionOptimisedMultimat(// in
					 int nx,
					 int ny,
					 RealType dt,
					 RealType dx,
					 RealType dy,
					 const RealType* RESTRICT in_mass_1,
					 const RealType* RESTRICT in_mass_2,
					 const RealType* RESTRICT in_energy_1,
					 const RealType* RESTRICT in_energy_2,
					 const RealType* RESTRICT in_pressure_1,
					 const RealType* RESTRICT in_pressure_2,
					 const RealType* RESTRICT in_pseudo_pressure,
					 const RealType* RESTRICT cell_volumic_fraction,
					 const RealType* RESTRICT in_velocity_x,
					 const RealType* RESTRICT in_velocity_y,
					 // out
					 RealType* RESTRICT out_energy,
					 RealType* RESTRICT out_energy_1,
					 RealType* RESTRICT out_energy_2) {

  const RealType one = 1.0;
  const RealType one_over_dx = one / dx;
  const RealType one_over_dy = one / dy;

#pragma omp parallel for
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx ; ++ix) {

      SHY_ASM_COMMENT("LagrangeCorrection -- INNER LOOP BEGIN");

      //DATA LOAD
      const int cell_ooo = nx * iy + ix;

      const int node_SW = CellNodeM1M1(cell_ooo,iy,nx);
      const int node_SE = CellNodeP1M1(cell_ooo,iy,nx);
      const int node_NW = CellNodeM1P1(cell_ooo,iy,nx);
      const int node_NE = CellNodeP1P1(cell_ooo,iy,nx);
      
      const RealType ux_sw = in_velocity_x[node_SW];
      const RealType ux_se = in_velocity_x[node_SE];
      const RealType ux_nw = in_velocity_x[node_NW];
      const RealType ux_ne = in_velocity_x[node_NE];

      const RealType uy_sw = in_velocity_y[node_SW];
      const RealType uy_se = in_velocity_y[node_SE];
      const RealType uy_nw = in_velocity_y[node_NW];
      const RealType uy_ne = in_velocity_y[node_NE];

      const RealType mass_1_ooo = in_mass_1[cell_ooo];
      const RealType mass_2_ooo = in_mass_2[cell_ooo];
      const RealType mass_ooo = mass_1_ooo + mass_2_ooo;

      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType e_1_ooo = in_energy_1[cell_ooo];
      const RealType e_2_ooo = in_energy_2[cell_ooo];

      const RealType p_1_ooo = in_pressure_1[cell_ooo];
      const RealType p_2_ooo = in_pressure_2[cell_ooo];

      const RealType q_ooo = in_pseudo_pressure[cell_ooo];

      const RealType half = 0.5;
      
      const RealType one_over_mass_1_ooo = one / mass_1_ooo; // 1 DIV
      const RealType one_over_mass_2_ooo = one / mass_2_ooo; // 1 DIV
      const RealType one_over_mass_ooo = one / mass_ooo; // 1 DIV

      const RealType div_u_ooo = half * // 3 MUL, 6 ADD
      (one_over_dx * (ux_se + ux_ne - ux_sw - ux_nw) +
       one_over_dy * (uy_nw + uy_ne - uy_sw - uy_se));

      RealType e_1_lag_ooo, e_2_lag_ooo;
      
      if (vol_fraction_1 == 0.0) {

	e_1_lag_ooo = 0.0;
	e_2_lag_ooo = e_2_ooo // 2 MUL, 4 ADD
	  - dt * vol_fraction_2 * (p_2_ooo + q_ooo) * div_u_ooo * dx * dy * one_over_mass_2_ooo; 

      } else if (vol_fraction_2 == 0.0) {

	e_1_lag_ooo = e_1_ooo // 2 MUL, 4 ADD
	  - dt * vol_fraction_1 * (p_1_ooo + q_ooo) * div_u_ooo * dx * dy * one_over_mass_1_ooo; 
	e_2_lag_ooo = 0.0;

      } else {

	e_1_lag_ooo = e_1_ooo // 2 MUL, 4 ADD
	  - dt * vol_fraction_1 * (p_1_ooo + q_ooo) * div_u_ooo * vol_fraction_1 * dx * dy * one_over_mass_1_ooo;      
	e_2_lag_ooo = e_2_ooo // 2 MUL, 4 ADD
	  - dt * vol_fraction_2 * (p_2_ooo + q_ooo) * div_u_ooo * vol_fraction_2 * dx * dy * one_over_mass_2_ooo;

      }

      const RealType e_lag_ooo = (mass_1_ooo * e_1_lag_ooo + mass_2_ooo * e_2_lag_ooo) * one_over_mass_ooo;

      out_energy[cell_ooo] = e_lag_ooo;
      out_energy_1[cell_ooo] = e_1_lag_ooo;
      out_energy_2[cell_ooo] = e_2_lag_ooo;

      // 10 ADD, 5 MUL, 1 DIV : approx 22 FLOPS
      SHY_ASM_COMMENT("LagrangeCorrection -- INNER LOOP END");
    }
  }
} // end LagrangeCorrectionOptimisedMultimat

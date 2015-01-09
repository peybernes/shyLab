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



index_t sign(RealType x) {

  return (x > 0.0) ? 1 : ((x < 0.0) ? -1 : 0);

}

RealType modulo(RealType x, RealType y) {
  
  RealType res = x - y * floor(x/y);
  return res;

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


RealType IntersectLineHorizontalLine(// in
				 const RealType* point1,
				 const RealType* point2,
				 const RealType y) {


  RealType a,b,c,x1,x2,y1,y2,x;

  x1 = point1[0];
  x2 = point2[0];
  y1 = point1[1];
  y2 = point2[1];

  // Compute the coefficients of the line
  a = y1-y2;
  b = -(x1-x2);
  c = -a*x1 - b*y1;

  if (a == 0.0) {
    x = 1.0e20;
  } else {
    x = - (b * y + c) / a;
  }

  return x;
  
} //end IntersectLineHorizontalLine



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



  // AdProjection
  // ------------


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



void ReconstructGradientMultimatDiag(index_t nx,
				     index_t ny,
				     RealType dx,
				     RealType dy,
				     RealType dt,
				     const RealType* RESTRICT predicted_u,
				     const RealType* RESTRICT predicted_v,
				     const RealType* RESTRICT in_cell_volumic_fraction,
				     const RealType* RESTRICT volume_fluxes_1_corner,
				     const RealType* RESTRICT volume_fluxes_2_corner,
				     const RealType* RESTRICT lag_variable_1,
				     const RealType* RESTRICT lag_variable_2,
				     RealType* RESTRICT gradient_variable_1_diag,
				     RealType* RESTRICT gradient_variable_2_diag,
				     RealType* RESTRICT gradient_variable_1_antidiag,
				     RealType* RESTRICT gradient_variable_2_antidiag) {


#pragma omp parallel for
  for (index_t iy = 1; iy < ny - 1; ++iy) {
    for (index_t ix = 1; ix < nx - 1; ++ix) {

      const index_t cell_ooo = iy * nx + ix;
      const index_t cell_m1m1 = CellCellM1M1(cell_ooo, nx);
      const index_t cell_m1p1 = CellCellM1P1(cell_ooo, nx);
      const index_t cell_p1m1 = CellCellP1M1(cell_ooo, nx);
      const index_t cell_p1p1 = CellCellP1P1(cell_ooo, nx);

      const index_t node_mm1mm1 = CellNodeM1M1(cell_m1m1, iy - 1, nx);
      const index_t node_pp1mm1 = CellNodeP1M1(cell_p1m1, iy - 1, nx);
      const index_t node_mm1pp1 = CellNodeM1P1(cell_m1p1, iy + 1, nx);
      const index_t node_pp1pp1 = CellNodeP1P1(cell_p1p1, iy + 1, nx);

#include "reconstruct_variable_gradient_multimat_diag.h"

      gradient_variable_1_diag[cell_ooo] = limited_grad_variable_1_diag; // 1 store
      gradient_variable_2_diag[cell_ooo] = limited_grad_variable_2_diag; // 1 store
      gradient_variable_1_antidiag[cell_ooo] = limited_grad_variable_1_antidiag; // 1 store
      gradient_variable_2_antidiag[cell_ooo] = limited_grad_variable_2_antidiag; // 1 store
     
    }
  }
 
}//end ReconstructGradientMultimatDiag




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

	if ((gradient_var_x == 0.0) && (gradient_var_y == 0.0))	{
	  interface_normal_x[cell_oooo] = 1.0 / sqrt(2);
	  interface_normal_y[cell_oooo] = 1.0 / sqrt(2);

	} else {	  
	  interface_normal_x[cell_oooo] = - gradient_var_x / norm_gradient_var;
	  interface_normal_y[cell_oooo] = - gradient_var_y / norm_gradient_var;
	}
      
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

      if ( (vol_proj_1 != 0) && (vol_proj_1 / (dx * dy) < 1.0e-10)) {

	RealType total_vol_flux_1 = prev_vol_flux_1 - next_vol_flux_1;
	RealType prev_fraction_1 = prev_vol_flux_1 / total_vol_flux_1;
	RealType next_fraction_1 = next_vol_flux_1 / total_vol_flux_1;
	RealType new_total_vol_flux_1 = - vol_1_lag;

	volume_fluxes_1[prev_face] = prev_fraction_1 * new_total_vol_flux_1;
	volume_fluxes_1[next_face] = next_fraction_1 * new_total_vol_flux_1;
     }

      if ( (vol_proj_2 != 0) && (vol_proj_2 / (dx * dy) < 1.0e-10) ) {

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

      if ( (vol_proj_1 != 0) && (vol_proj_1 / (dx * dy) < 1.0e-10)) {

	RealType total_vol_flux_1 = prev_vol_flux_1 - next_vol_flux_1;
	RealType prev_fraction_1 = prev_vol_flux_1 / total_vol_flux_1;
	RealType next_fraction_1 = next_vol_flux_1 / total_vol_flux_1;
	RealType new_total_vol_flux_1 = - vol_1_lag;

	volume_fluxes_1_y[prev_face] = prev_fraction_1 * new_total_vol_flux_1;
	volume_fluxes_1_y[next_face] = next_fraction_1 * new_total_vol_flux_1;
      }

      if ( (vol_proj_2 != 0) && (vol_proj_2 / (dx * dy) < 1.0e-10) ) {

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
 
      if ((out_vol_fraction_1 < 1.0e-12) || (out_mass_1 == 0.0)) {
	
	out_vol_fraction_1 = 0.0;
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_cell_mass[cell_ooo];
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if ((out_vol_fraction_2 < 1.0e-12) || (out_mass_2 == 0.0)) {
	
	out_vol_fraction_1 = 1.0;
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
 
      if ((out_vol_fraction_1 < 1.0e-12) || (out_mass_1 == 0.0)) {
	
	out_vol_fraction_1 = 0.0;
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_cell_mass[cell_ooo];
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if ((out_vol_fraction_2 < 1.0e-12) || (out_mass_2 == 0.0)) {
	
	out_vol_fraction_1 = 1.0;
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





  // DirectProjection
  // ----------------


void CheckDirectionalVolumeFluxesXY(//in
				    index_t nx, 
				    index_t ny, 
				    RealType dx,
				    RealType dy,
				    const RealType* RESTRICT volume_fluxes_x,
				    const RealType* RESTRICT volume_fluxes_y,
				    const RealType* RESTRICT cell_volumic_fraction,
				    // out
				    RealType* RESTRICT bool_check_fluxes_x,
				    RealType* RESTRICT bool_check_fluxes_y,
				    RealType* RESTRICT volume_fluxes_1_x,
				    RealType* RESTRICT volume_fluxes_1_y,
				    RealType* RESTRICT volume_fluxes_2_x,
				    RealType* RESTRICT volume_fluxes_2_y) {
				   
  #pragma omp parallel for
  for (index_t face_x = 0; face_x < (nx + 1) * ny; ++face_x) {
    bool_check_fluxes_x[face_x] = 0.0;
  }

  #pragma omp parallel for
  for (index_t face_y = 0; face_y < (ny + 1) * nx; ++face_y) {
    bool_check_fluxes_y[face_y] = 0.0;
  }

  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_lag = dx * dy + volume_fluxes_x[prev_face_x] - volume_fluxes_x[next_face_x] + volume_fluxes_y[prev_face_y] - volume_fluxes_y[next_face_y];
      const RealType vol_1_lag = vol_fraction_1 * vol_lag;
      const RealType vol_2_lag = vol_fraction_2 * vol_lag;

      const RealType prev_vol_flux_1_x = volume_fluxes_1_x[prev_face_x];
      const RealType next_vol_flux_1_x = volume_fluxes_1_x[next_face_x];
      const RealType prev_vol_flux_2_x = volume_fluxes_2_x[prev_face_x];
      const RealType next_vol_flux_2_x = volume_fluxes_2_x[next_face_x];  
      const RealType prev_vol_flux_1_y = volume_fluxes_1_y[prev_face_y];
      const RealType next_vol_flux_1_y = volume_fluxes_1_y[next_face_y];
      const RealType prev_vol_flux_2_y = volume_fluxes_2_y[prev_face_y];
      const RealType next_vol_flux_2_y = volume_fluxes_2_y[next_face_y];

      const RealType prev_bool_check_fluxes_x = bool_check_fluxes_x[prev_face_x];
      const RealType next_bool_check_fluxes_x = bool_check_fluxes_x[next_face_x];
      const RealType prev_bool_check_fluxes_y = bool_check_fluxes_y[prev_face_y];
      const RealType next_bool_check_fluxes_y = bool_check_fluxes_y[next_face_y];

      const RealType vol_proj_1 = vol_1_lag + prev_vol_flux_1_x - next_vol_flux_1_x + prev_vol_flux_1_y - next_vol_flux_1_y;
      const RealType vol_proj_2 = vol_2_lag + prev_vol_flux_2_x - next_vol_flux_2_x + prev_vol_flux_2_y - next_vol_flux_2_y;
  
      if ( (vol_proj_1 != 0) && (vol_proj_1 / (dx * dy) < 1.0e-10)) {

	RealType total_vol_flux_1 = prev_vol_flux_1_x - next_vol_flux_1_x + prev_vol_flux_1_y - next_vol_flux_1_y;
	RealType prev_fraction_1_x = prev_vol_flux_1_x / total_vol_flux_1;
	RealType next_fraction_1_x = next_vol_flux_1_x / total_vol_flux_1;
	RealType prev_fraction_1_y = prev_vol_flux_1_y / total_vol_flux_1;
	RealType next_fraction_1_y = next_vol_flux_1_y / total_vol_flux_1;
	RealType new_total_vol_flux_1 = - vol_1_lag;
	RealType remaining_vol_flux_1, new_remaining_vol_flux_1, new_prev_fraction_1_x, new_next_fraction_1_x, new_prev_fraction_1_y, new_next_fraction_1_y;

	if (prev_bool_check_fluxes_x == 1.0) {

	  if (prev_bool_check_fluxes_y == 1.0) {

	    if (next_bool_check_fluxes_x == 1.0) {

	      volume_fluxes_1_y[next_face_y] = new_total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_x + prev_vol_flux_1_y);
	      bool_check_fluxes_y[next_face_y] = 1.0;

	    } else if (next_bool_check_fluxes_y == 1.0) {
	      
	      volume_fluxes_1_x[next_face_x] = new_total_vol_flux_1 - (prev_vol_flux_1_x + prev_vol_flux_1_y - next_vol_flux_1_y);
	      bool_check_fluxes_x[next_face_x] = 1.0;

	    } else {

	      remaining_vol_flux_1 = total_vol_flux_1 - (prev_vol_flux_1_x + prev_vol_flux_1_y);
	      new_remaining_vol_flux_1 = new_total_vol_flux_1 - (prev_vol_flux_1_x + prev_vol_flux_1_y);
	      new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1; 
	      new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	      volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	      volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	      bool_check_fluxes_x[next_face_x] = 1.0;
	      bool_check_fluxes_y[next_face_y] = 1.0;	      
	    }

	  } else if (next_bool_check_fluxes_y == 1.0) {

	    if (next_bool_check_fluxes_x == 1.0) {

	      volume_fluxes_1_y[prev_face_y] = new_total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_x - next_vol_flux_1_y);
	      bool_check_fluxes_y[prev_face_y] = 1.0;

	    } else {

	      remaining_vol_flux_1 = total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_y);
	      new_remaining_vol_flux_1 = new_total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_y);
	      new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1;
	      new_prev_fraction_1_y = prev_vol_flux_1_y / remaining_vol_flux_1;
	      volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	      volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	      bool_check_fluxes_x[next_face_x] = 1.0;
	      bool_check_fluxes_y[prev_face_y] = 1.0;
	    }

	  } else if (next_bool_check_fluxes_x == 1.0) {

	    remaining_vol_flux_1 = total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_x);
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 - (prev_vol_flux_1_x - next_vol_flux_1_x);
	    new_prev_fraction_1_y = prev_vol_flux_1_x / remaining_vol_flux_1;
	    new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	    volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;

	  } else {

	    remaining_vol_flux_1 = total_vol_flux_1 - prev_vol_flux_1_x;
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 - prev_vol_flux_1_x;
	    new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1;
	    new_prev_fraction_1_y = prev_vol_flux_1_y / remaining_vol_flux_1;
	    new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	    volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	    bool_check_fluxes_x[next_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }

	} else if (prev_bool_check_fluxes_y == 1.0) {

	  if (next_bool_check_fluxes_x == 1.0) {
	    
	    if (next_bool_check_fluxes_y == 1.0) {
	      
	      volume_fluxes_1_x[prev_face_x] = new_total_vol_flux_1 - (- next_vol_flux_1_x + prev_vol_flux_1_y - next_vol_flux_1_y);
	      bool_check_fluxes_x[prev_face_x] = 1.0;

	    } else {
	      
	      remaining_vol_flux_1 = total_vol_flux_1 - (prev_vol_flux_1_y - next_vol_flux_1_x);
	      new_remaining_vol_flux_1 = new_total_vol_flux_1 - (prev_vol_flux_1_y - next_vol_flux_1_x);
	      new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	      new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	      volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	      volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	      bool_check_fluxes_x[prev_face_x] = 1.0;
	      bool_check_fluxes_y[next_face_y] = 1.0;
	    }

	  } else if (next_bool_check_fluxes_y == 1.0) {
	    
	    remaining_vol_flux_1 = total_vol_flux_1 - (prev_vol_flux_1_y - next_vol_flux_1_y);
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 - (prev_vol_flux_1_y - next_vol_flux_1_y);
	    new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	    new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1;
	    volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_x[next_face_x] = 1.0;

	  } else {
	    
	    remaining_vol_flux_1 = total_vol_flux_1 - prev_vol_flux_1_y;
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 - prev_vol_flux_1_y;
	    new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1;
	    new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	    new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	    volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	    bool_check_fluxes_x[next_face_x] = 1.0;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }
	  
	} else if (next_bool_check_fluxes_x == 1.0) {
	  
	  if (next_bool_check_fluxes_y == 1.0) {
	    
	    remaining_vol_flux_1 = total_vol_flux_1 - (- next_vol_flux_1_x - next_vol_flux_1_y);
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 - (- next_vol_flux_1_x - next_vol_flux_1_y);
	    new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	    new_prev_fraction_1_y = prev_vol_flux_1_y / remaining_vol_flux_1;
	    volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;

	  } else {
	    
	    remaining_vol_flux_1 = total_vol_flux_1 + next_vol_flux_1_x;
	    new_remaining_vol_flux_1 = new_total_vol_flux_1 + next_vol_flux_1_x;
	    new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	    new_prev_fraction_1_y = prev_vol_flux_1_y / remaining_vol_flux_1;
	    new_next_fraction_1_y = next_vol_flux_1_y / remaining_vol_flux_1;
	    volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	    volume_fluxes_1_y[next_face_y] = new_next_fraction_1_y * new_remaining_vol_flux_1;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }
	  
	} else if (next_bool_check_fluxes_y == 1.0) {
	  
	  remaining_vol_flux_1 = total_vol_flux_1 + next_vol_flux_1_y;
	  new_remaining_vol_flux_1 = new_total_vol_flux_1 + next_vol_flux_1_y;
	  new_next_fraction_1_x = next_vol_flux_1_x / remaining_vol_flux_1;
	  new_prev_fraction_1_x = prev_vol_flux_1_x / remaining_vol_flux_1;
	  new_prev_fraction_1_y = prev_vol_flux_1_y / remaining_vol_flux_1;
	  volume_fluxes_1_x[next_face_x] = new_next_fraction_1_x * new_remaining_vol_flux_1;
	  volume_fluxes_1_x[prev_face_x] = new_prev_fraction_1_x * new_remaining_vol_flux_1;
	  volume_fluxes_1_y[prev_face_y] = new_prev_fraction_1_y * new_remaining_vol_flux_1;
	  bool_check_fluxes_x[next_face_x] = 1.0;
	  bool_check_fluxes_x[prev_face_x] = 1.0;
	  bool_check_fluxes_y[prev_face_y] = 1.0;

	} else {
	  
	  volume_fluxes_1_x[prev_face_x] = prev_fraction_1_x * new_total_vol_flux_1;
	  volume_fluxes_1_x[next_face_x] = next_fraction_1_x * new_total_vol_flux_1;
	  volume_fluxes_1_y[prev_face_y] = prev_fraction_1_y * new_total_vol_flux_1;
	  volume_fluxes_1_y[next_face_y] = next_fraction_1_y * new_total_vol_flux_1;
	  bool_check_fluxes_x[prev_face_x] = 1.0;
	  bool_check_fluxes_x[next_face_x] = 1.0;
	  bool_check_fluxes_y[prev_face_y] = 1.0;
	  bool_check_fluxes_y[next_face_y] = 1.0;
	}

	volume_fluxes_2_x[prev_face_x] = volume_fluxes_x[prev_face_x] - volume_fluxes_1_x[prev_face_x];
	volume_fluxes_2_x[next_face_x] = volume_fluxes_x[next_face_x] - volume_fluxes_1_x[next_face_x];
	volume_fluxes_2_y[prev_face_y] = volume_fluxes_y[prev_face_y] - volume_fluxes_1_y[prev_face_y];
	volume_fluxes_2_y[next_face_y] = volume_fluxes_y[next_face_y] - volume_fluxes_1_y[next_face_y];

	//We take the boundaries into account like if they were periodic (ok for wall conditions as well since flux = flux_sym = 0)
	//X min
	if (ix == 0) {
	  index_t face_sym_x = prev_face_x + nx;
	  volume_fluxes_1_x[face_sym_x] = volume_fluxes_1_x[prev_face_x];
	  volume_fluxes_2_x[face_sym_x] = volume_fluxes_2_x[prev_face_x];
	  bool_check_fluxes_x[face_sym_x] = 1.0;
	}	
	//X max
	if (ix == nx - 1) {
	  index_t face_sym_x = next_face_x - nx;
	  volume_fluxes_1_x[face_sym_x] = volume_fluxes_1_x[next_face_x];
	  volume_fluxes_2_x[face_sym_x] = volume_fluxes_2_x[next_face_x];
	  bool_check_fluxes_x[face_sym_x] = 1.0;
	}
	//Y min
	if (iy == 0) {
	  index_t face_sym_y = ny * nx + ix;
	  volume_fluxes_1_y[face_sym_y] = volume_fluxes_1_y[prev_face_y];
	  volume_fluxes_2_y[face_sym_y] = volume_fluxes_2_y[prev_face_y];
	  bool_check_fluxes_y[face_sym_y] = 1.0;
	}	
	//Y max
	if (iy == ny - 1) {
	  index_t face_sym_y = ix;
	  volume_fluxes_1_y[face_sym_y] = volume_fluxes_1_y[next_face_y];
	  volume_fluxes_2_y[face_sym_y] = volume_fluxes_2_y[next_face_y];
	  bool_check_fluxes_y[face_sym_y] = 1.0;
	}
      }

      if ( (vol_proj_2 != 0) && (vol_proj_2  / (dx * dy) < 1.0e-10) ) {

	RealType total_vol_flux_2 = prev_vol_flux_2_x - next_vol_flux_2_x + prev_vol_flux_2_y - next_vol_flux_2_y;
	RealType prev_fraction_2_x = prev_vol_flux_2_x / total_vol_flux_2;
	RealType next_fraction_2_x = next_vol_flux_2_x / total_vol_flux_2;
	RealType prev_fraction_2_y = prev_vol_flux_2_y / total_vol_flux_2;
	RealType next_fraction_2_y = next_vol_flux_2_y / total_vol_flux_2;
	RealType new_total_vol_flux_2 = - vol_2_lag;
	RealType remaining_vol_flux_2, new_remaining_vol_flux_2, new_prev_fraction_2_x, new_next_fraction_2_x, new_prev_fraction_2_y, new_next_fraction_2_y;

	if (prev_bool_check_fluxes_x == 1.0) {

	  if (prev_bool_check_fluxes_y == 1.0) {

	    if (next_bool_check_fluxes_x == 1.0) {

	      volume_fluxes_2_y[next_face_y] = new_total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_x + prev_vol_flux_2_y);
	      bool_check_fluxes_y[next_face_y] = 1.0;

	    } else if (next_bool_check_fluxes_y == 1.0) {
	      
	      volume_fluxes_2_x[next_face_x] = new_total_vol_flux_2 - (prev_vol_flux_2_x + prev_vol_flux_2_y - next_vol_flux_2_y);
	      bool_check_fluxes_x[next_face_x] = 1.0;

	    } else {

	      remaining_vol_flux_2 = total_vol_flux_2 - (prev_vol_flux_2_x + prev_vol_flux_2_y);
	      new_remaining_vol_flux_2 = new_total_vol_flux_2 - (prev_vol_flux_2_x + prev_vol_flux_2_y);
	      new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2; 
	      new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	      volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	      volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	      bool_check_fluxes_x[next_face_x] = 1.0;
	      bool_check_fluxes_y[next_face_y] = 1.0;	      
	    }

	  } else if (next_bool_check_fluxes_y == 1.0) {

	    if (next_bool_check_fluxes_x == 1.0) {

	      volume_fluxes_2_y[prev_face_y] = new_total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_x - next_vol_flux_2_y);
	      bool_check_fluxes_y[prev_face_y] = 1.0;

	    } else {

	      remaining_vol_flux_2 = total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_y);
	      new_remaining_vol_flux_2 = new_total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_y);
	      new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2;
	      new_prev_fraction_2_y = prev_vol_flux_2_y / remaining_vol_flux_2;
	      volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	      volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	      bool_check_fluxes_x[next_face_x] = 1.0;
	      bool_check_fluxes_y[prev_face_y] = 1.0;
	    }

	  } else if (next_bool_check_fluxes_x == 1.0) {

	    remaining_vol_flux_2 = total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_x);
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 - (prev_vol_flux_2_x - next_vol_flux_2_x);
	    new_prev_fraction_2_y = prev_vol_flux_2_x / remaining_vol_flux_2;
	    new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	    volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;

	  } else {

	    remaining_vol_flux_2 = total_vol_flux_2 - prev_vol_flux_2_x;
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 - prev_vol_flux_2_x;
	    new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2;
	    new_prev_fraction_2_y = prev_vol_flux_2_y / remaining_vol_flux_2;
	    new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	    volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	    bool_check_fluxes_x[next_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }

	} else if (prev_bool_check_fluxes_y == 1.0) {

	  if (next_bool_check_fluxes_x == 1.0) {
	    
	    if (next_bool_check_fluxes_y == 1.0) {
	      
	      volume_fluxes_2_x[prev_face_x] = new_total_vol_flux_2 - (- next_vol_flux_2_x + prev_vol_flux_2_y - next_vol_flux_2_y);
	      bool_check_fluxes_x[prev_face_x] = 1.0;

	    } else {
	      
	      remaining_vol_flux_2 = total_vol_flux_2 - (prev_vol_flux_2_y - next_vol_flux_2_x);
	      new_remaining_vol_flux_2 = new_total_vol_flux_2 - (prev_vol_flux_2_y - next_vol_flux_2_x);
	      new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	      new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	      volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	      volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	      bool_check_fluxes_x[prev_face_x] = 1.0;
	      bool_check_fluxes_y[next_face_y] = 1.0;
	    }

	  } else if (next_bool_check_fluxes_y == 1.0) {
	    
	    remaining_vol_flux_2 = total_vol_flux_2 - (prev_vol_flux_2_y - next_vol_flux_2_y);
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 - (prev_vol_flux_2_y - next_vol_flux_2_y);
	    new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	    new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2;
	    volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_x[next_face_x] = 1.0;

	  } else {
	    
	    remaining_vol_flux_2 = total_vol_flux_2 - prev_vol_flux_2_y;
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 - prev_vol_flux_2_y;
	    new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2;
	    new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	    new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	    volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	    bool_check_fluxes_x[next_face_x] = 1.0;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }
	  
	} else if (next_bool_check_fluxes_x == 1.0) {
	  
	  if (next_bool_check_fluxes_y == 1.0) {
	    
	    remaining_vol_flux_2 = total_vol_flux_2 - (- next_vol_flux_2_x - next_vol_flux_2_y);
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 - (- next_vol_flux_2_x - next_vol_flux_2_y);
	    new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	    new_prev_fraction_2_y = prev_vol_flux_2_y / remaining_vol_flux_2;
	    volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;

	  } else {
	    
	    remaining_vol_flux_2 = total_vol_flux_2 + next_vol_flux_2_x;
	    new_remaining_vol_flux_2 = new_total_vol_flux_2 + next_vol_flux_2_x;
	    new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	    new_prev_fraction_2_y = prev_vol_flux_2_y / remaining_vol_flux_2;
	    new_next_fraction_2_y = next_vol_flux_2_y / remaining_vol_flux_2;
	    volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	    volume_fluxes_2_y[next_face_y] = new_next_fraction_2_y * new_remaining_vol_flux_2;
	    bool_check_fluxes_x[prev_face_x] = 1.0;
	    bool_check_fluxes_y[prev_face_y] = 1.0;
	    bool_check_fluxes_y[next_face_y] = 1.0;
	  }
	  
	} else if (next_bool_check_fluxes_y == 1.0) {
	  
	  remaining_vol_flux_2 = total_vol_flux_2 + next_vol_flux_2_y;
	  new_remaining_vol_flux_2 = new_total_vol_flux_2 + next_vol_flux_2_y;
	  new_next_fraction_2_x = next_vol_flux_2_x / remaining_vol_flux_2;
	  new_prev_fraction_2_x = prev_vol_flux_2_x / remaining_vol_flux_2;
	  new_prev_fraction_2_y = prev_vol_flux_2_y / remaining_vol_flux_2;
	  volume_fluxes_2_x[next_face_x] = new_next_fraction_2_x * new_remaining_vol_flux_2;
	  volume_fluxes_2_x[prev_face_x] = new_prev_fraction_2_x * new_remaining_vol_flux_2;
	  volume_fluxes_2_y[prev_face_y] = new_prev_fraction_2_y * new_remaining_vol_flux_2;
	  bool_check_fluxes_x[next_face_x] = 1.0;
	  bool_check_fluxes_x[prev_face_x] = 1.0;
	  bool_check_fluxes_y[prev_face_y] = 1.0;

	} else {
	  
	  volume_fluxes_2_x[prev_face_x] = prev_fraction_2_x * new_total_vol_flux_2;
	  volume_fluxes_2_x[next_face_x] = next_fraction_2_x * new_total_vol_flux_2;
	  volume_fluxes_2_y[prev_face_y] = prev_fraction_2_y * new_total_vol_flux_2;
	  volume_fluxes_2_y[next_face_y] = next_fraction_2_y * new_total_vol_flux_2;
	  bool_check_fluxes_x[prev_face_x] = 1.0;
	  bool_check_fluxes_x[next_face_x] = 1.0;
	  bool_check_fluxes_y[prev_face_y] = 1.0;
	  bool_check_fluxes_y[next_face_y] = 1.0;
	}

	volume_fluxes_1_x[prev_face_x] = volume_fluxes_x[prev_face_x] - volume_fluxes_2_x[prev_face_x];
	volume_fluxes_1_x[next_face_x] = volume_fluxes_x[next_face_x] - volume_fluxes_2_x[next_face_x];
	volume_fluxes_1_y[prev_face_y] = volume_fluxes_y[prev_face_y] - volume_fluxes_2_y[prev_face_y];
	volume_fluxes_1_y[next_face_y] = volume_fluxes_y[next_face_y] - volume_fluxes_2_y[next_face_y];

	//We take the boundaries into account just like if they were periodic (ok for wall conditions as well since flux = flux_sym = 0)
	//X min
	if (ix == 0) {
	  index_t face_sym_x = prev_face_x + nx;
	  volume_fluxes_1_x[face_sym_x] = volume_fluxes_1_x[prev_face_x];
	  volume_fluxes_2_x[face_sym_x] = volume_fluxes_2_x[prev_face_x];
	  bool_check_fluxes_x[face_sym_x] = 1.0;
	}	
	//X max
	if (ix == nx - 1) {
	  index_t face_sym_x = next_face_x - nx;
	  volume_fluxes_1_x[face_sym_x] = volume_fluxes_1_x[next_face_x];
	  volume_fluxes_2_x[face_sym_x] = volume_fluxes_2_x[next_face_x];
	  bool_check_fluxes_x[face_sym_x] = 1.0;
	}
	//Y min
	if (iy == 0) {
	  index_t face_sym_y = ny * nx + ix;
	  volume_fluxes_1_y[face_sym_y] = volume_fluxes_1_y[prev_face_y];
	  volume_fluxes_2_y[face_sym_y] = volume_fluxes_2_y[prev_face_y];
	  bool_check_fluxes_y[face_sym_y] = 1.0;
	}	
	//Y max
	if (iy == ny - 1) {
	  index_t face_sym_y = ix;
	  volume_fluxes_1_y[face_sym_y] = volume_fluxes_1_y[next_face_y];
	  volume_fluxes_2_y[face_sym_y] = volume_fluxes_2_y[next_face_y];
	  bool_check_fluxes_y[face_sym_y] = 1.0;
	}
      }      

    }
  }
} // end CheckDirectionalVolumeFluxesXY



void CheckDirectionalMassFluxXY(//in
				index_t nx, 
				index_t ny, 
				RealType dx,
				RealType dy,
				const RealType* RESTRICT cell_mass_1,
				const RealType* RESTRICT cell_mass_2,
				const RealType* RESTRICT volume_fluxes_x,
				const RealType* RESTRICT volume_fluxes_y,
				const RealType* RESTRICT cell_volumic_fraction,
				const RealType* RESTRICT volume_fluxes_1_x,
				const RealType* RESTRICT volume_fluxes_1_y,
				const RealType* RESTRICT volume_fluxes_2_x,
				const RealType* RESTRICT volume_fluxes_2_y,
				// out
				RealType* RESTRICT mass_flux_1_x,
				RealType* RESTRICT mass_flux_1_y,
				RealType* RESTRICT mass_flux_2_x,
				RealType* RESTRICT mass_flux_2_y) {
				   
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_lag = dx * dy + volume_fluxes_x[prev_face_x] - volume_fluxes_x[next_face_x] + volume_fluxes_y[prev_face_y] - volume_fluxes_y[next_face_y];
      const RealType vol_1_lag = vol_fraction_1 * vol_lag;
      const RealType vol_2_lag = vol_fraction_2 * vol_lag;

      const RealType prev_vol_flux_1_x = volume_fluxes_1_x[prev_face_x];
      const RealType next_vol_flux_1_x = volume_fluxes_1_x[next_face_x];
      const RealType prev_vol_flux_2_x = volume_fluxes_2_x[prev_face_x];
      const RealType next_vol_flux_2_x = volume_fluxes_2_x[next_face_x];  
      const RealType prev_vol_flux_1_y = volume_fluxes_1_y[prev_face_y];
      const RealType next_vol_flux_1_y = volume_fluxes_1_y[next_face_y];
      const RealType prev_vol_flux_2_y = volume_fluxes_2_y[prev_face_y];
      const RealType next_vol_flux_2_y = volume_fluxes_2_y[next_face_y]; 

      if (prev_vol_flux_1_x == - vol_1_lag) {
	mass_flux_1_x[prev_face_x] = - cell_mass_1[cell_ooo];
      }

      if (next_vol_flux_1_x == vol_1_lag) {
	mass_flux_1_x[next_face_x] = cell_mass_1[cell_ooo];
      }
     
      if (prev_vol_flux_2_x == - vol_2_lag) {
	mass_flux_2_x[prev_face_x] = - cell_mass_2[cell_ooo];
      }

      if (next_vol_flux_2_x == vol_2_lag) {
	mass_flux_2_x[next_face_x] = cell_mass_2[cell_ooo];
      }

      if (prev_vol_flux_1_y == - vol_1_lag) {
	mass_flux_1_y[prev_face_y] = - cell_mass_1[cell_ooo];
      }

      if (next_vol_flux_1_y == vol_1_lag) {
	mass_flux_1_y[next_face_y] = cell_mass_1[cell_ooo];
      }
     
      if (prev_vol_flux_2_y == - vol_2_lag) {
	mass_flux_2_y[prev_face_y] = - cell_mass_2[cell_ooo];
      }

      if (next_vol_flux_2_y == vol_2_lag) {
	mass_flux_2_y[next_face_y] = cell_mass_2[cell_ooo];
      }

    }
  }
} // end CheckDirectionalMassFluxXY



void ProjectDensityMultimatDirect(// in
				  index_t nx, 
				  index_t ny, 
				  RealType dx,
				  RealType dy,
				  const RealType* RESTRICT volume_fluxes_x,
				  const RealType* RESTRICT volume_fluxes_y,
				  const RealType* RESTRICT in_cell_volumic_fraction,
				  const RealType* RESTRICT volume_fluxes_1_x,			     
				  const RealType* RESTRICT volume_fluxes_1_y,			     
				  const RealType* RESTRICT volume_fluxes_2_x,
				  const RealType* RESTRICT volume_fluxes_2_y,
				  const RealType* RESTRICT in_rho_1,
				  const RealType* RESTRICT in_rho_2,
				  const RealType* RESTRICT in_cell_mass,
				  const RealType* RESTRICT mass_flux_x,
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

      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);

      const RealType out_mass_tot = in_cell_mass[cell_ooo] + mass_flux_x[prev_face_x] - mass_flux_x[next_face_x] + mass_flux_y[prev_face_y] - mass_flux_y[next_face_y];
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

      const RealType lagrangian_volume = dx * dy + volume_fluxes_x[prev_face_x] - volume_fluxes_x[next_face_x] + volume_fluxes_y[prev_face_y] - volume_fluxes_y[next_face_y];
      const RealType volume_flux_1_prev_x = volume_fluxes_1_x[prev_face_x];
      const RealType volume_flux_1_next_x = volume_fluxes_1_x[next_face_x];
      const RealType volume_flux_2_prev_x = volume_fluxes_2_x[prev_face_x];
      const RealType volume_flux_2_next_x = volume_fluxes_2_x[next_face_x];
      const RealType volume_flux_1_prev_y = volume_fluxes_1_y[prev_face_y];
      const RealType volume_flux_1_next_y = volume_fluxes_1_y[next_face_y];
      const RealType volume_flux_2_prev_y = volume_fluxes_2_y[prev_face_y];
      const RealType volume_flux_2_next_y = volume_fluxes_2_y[next_face_y];

      RealType out_volume_11 = vol_fraction_1 * lagrangian_volume + volume_flux_1_prev_x - volume_flux_1_next_x + volume_flux_1_prev_y - volume_flux_1_next_y;
      RealType out_volume_22 = vol_fraction_2 * lagrangian_volume + volume_flux_2_prev_x - volume_flux_2_next_x + volume_flux_2_prev_y - volume_flux_2_next_y;
      RealType out_vol_fraction_1 =  out_volume_11 / (out_volume_11 + out_volume_22);
      RealType out_vol_fraction_2 =  out_volume_22 / (out_volume_11 + out_volume_22);      

      RealType out_volume_1 = out_vol_fraction_1 * dx * dy;
      RealType out_volume_2 = out_vol_fraction_2 * dx * dy;

      RealType out_density_1, out_density_2;
 
      if ((out_vol_fraction_1 < 1.0e-12) || (out_mass_1 == 0.0)) {
	
	out_vol_fraction_1 = 0.0;
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_cell_mass[cell_ooo];
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if ((out_vol_fraction_2 < 1.0e-12) || (out_mass_2 == 0.0)) {
	
	out_vol_fraction_1 = 1.0;
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

} // end ProjectDensityMultimatDirect




  // DirectProjectionCornerFluxes
  // ----------------------------


void ComputeDirectionalLagrangianQuantitiesCornerFluxesMultimatXY(index_t nx, 
								  index_t ny, 
								  RealType dt,
								  RealType dx,
								  RealType dy,
								  const RealType* RESTRICT in_vx,
								  const RealType* RESTRICT in_vy,
								  const RealType* RESTRICT cell_mass_1,
								  const RealType* RESTRICT cell_mass_2,
								  const RealType* RESTRICT cell_volumic_fraction,
								  RealType* RESTRICT volume_fluxes,
								  RealType* RESTRICT volume_fluxes_y,
								  index_t* RESTRICT sign_x_corner_fluxes,
								  index_t* RESTRICT sign_y_corner_fluxes,
								  RealType* RESTRICT directional_lagrangian_volume,
								  RealType* RESTRICT directional_lagrangian_density_1,
								  RealType* RESTRICT directional_lagrangian_density_2) {

#pragma omp parallel for
  for (index_t iy = 0; iy < ny + 1; ++iy) {
    for (index_t ix = 0; ix < nx + 1; ++ix) {

      const index_t node_ooo = (nx + 1) * iy + ix;
      sign_x_corner_fluxes[node_ooo] = sign(in_vx[node_ooo]);
      sign_y_corner_fluxes[node_ooo] = sign(in_vy[node_ooo]);
    }
  }

#pragma omp parallel 
  {

    // Faces X
#pragma omp  for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianXpart1");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx + 1; ++ix) {
      
	const index_t face_ooo = ((nx + 1) * iy) + ix;
      
	const index_t prev_node = FaceXNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceXNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];                   
	const RealType vx_next = in_vx[next_node];                   
	const RealType vy_prev = in_vy[prev_node];                   
	const RealType vy_next = in_vy[next_node];

	RealType x = ix * dx;
	RealType ymin = iy * dy;
	RealType ymax = (iy + 1) * dy;

	RealType ymin_lag = ymin + vy_prev * dt;
	RealType ymax_lag = ymax + vy_next * dt;
	RealType xmin_lag = x + vx_prev * dt;
	RealType xmax_lag = x + vx_next * dt;

	// We compute x_ymin, x_ymax, ymin, ymax : the coordinates of the diagonal points of the volume flux

	RealType x_ymin, x_ymax;

	ymin = std::max(ymin, ymin_lag);
	ymax = std::min(ymax, ymax_lag);

	//Computation of coefficients of the line joining the two lagrangian vertices
	RealType a = ymax_lag - ymin_lag;
	RealType b = - (xmax_lag - xmin_lag);
	RealType c = - xmin_lag * a - ymin_lag * b;

	//Computation of x_ymin, x_ymax
	if (b == 0.0) {
	  x_ymin = xmin_lag;
	  x_ymax = xmin_lag;
	} else {
	  x_ymin = - (b * ymin + c) / a;
	  x_ymax = - (b * ymax + c) / a;
	}

	const RealType half = 0.5;
	RealType volume_flux = (ymax - ymin) * (half * (x_ymin + x_ymax) - x);
	if (fabs(volume_flux) / (dx * dy) < epsilon_face) {
	  volume_flux = 0.0;
	}

	volume_fluxes[face_ooo] = volume_flux;	

      }
      //likwid_markerStopRegion("directionnnal_lagrangianXpart1");
    }  

    // Faces Y    
#pragma omp for
    for (index_t iy = 0; iy < ny + 1; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianYpart1");
       //#pragma unroll (UnrollFactor)
#pragma simd
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t face_ooo = (nx * iy) + ix;
      
	const index_t prev_node = FaceYNodeM1(face_ooo, iy, nx);
	const index_t next_node = FaceYNodeP1(face_ooo, iy, nx);

	const RealType vx_prev = in_vx[prev_node];                   
	const RealType vx_next = in_vx[next_node];                   
	const RealType vy_prev = in_vy[prev_node];                   
	const RealType vy_next = in_vy[next_node];

	RealType y = iy * dy;
	RealType xmin = ix * dx;
	RealType xmax = (ix + 1) * dx;

	RealType xmin_lag = xmin + vx_prev * dt;
	RealType xmax_lag = xmax + vx_next * dt;
	RealType ymin_lag = y + vy_prev * dt;
	RealType ymax_lag = y + vy_next * dt;

	// We compute xmin, xmax, y_xmin, y_xmax : the coordinates of the diagonal points of the volume flux

	RealType y_xmin, y_xmax;

	xmin = std::max(xmin, xmin_lag);
	xmax = std::min(xmax, xmax_lag);

	//Computation of coefficients of the line joining the two lagrangian vertices
	RealType a = ymax_lag - ymin_lag;
	RealType b = - (xmax_lag - xmin_lag);
	RealType c = - xmin_lag * a - ymin_lag * b;

	//Computation of y_xmin, y_xmax
	if (a == 0.0) {
	  y_xmin = ymin_lag;
	  y_xmax = ymin_lag;
	} else {
	  y_xmin = - (a * xmin + c) / b;
	  y_xmax = - (a * xmax + c) / b;
	}

	const RealType half = 0.5;
	RealType volume_flux = (xmax - xmin) * (half * (y_xmin + y_xmax) - y);
	if (fabs(volume_flux) / (dx * dy) < epsilon_face) {
	  volume_flux = 0.0;
	}
      
	volume_fluxes_y[face_ooo] = volume_flux;     // 1 Store

     }
      //likwid_markerStopRegion("directionnnal_lagrangianYpart1");
    }  
     
#pragma omp for
    for (index_t iy = 0; iy < ny; ++iy) {
      //likwid_markerStartRegion("directionnnal_lagrangianXpart2");
#pragma unroll (UnrollFactor)
      for (index_t ix = 0; ix < nx; ++ix) {

	const index_t cell_ooo = (nx * iy) + ix;

	const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
	const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);

	const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

	const RealType volume_flux_prev_x = volume_fluxes[prev_face_x];    //1 Load
	const RealType volume_flux_next_x = volume_fluxes[next_face_x];     // 0 load (in cache?)
	const RealType volume_flux_prev_y = volume_fluxes_y[prev_face_y];    //1 Load
	const RealType volume_flux_next_y = volume_fluxes_y[next_face_y];

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	const RealType current_mass_1 = cell_mass_1[cell_ooo];
	const RealType current_mass_2 = cell_mass_2[cell_ooo];

	//Compute algebraic corner fluxes
	RealType volume_flux_prev_prev = 0.0;
	RealType volume_flux_prev_next = 0.0;
	RealType volume_flux_next_prev = 0.0;
	RealType volume_flux_next_next = 0.0;

	if ((in_vx[node_m1m1] > 0.0) && (in_vy[node_m1m1] > 0.0)) {
	  volume_flux_prev_prev = - fabs(in_vx[node_m1m1]) * dt * fabs(in_vy[node_m1m1]) * dt;
	} else if ((in_vx[node_m1m1] < 0.0) && (in_vy[node_m1m1] < 0.0)) {
	  volume_flux_prev_prev = fabs(in_vx[node_m1m1]) * dt * fabs(in_vy[node_m1m1]) * dt;
	}
	if (fabs(volume_flux_prev_prev) / (dx * dy) < epsilon_corner) {
	  volume_flux_prev_prev = 0.0;
	}

	if ((in_vx[node_m1p1] > 0.0) && (in_vy[node_m1p1] < 0.0)) {
	  volume_flux_prev_next = - fabs(in_vx[node_m1p1]) * dt * fabs(in_vy[node_m1p1]) * dt;
	} else if ((in_vx[node_m1p1] < 0.0) && (in_vy[node_m1p1] > 0.0)) {
	  volume_flux_prev_next = fabs(in_vx[node_m1p1]) * dt * fabs(in_vy[node_m1p1]) * dt;
	}
	if (fabs(volume_flux_prev_next) / (dx * dy) < epsilon_corner) {
	  volume_flux_prev_next = 0.0;
	}

	if ((in_vx[node_p1m1] < 0.0) && (in_vy[node_p1m1] > 0.0)) {
	  volume_flux_next_prev = - fabs(in_vx[node_p1m1]) * dt * fabs(in_vy[node_p1m1]) * dt;
	} else if ((in_vx[node_p1m1] > 0.0) && (in_vy[node_p1m1] < 0.0)) {
	  volume_flux_next_prev = fabs(in_vx[node_p1m1]) * dt * fabs(in_vy[node_p1m1]) * dt;
	}
	if (fabs(volume_flux_next_prev) / (dx * dy) < epsilon_corner) {
	  volume_flux_next_prev = 0.0;
	}

	if ((in_vx[node_p1p1] < 0.0) && (in_vy[node_p1p1] < 0.0)) {
	  volume_flux_next_next = - fabs(in_vx[node_p1p1]) * dt * fabs(in_vy[node_p1p1]) * dt;
	} else if ((in_vx[node_p1p1] > 0.0) && (in_vy[node_p1p1] > 0.0)) {
	  volume_flux_next_next = fabs(in_vx[node_p1p1]) * dt * fabs(in_vy[node_p1p1]) * dt;
	}
	if (fabs(volume_flux_next_next) / (dx * dy) < epsilon_corner) {
	  volume_flux_next_next = 0.0;
	}

	const RealType current_directional_lagrangian_volume =         // 1 FMA   1 ADD
	  (dx * dy) - volume_flux_prev_x + volume_flux_next_x - volume_flux_prev_y + volume_flux_next_y
	  + volume_flux_prev_prev + volume_flux_prev_next + volume_flux_next_prev + volume_flux_next_next;

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
} // end ComputeDirectionalLagrangianQuantitiesCornerFluxesMultimatXY



void ComputeDirectionalLagrangianFractionalVolumeFluxesXYCorner(// in
								index_t nx, 
								index_t ny, 
								RealType dx,
								RealType dy,
								RealType dt,
								const RealType* RESTRICT volume_fluxes,
								const RealType* RESTRICT volume_fluxes_y,
								const RealType* RESTRICT directional_lagrangian_volume,
								const RealType* RESTRICT cell_volumic_fraction,
								const RealType* RESTRICT interface_normal_x,
								const RealType* RESTRICT interface_normal_y,
								const RealType* RESTRICT in_vx,
								const RealType* RESTRICT in_vy,
								// out
								RealType* RESTRICT volume_fluxes_1,
								RealType* RESTRICT volume_fluxes_1_y,
								RealType* RESTRICT volume_fluxes_1_corner,
								RealType* RESTRICT volume_fluxes_2,
								RealType* RESTRICT volume_fluxes_2_y,
								RealType* RESTRICT volume_fluxes_2_corner) {
  RealType pointA[nx * ny][2];
  RealType pointB[nx * ny][2];

  // Place Interface on each reconstructed lagragian cell  
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo,iy,nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo,iy,nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo,iy,nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo,iy,nx);

      const RealType normal_x = interface_normal_x[cell_ooo];
      const RealType normal_y = interface_normal_y[cell_ooo];
     
      RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      if ((vol_fraction_1 != 0.0) && (vol_fraction_2 != 0.0)) { 
     
	const RealType xmin_lag = ix * dx + 0.5 * (in_vx[node_m1m1] + in_vx[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (in_vx[node_p1m1] + in_vx[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (in_vy[node_m1m1] + in_vy[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (in_vy[node_p1m1] + in_vy[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA[cell_ooo], pointB[cell_ooo]);

      } else {
	
	pointA[cell_ooo][0] = 1.0e19;
	pointA[cell_ooo][1] = 1.0e19;
	pointB[cell_ooo][0] = 1.0e19;
	pointB[cell_ooo][1] = 1.0e19;
      }
    }
  }

  // Faces X
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

	    const index_t prev_face = CellFaceM1O(cell_ooo, iy, nx);
	    const index_t next_face = CellFaceP1O(cell_ooo, iy, nx);
	    
	    const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	    const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	    const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	    const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);
	    
	    const RealType normal_x = interface_normal_x[cell_ooo];
	    const RealType normal_y = interface_normal_y[cell_ooo];

	    // Algebraic values of left and right displacement
	    RealType ymin = (cell_ooo / nx) * dy;
	    RealType ymax = (cell_ooo / nx + 1) * dy;
	    RealType ymin_lag_prev = ymin + in_vy[node_m1m1] * dt;
	    RealType ymax_lag_prev = ymax + in_vy[node_m1p1] * dt;
	    RealType ymin_lag_next = ymin + in_vy[node_p1m1] * dt;
	    RealType ymax_lag_next = ymax + in_vy[node_p1p1] * dt;
	    RealType ymin_prev = std::max(ymin_lag_prev,ymin);
	    RealType ymin_next = std::max(ymin_lag_next,ymin);
	    RealType ymax_prev = std::min(ymax_lag_prev,ymax);
	    RealType ymax_next = std::min(ymax_lag_next,ymax); 
	    RealType dy_lag_prev = ymax_prev - ymin_prev;
	    RealType dy_lag_next = ymax_next - ymin_next;
	    
	    RealType delta_x_prev = volume_fluxes[prev_face] / dy_lag_prev;
	    RealType delta_x_next = volume_fluxes[next_face] / dy_lag_next;
	  
	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = ix * dx;
	    RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * x + (1 - sign(vol_flux)) * (x + delta_x_prev) );
	    RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * (x + delta_x_next) + (1 - sign(vol_flux)) * x );
	    RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * ymin_next + (1 - sign(vol_flux)) * ymin_prev );
	    RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * ymax_next + (1 - sign(vol_flux)) * ymax_prev );

	    IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[cell_ooo][0] + pointB[cell_ooo][0]) / 2.;
	    RealType yI = (pointA[cell_ooo][1] + pointB[cell_ooo][1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      volume_fluxes_1[face_ooo] = vol_flux;
	      volume_fluxes_2[face_ooo] = 0.0;
	    } else {
	      volume_fluxes_1[face_ooo] = 0.0;
	      volume_fluxes_2[face_ooo] = vol_flux;
	    }	      

	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
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
  } //end Faces X


  // Faces Y
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

	    const index_t prev_face = CellFaceOM1(cell_ooo, iy, nx);
	    const index_t next_face = CellFaceOP1(cell_ooo, iy, nx);

	    const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	    const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	    const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	    const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

	    const RealType normal_x = interface_normal_x[cell_ooo];
	    const RealType normal_y = interface_normal_y[cell_ooo];

	    // Algebraic values of up and down displacement
	    RealType xmin = (cell_ooo - cell_ooo / nx * nx) * dx;
	    RealType xmax = (cell_ooo - cell_ooo / nx * nx + 1) * dx;
	    RealType xmin_lag_prev = xmin + in_vx[node_m1m1] * dt;
	    RealType xmax_lag_prev = xmax + in_vx[node_p1m1] * dt;
	    RealType xmin_lag_next = xmin + in_vx[node_m1p1] * dt;
	    RealType xmax_lag_next = xmax + in_vx[node_p1p1] * dt;
	    RealType xmin_prev = std::max(xmin_lag_prev,xmin);
	    RealType xmin_next = std::max(xmin_lag_next,xmin);
	    RealType xmax_prev = std::min(xmax_lag_prev,xmax);
	    RealType xmax_next = std::min(xmax_lag_next,xmax); 
	    RealType dx_lag_prev = xmax_prev - xmin_prev;
	    RealType dx_lag_next = xmax_next - xmin_next;

	    RealType delta_y_prev = volume_fluxes_y[prev_face] / dx_lag_prev;
	    RealType delta_y_next = volume_fluxes_y[next_face] / dx_lag_next;

	    // Compute the vertices of lagrangian dVol at face_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType y = iy * dy;
	    RealType ymin_d = 0.5 * ( (1 + sign(vol_flux)) * y + (1 - sign(vol_flux)) * (y + delta_y_prev) );
	    RealType ymax_d = 0.5 * ( (1 + sign(vol_flux)) * (y + delta_y_next) + (1 - sign(vol_flux)) * y );
	    RealType xmin_d = 0.5 * ( (1 + sign(vol_flux)) * xmin_next + (1 - sign(vol_flux)) * xmin_prev );
	    RealType xmax_d = 0.5 * ( (1 + sign(vol_flux)) * xmax_next + (1 - sign(vol_flux)) * xmax_prev );

	    IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[cell_ooo][0] + pointB[cell_ooo][0]) / 2.;
	    RealType yI = (pointA[cell_ooo][1] + pointB[cell_ooo][1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      volume_fluxes_1_y[face_ooo] = vol_flux;
	      volume_fluxes_2_y[face_ooo] = 0.0;
	    } else {
	      volume_fluxes_1_y[face_ooo] = 0.0;
	      volume_fluxes_2_y[face_ooo] = vol_flux;
	    }

	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
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
  } //end Faces Y 


  // Corners
  #pragma omp parallel for
  for (index_t iy = 1; iy < ny; ++iy) {
    for (index_t ix = 1; ix < nx; ++ix) {

      const index_t node_ooo = (nx + 1) * iy + ix;

      const RealType delta_x = in_vx[node_ooo] * dt;
      const RealType delta_y = in_vy[node_ooo] * dt;
      const RealType vol_flux = fabs(delta_x) * fabs(delta_y);

      if (vol_flux / (dx * dy) < epsilon_corner) {

	volume_fluxes_1_corner[node_ooo] = 0.0;
	volume_fluxes_2_corner[node_ooo] = 0.0;
	
      } else {
	
	index_t cell_ooo = NodeCellM1M1(node_ooo,iy,nx);
	
	if (delta_x  < 0.0) {
	  if (delta_y > 0.0) {
	    cell_ooo = NodeCellP1M1(node_ooo, iy, nx);
	  } else {
	    cell_ooo = NodeCellP1P1(node_ooo, iy, nx);
	  }
	} else if (delta_y < 0.0) {
	  cell_ooo = NodeCellM1P1(node_ooo, iy, nx);
	}

	const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
	const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
	const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
	const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

	RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
	RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

	// If only one material in cell_ooo
	if ( (vol_fraction_1 == 0.0) || (vol_fraction_2 == 0.0) ) {

	  volume_fluxes_1_corner[node_ooo] = vol_fraction_1 * vol_flux;
	  volume_fluxes_2_corner[node_ooo] = vol_fraction_2 * vol_flux;

	} else { 

	  const RealType normal_x = interface_normal_x[cell_ooo];
	  const RealType normal_y = interface_normal_y[cell_ooo];

	  // Compute the vertices of lagrangian dVol at node_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	  RealType pointA_d[2];
	  RealType pointB_d[2];
	  RealType x = ix * dx;
	  RealType x_lag = x + delta_x;
	  RealType y = iy * dy;
	  RealType y_lag = y + delta_y;
	  RealType xmin_d = std::min(x,x_lag);
	  RealType xmax_d = std::max(x,x_lag);
	  RealType ymin_d = std::min(y,y_lag);
	  RealType ymax_d = std::max(y,y_lag);

	  IntersectLineRectangle(pointA[cell_ooo], pointB[cell_ooo], xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);
	    
	  // Default value : only one material in dVol
	  RealType x_inter =  IntersectLineHorizontalLine(pointA[cell_ooo], pointB[cell_ooo], y_lag);

	  // If the interface is not horinzontal
	  if (x_inter < 1.0e19) {
	    // If the corner is on the right side of the lagrangian interface
	    if (x_inter < x_lag) {
	      if (normal_x < 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    // Else it is on the left side of the lagrangian interface
	    } else {
	      if (normal_x > 0.0) {
		volume_fluxes_1_corner[node_ooo] = vol_flux;
		volume_fluxes_2_corner[node_ooo] = 0.0;
	      } else { 
		volume_fluxes_1_corner[node_ooo] = 0.0;
		volume_fluxes_2_corner[node_ooo] = vol_flux;
	      }
	    }
	  // If the interface is horizontal
	  } else {
	    if (normal_y * delta_y < 0.0) {
	      volume_fluxes_1_corner[node_ooo] = vol_flux;
	      volume_fluxes_2_corner[node_ooo] = 0.0;
	    } else { 
	      volume_fluxes_1_corner[node_ooo] = 0.0;
	      volume_fluxes_2_corner[node_ooo] = vol_flux;
	    }
	  }	  	    
	    
	  // If the interface intersects dVol, we compute the volumic fraction in dVol
	  if ( pointA_d[0] < 1.0e19 ) {

	      // Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at node_ooo
	      RealType vol_fraction_1_d, vol_fraction_2_d;

	      vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

	      if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 0.0;
		vol_fraction_2_d = 1.0;
		
	      } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 1.0;
		vol_fraction_2_d = 0.0;
		
	      } else {
	      
		RealType vol_fraction_11_d = vol_fraction_1_d;
		RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
	      }

	      volume_fluxes_1_corner[node_ooo] = vol_fraction_1_d * vol_flux;
	      volume_fluxes_2_corner[node_ooo] = vol_fraction_2_d * vol_flux;	    	       
	  }

	}
      }  
    }
  } //end Corners 

} // end ComputeDirectionalLagrangianFractionalVolumeFluxesXYCorner



void CheckLagrangianVolumeFraction(//in
				   index_t nx, 
				   index_t ny, 
				   RealType dx,
				   RealType dy,
				   RealType dt,
				   const RealType* RESTRICT directional_lagrangian_volume,
				   const RealType* RESTRICT cell_volumic_fraction,
				   const RealType* RESTRICT volume_fluxes_x,
				   const RealType* RESTRICT volume_fluxes_y,
				   const RealType* RESTRICT interface_normal_x,
				   const RealType* RESTRICT interface_normal_y,
				   const RealType* RESTRICT predicted_u,
				   const RealType* RESTRICT predicted_v,
				   const index_t* RESTRICT sign_x_corner_fluxes,
				   const index_t* RESTRICT sign_y_corner_fluxes,
				   // out
				   RealType* RESTRICT volume_fluxes_1,
				   RealType* RESTRICT volume_fluxes_1_y,
				   RealType* RESTRICT volume_fluxes_1_corner,
				   RealType* RESTRICT volume_fluxes_2,
				   RealType* RESTRICT volume_fluxes_2_y,
				   RealType* RESTRICT volume_fluxes_2_corner) {
				   
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;

      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      if ((vol_fraction_1 != 0.) && (vol_fraction_2 != 0.)) {

	const RealType vol_1 = vol_fraction_1 * dx * dy;
	const RealType vol_2 = vol_fraction_2 * dx * dy;

	const RealType normal_x = interface_normal_x[cell_ooo];
	const RealType normal_y = interface_normal_y[cell_ooo];

	const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
	const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
	const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
	const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);
      
	const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
	const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
	const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
	const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

	const RealType prev_vol_flux_x = volume_fluxes_x[prev_face_x];
	const RealType next_vol_flux_x = volume_fluxes_x[next_face_x];
	const RealType prev_vol_flux_y = volume_fluxes_y[prev_face_y];
	const RealType next_vol_flux_y = volume_fluxes_y[next_face_y];

	const RealType vol_flux_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * (volume_fluxes_1_corner[node_m1m1] + volume_fluxes_2_corner[node_m1m1]);
	const RealType vol_flux_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * (volume_fluxes_1_corner[node_m1p1] + volume_fluxes_2_corner[node_m1p1]);
	const RealType vol_flux_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * (volume_fluxes_1_corner[node_p1m1] + volume_fluxes_2_corner[node_p1m1]);
	const RealType vol_flux_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * (volume_fluxes_1_corner[node_p1p1] + volume_fluxes_2_corner[node_p1p1]);

	const RealType prev_vol_flux_1_x = volume_fluxes_1[prev_face_x];
	const RealType next_vol_flux_1_x = volume_fluxes_1[next_face_x];
	const RealType prev_vol_flux_2_x = volume_fluxes_2[prev_face_x];
	const RealType next_vol_flux_2_x = volume_fluxes_2[next_face_x];  
	const RealType prev_vol_flux_1_y = volume_fluxes_1_y[prev_face_y];
	const RealType next_vol_flux_1_y = volume_fluxes_1_y[next_face_y];
	const RealType prev_vol_flux_2_y = volume_fluxes_2_y[prev_face_y];
	const RealType next_vol_flux_2_y = volume_fluxes_2_y[next_face_y];

	const RealType vol_flux_1_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_1_corner[node_m1m1];
	const RealType vol_flux_1_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_1_corner[node_m1p1];
	const RealType vol_flux_1_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_1_corner[node_p1m1];
	const RealType vol_flux_1_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_1_corner[node_p1p1];
	const RealType vol_flux_2_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_2_corner[node_m1m1];
	const RealType vol_flux_2_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_2_corner[node_m1p1];
	const RealType vol_flux_2_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_2_corner[node_p1m1];
	const RealType vol_flux_2_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_2_corner[node_p1p1];


	// Placement of the interface on the rectangular approximation of Lagrangian cell
	
	RealType pointA[2], pointB[2];
	const RealType xmin_lag = ix * dx + 0.5 * (predicted_u[node_m1m1] + predicted_u[node_m1p1]) * dt;
	const RealType xmax_lag = (ix + 1) * dx + 0.5 * (predicted_u[node_p1m1] + predicted_u[node_p1p1]) * dt;
	const RealType ymin_lag = iy * dy + 0.5 * (predicted_v[node_m1m1] + predicted_v[node_p1m1]) * dt;
	const RealType ymax_lag = (iy + 1) * dy + 0.5 * (predicted_v[node_p1m1] + predicted_v[node_p1p1]) * dt;
	
	PlaceInterface(xmin_lag, xmax_lag, ymin_lag, ymax_lag, normal_x, normal_y, vol_fraction_1, pointA, pointB);

	// Computation of Vol_1_lag_remap = Lagrangian volume of material 1 considered by the remap

	RealType prev_vol_flux_1_x_lag = prev_vol_flux_1_x;
	RealType prev_vol_flux_1_y_lag = prev_vol_flux_1_y;
	RealType next_vol_flux_1_x_lag = next_vol_flux_1_x;
	RealType next_vol_flux_1_y_lag = next_vol_flux_1_y;
	RealType vol_flux_1_m1m1_lag = vol_flux_1_m1m1;
	RealType vol_flux_1_m1p1_lag = vol_flux_1_m1p1;
	RealType vol_flux_1_p1m1_lag = vol_flux_1_p1m1;
	RealType vol_flux_1_p1p1_lag = vol_flux_1_p1p1;
	
	// X Faces
	// Left
	if (prev_vol_flux_x > 0.) {
	  
	    RealType ymin = iy * dy;
	    RealType ymax = (iy + 1) * dy;
	    RealType ymin_lag_prev = ymin + predicted_v[node_m1m1] * dt;
	    RealType ymax_lag_prev = ymax + predicted_v[node_m1p1] * dt;
	    RealType ymin_prev = std::max(ymin_lag_prev,ymin);
	    RealType ymax_prev = std::min(ymax_lag_prev,ymax);
	    RealType dy_lag_prev = ymax_prev - ymin_prev;
	    
	    RealType delta_x_prev = prev_vol_flux_x / dy_lag_prev;
	  
	    // Compute the vertices of lagrangian dVol at left face (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = ix * dx;
	    RealType xmin_d = 0.5 * ( (1 + sign(prev_vol_flux_x)) * x + (1 - sign(prev_vol_flux_x)) * (x + delta_x_prev) );
	    RealType xmax_d = 0.5 * ( (1 + sign(prev_vol_flux_x)) * (x + delta_x_prev) + (1 - sign(prev_vol_flux_x)) * x );
	    RealType ymin_d = ymin_prev;
	    RealType ymax_d = ymax_prev;

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[0] + pointB[0]) / 2.;
	    RealType yI = (pointA[1] + pointB[1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      prev_vol_flux_1_x_lag = prev_vol_flux_x;
	    }	      

	    // If the interface intersects dVol, we compute the volume fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volume fraction of material 1 in dVol at face_ooo
	      RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		}

		prev_vol_flux_1_x_lag = vol_fraction_1_d * prev_vol_flux_x;
	    }
	}
	  
	// Right
	if (next_vol_flux_x < 0.) {
	  
	    RealType ymin = iy * dy;
	    RealType ymax = (iy + 1) * dy;
	    RealType ymin_lag_next = ymin + predicted_v[node_p1m1] * dt;
	    RealType ymax_lag_next = ymax + predicted_v[node_p1p1] * dt;
	    RealType ymin_next = std::max(ymin_lag_next,ymin);
	    RealType ymax_next = std::min(ymax_lag_next,ymax);
	    RealType dy_lag_next = ymax_next - ymin_next;
	    
	    RealType delta_x_next = next_vol_flux_x / dy_lag_next;
	  
	    // Compute the vertices of lagrangian dVol at left face (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = (ix + 1) * dx;
	    RealType xmin_d = 0.5 * ( (1 + sign(next_vol_flux_x)) * x + (1 - sign(next_vol_flux_x)) * (x + delta_x_next) );
	    RealType xmax_d = 0.5 * ( (1 + sign(next_vol_flux_x)) * (x + delta_x_next) + (1 - sign(next_vol_flux_x)) * x );
	    RealType ymin_d = ymin_next;
	    RealType ymax_d = ymax_next;

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[0] + pointB[0]) / 2.;
	    RealType yI = (pointA[1] + pointB[1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      next_vol_flux_1_x_lag = next_vol_flux_x;
	    }	      

	    // If the interface intersects dVol, we compute the volume fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volume fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		}

		next_vol_flux_1_x_lag = vol_fraction_1_d * next_vol_flux_x;
	    }
	}

	// Y Faces
	// Bottom
	if (prev_vol_flux_y > 0.) {
	  
	    RealType xmin = ix * dx;
	    RealType xmax = (ix + 1) * dx;
	    RealType xmin_lag_prev = xmin + predicted_u[node_m1m1] * dt;
	    RealType xmax_lag_prev = xmax + predicted_u[node_p1m1] * dt;
	    RealType xmin_prev = std::max(xmin_lag_prev,xmin);
	    RealType xmax_prev = std::min(xmax_lag_prev,xmax);
	    RealType dx_lag_prev = xmax_prev - xmin_prev;
	    
	    RealType delta_y_prev = prev_vol_flux_y / dx_lag_prev;
	  
	    // Compute the vertices of lagrangian dVol at left face (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType y = iy * dy;
	    RealType ymin_d = 0.5 * ( (1 + sign(prev_vol_flux_y)) * y + (1 - sign(prev_vol_flux_y)) * (y + delta_y_prev) );
	    RealType ymax_d = 0.5 * ( (1 + sign(prev_vol_flux_y)) * (y + delta_y_prev) + (1 - sign(prev_vol_flux_y)) * y );
	    RealType xmin_d = xmin_prev;
	    RealType xmax_d = xmax_prev;

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[0] + pointB[0]) / 2.;
	    RealType yI = (pointA[1] + pointB[1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      prev_vol_flux_1_y_lag = prev_vol_flux_y;
	    }	      

	    // If the interface intersects dVol, we compute the volume fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volume fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		}

		prev_vol_flux_1_y_lag = vol_fraction_1_d * prev_vol_flux_y;
	    }
	}
	  
	// Top
	if (next_vol_flux_y < 0.) {
	  
	    RealType xmin = ix * dx;
	    RealType xmax = (ix + 1) * dx;
	    RealType xmin_lag_next = xmin + predicted_u[node_m1p1] * dt;
	    RealType xmax_lag_next = xmax + predicted_u[node_p1p1] * dt;
	    RealType xmin_next = std::max(xmin_lag_next,xmin);
	    RealType xmax_next = std::min(xmax_lag_next,xmax);
	    RealType dx_lag_next = xmax_next - xmin_next;
	    
	    RealType delta_y_next = next_vol_flux_y / dx_lag_next;
	  
	    // Compute the vertices of lagrangian dVol at left face (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType y = (iy + 1) * dy;
	    RealType ymin_d = 0.5 * ( (1 + sign(next_vol_flux_y)) * y + (1 - sign(next_vol_flux_y)) * (y + delta_y_next) );
	    RealType ymax_d = 0.5 * ( (1 + sign(next_vol_flux_y)) * (y + delta_y_next) + (1 - sign(next_vol_flux_y)) * y );
	    RealType xmin_d = xmin_next;
	    RealType xmax_d = xmax_next;

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[0] + pointB[0]) / 2.;
	    RealType yI = (pointA[1] + pointB[1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      next_vol_flux_1_y_lag = next_vol_flux_y;
	    }	      

	    // If the interface intersects dVol, we compute the volume fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

		// Compute vol_fraction_1_d the volume fraction of material 1 in dVol at face_ooo
		RealType vol_fraction_1_d, vol_fraction_2_d;

		vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

		if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 0.0;
		  vol_fraction_2_d = 1.0;
		
		} else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		  vol_fraction_1_d = 1.0;
		  vol_fraction_2_d = 0.0;
		
		} else {
	      
		  RealType vol_fraction_11_d = vol_fraction_1_d;
		  RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		  vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		  vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
		  assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
		}

		next_vol_flux_1_y_lag = vol_fraction_1_d * next_vol_flux_y;
	    }
	}

	// Corners
	for (index_t corner = 0; corner < 4; corner++) {

	  index_t node_ooo;
	  RealType vol_flux;
	  if (corner == 0) { 
	    node_ooo = node_m1m1;
	    vol_flux = vol_flux_m1m1;
	  }
	  if (corner == 1) { 
	    node_ooo = node_m1p1;
	    vol_flux = vol_flux_m1p1;
	  }
	  if (corner == 2) { 
	    node_ooo = node_p1m1;
	    vol_flux = vol_flux_p1m1;
	  }
	  if (corner == 3) { 
	    node_ooo = node_p1p1;
	    vol_flux = vol_flux_p1p1;
	  }

	  if (vol_flux > 0.) {

	    const RealType delta_x = predicted_u[node_ooo] * dt;
	    const RealType delta_y = predicted_v[node_ooo] * dt;

	    // Compute the vertices of lagrangian dVol at node_ooo (variables noted _d refer to dVol) and init of variables for intersection interface--dVol
	    RealType pointA_d[2];
	    RealType pointB_d[2];
	    RealType x = (node_ooo - node_ooo / (nx + 1) * (nx + 1)) * dx;
	    RealType x_lag = x + delta_x;
	    RealType y = node_ooo / (nx + 1) * dy;
	    RealType y_lag = y + delta_y;
	    RealType xmin_d = std::min(x,x_lag);
	    RealType xmax_d = std::max(x,x_lag);
	    RealType ymin_d = std::min(y,y_lag);
	    RealType ymax_d = std::max(y,y_lag);

	    IntersectLineRectangle(pointA, pointB, xmin_d, xmax_d, ymin_d, ymax_d, pointA_d, pointB_d);

	    // Default value : only one material in dVol
	    // Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	    RealType xI = (pointA[0] + pointB[0]) / 2.;
	    RealType yI = (pointA[1] + pointB[1]) / 2.;
	    RealType xO = (xmin_d + xmax_d) / 2.;
	    RealType yO = (ymin_d + ymax_d) / 2.;
	    RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	    if (ps > 0.0) {
	      if (corner == 0) { 
		vol_flux_1_m1m1_lag = vol_flux;
	      }
	      if (corner == 1) { 
		vol_flux_1_m1p1_lag = vol_flux;
	      }
	      if (corner == 2) { 
		vol_flux_1_p1m1_lag = vol_flux;
	      }
	      if (corner == 3) { 
		vol_flux_1_p1p1_lag = vol_flux;
	      }	      
	    }
	    
	    // If the interface intersects dVol, we compute the volumic fraction in dVol
	    if ( pointA_d[0] < 1.0e19 ) {

	      // Compute vol_fraction_1_d the volumic fraction of material 1 in dVol at node_ooo
	      RealType vol_fraction_1_d, vol_fraction_2_d;

	      vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin_d, xmax_d, ymin_d, ymax_d, normal_x, normal_y);

	      if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 0.0;
		vol_fraction_2_d = 1.0;
		
	      } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
		vol_fraction_1_d = 1.0;
		vol_fraction_2_d = 0.0;
		
	      } else {
	      
		RealType vol_fraction_11_d = vol_fraction_1_d;
		RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
		vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
		vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);

		assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));		
	      }

	      if (corner == 0) { 
		vol_flux_1_m1m1_lag = vol_fraction_1_d * vol_flux;
	      }
	      if (corner == 1) { 
		vol_flux_1_m1p1_lag = vol_fraction_1_d * vol_flux;
	      }
	      if (corner == 2) { 
		vol_flux_1_p1m1_lag = vol_fraction_1_d * vol_flux;
	      }
	      if (corner == 3) { 
		vol_flux_1_p1p1_lag = vol_fraction_1_d * vol_flux;
	      }	
	    	       
	    }
	  }
	}
	
	// Computation of Vol_1_star = volume of material 1 when computing the intersection between Lagrangian interface and Eulerian cell
	RealType vol_1_star;

	// Intersection between Lagrangian interface and Eulerian cell
	RealType pointA_d[2];
	RealType pointB_d[2];
	RealType xmin = ix * dx;
	RealType xmax = (ix + 1) * dx;
	RealType ymin = iy * dy;
	RealType ymax = (iy + 1) * dy;

	IntersectLineRectangle(pointA, pointB, xmin, xmax, ymin, ymax, pointA_d, pointB_d);

	// Default value : only one material in Vol
	// Compute the scalar product between OI (O center of dVol, I middle of the interface) and the normal to the interface
	RealType xI = (pointA[0] + pointB[0]) / 2.;
	RealType yI = (pointA[1] + pointB[1]) / 2.;
	RealType xO = (xmin + xmax) / 2.;
	RealType yO = (ymin + ymax) / 2.;
	RealType ps = (xI - xO) * normal_x + (yI - yO) * normal_y;	    

	if (ps > 0.0) {
	  vol_1_star = dx * dy;
	} else {
	  vol_1_star = 0.;
	}

	// If the interface intersects dVol, we compute the volume fraction in dVol
	if ( pointA_d[0] < 1.0e19 ) {

	  // Compute vol_fraction_1_d the volume fraction of material 1 in dVol at face_ooo
	  RealType vol_fraction_1_d, vol_fraction_2_d;

	  vol_fraction_1_d = ComputePartialVolume(pointA_d, pointB_d, xmin, xmax, ymin, ymax, normal_x, normal_y);

	  if (vol_fraction_1_d < epsilon_dvol_fraction) {
		  
	    vol_fraction_1_d = 0.0;
	    vol_fraction_2_d = 1.0;
		
	  } else if (1.0 - vol_fraction_1_d < epsilon_dvol_fraction) {
		  
	    vol_fraction_1_d = 1.0;
	    vol_fraction_2_d = 0.0;
		
	  } else {
	      
	    RealType vol_fraction_11_d = vol_fraction_1_d;
	    RealType vol_fraction_22_d = 1.0 - vol_fraction_11_d;
	    vol_fraction_1_d = vol_fraction_11_d / (vol_fraction_11_d + vol_fraction_22_d);
	    vol_fraction_2_d = vol_fraction_22_d / (vol_fraction_11_d + vol_fraction_22_d);
		
	    assert((vol_fraction_1_d > 0.0) && (vol_fraction_2_d > 0.0));
	  }

	  vol_1_star = vol_fraction_1_d * dx * dy;
	}

	const RealType vol_1_lag_remap = vol_1_star - (prev_vol_flux_1_x_lag - next_vol_flux_1_x_lag + prev_vol_flux_1_y_lag - next_vol_flux_1_y_lag + vol_flux_1_m1m1_lag + vol_flux_1_m1p1_lag + vol_flux_1_p1m1_lag + vol_flux_1_p1p1_lag);

	// Real Lagrangian volume
	const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];

	// Error between the real Lagrangian volume and the one considered by the remap
	const RealType err_vol_1_lag = vol_1_lag - vol_1_lag_remap; 


	// If the error is not zero, we adjust the biggest outter volume flux of material 1
	if (err_vol_1_lag != 0) {
	  /*
	  printf("++++++++++++++++++++++++---------------- LET'S CHECK THE ERROR OF THE LAGRANGIAN VOLUME ----------------++++++++++++++++++++++++ \n");
	  printf("In cell %d,  vol_total = %f e-6,  vol_fraction_1 = %f,   vol_1 = %f e-6,   vol_lag_total = %f e-6,  vol_1_lag = %f e-6,  err_vol_1_lag = %f e-6 \n\n\n",cell_ooo,dx*dy*1.0e6,vol_fraction_1,vol_1*1.0e6,directional_lagrangian_volume[cell_ooo]*1.0e6,vol_1_lag*1.0e6,err_vol_1_lag*1.0e6);
	  printf("vol_1_star = %f e-6, prev_vol_flux_1_x_lag = %f, next_vol_flux_1_x_lag = %f, prev_vol_flux_1_y_lag = %f, next_vol_flux_1_y_lag = %f, vol_flux_1_m1m1_lag = %f, vol_flux_1_m1p1_lag = %f, vol_flux_1_p1m1_lag = %f, vol_flux_1_p1p1_lag = %f \n\n", vol_1_star*1.0e6, prev_vol_flux_1_x_lag*1.0e6,next_vol_flux_1_x_lag*1.0e6,prev_vol_flux_1_y_lag*1.0e6, next_vol_flux_1_y_lag*1.0e6, vol_flux_1_m1m1_lag*1.0e6, vol_flux_1_m1p1_lag*1.0e6, vol_flux_1_p1m1_lag*1.0e6, vol_flux_1_p1p1_lag*1.0e6);
	  printf("BEFORE : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n",volume_fluxes_1[prev_face_x]*1.0e6,volume_fluxes_1[next_face_x]*1.0e6,volume_fluxes_1_y[prev_face_y]*1.0e6,volume_fluxes_1_y[next_face_y]*1.0e6,volume_fluxes_1_corner[node_m1m1]*1.0e6,volume_fluxes_1_corner[node_m1p1]*1.0e6,volume_fluxes_1_corner[node_p1m1]*1.0e6,volume_fluxes_1_corner[node_p1p1]*1.0e6);
	  */
	  RealType min_vol_flux_1 = std::min(prev_vol_flux_1_x, std::min(- next_vol_flux_1_x, std::min(prev_vol_flux_1_y, std::min(- next_vol_flux_1_y, std::min(vol_flux_1_m1m1, std::min(vol_flux_1_m1p1, std::min(vol_flux_1_p1m1, vol_flux_1_p1p1))))))); 

	  if (min_vol_flux_1 == 0.0) {
	    min_vol_flux_1 = 1.0e20;
	  }

	  if (min_vol_flux_1 == prev_vol_flux_1_x) {
	    volume_fluxes_1[prev_face_x] -= err_vol_1_lag;
	    volume_fluxes_2[prev_face_x] += err_vol_1_lag;
	    if (ix == 0) {
	      volume_fluxes_1[prev_face_x + nx] = volume_fluxes_1[prev_face_x];
	      volume_fluxes_2[prev_face_x + nx] = volume_fluxes_2[prev_face_x];
	    }
	  }
	  if (min_vol_flux_1 == - next_vol_flux_1_x) {
	    volume_fluxes_1[next_face_x] += err_vol_1_lag;
	    volume_fluxes_2[next_face_x] -= err_vol_1_lag;
	    if (ix == nx) {
	      volume_fluxes_1[next_face_x - nx] = volume_fluxes_1[next_face_x];
	      volume_fluxes_2[next_face_x - nx] = volume_fluxes_2[next_face_x];
	    }
	  }
	  if (min_vol_flux_1 == prev_vol_flux_1_y) {
	    volume_fluxes_1_y[prev_face_y] -= err_vol_1_lag;
	    volume_fluxes_2_y[prev_face_y] += err_vol_1_lag;
	    if (iy == 0) {
	      volume_fluxes_1_y[nx * ny + prev_face_y] = volume_fluxes_1_y[prev_face_y];
	      volume_fluxes_2_y[nx * ny + prev_face_y] = volume_fluxes_2_y[prev_face_y];
	    }
	  }
	  if (min_vol_flux_1 == - next_vol_flux_1_y) {
	    volume_fluxes_1_y[next_face_y] += err_vol_1_lag;
	    volume_fluxes_2_y[next_face_y] -= err_vol_1_lag;
	    if (iy == ny) {
	      volume_fluxes_1_y[next_face_y - nx * ny] = volume_fluxes_1_y[next_face_y];
	      volume_fluxes_2_y[next_face_y - nx * ny] = volume_fluxes_2_y[next_face_y];
	    }
	  }
	  if (min_vol_flux_1 == vol_flux_1_m1m1) {
	    volume_fluxes_1_corner[node_m1m1] = fabs(vol_flux_1_m1m1 - err_vol_1_lag);
	    volume_fluxes_1_corner[node_m1m1] = std::max(volume_fluxes_1_corner[node_m1m1],0.0);
	    volume_fluxes_2_corner[node_m1m1] = fabs(vol_flux_2_m1m1 + err_vol_1_lag);
	    volume_fluxes_2_corner[node_m1m1] = std::max(volume_fluxes_2_corner[node_m1m1],0.0);
	    if (cell_ooo == 0) {
	      volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1m1];
	      volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_m1m1];
	      volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1m1];
	      volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1m1];
	      volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_m1m1];
	      volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1m1];
	    }
	  }
	  if (min_vol_flux_1 == vol_flux_1_p1m1) {
	    volume_fluxes_1_corner[node_p1m1] = fabs(vol_flux_1_p1m1 - err_vol_1_lag);
	    volume_fluxes_1_corner[node_p1m1] = std::max(volume_fluxes_1_corner[node_p1m1],0.0);
	    volume_fluxes_2_corner[node_p1m1] = fabs(vol_flux_2_p1m1 + err_vol_1_lag);
	    volume_fluxes_2_corner[node_p1m1] = std::max(volume_fluxes_2_corner[node_p1m1],0.0);
	    if (cell_ooo == nx - 1) {
	      volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1m1];
	      volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1m1];
	      volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_p1m1];
	      volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1m1];
	      volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1m1];
	      volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_p1m1];
	    }
	  }
	  if (min_vol_flux_1 == vol_flux_1_m1p1) {
	    volume_fluxes_1_corner[node_m1p1] = fabs(vol_flux_1_m1p1 - err_vol_1_lag);
	    volume_fluxes_1_corner[node_m1p1] = std::max(volume_fluxes_1_corner[node_m1p1],0.0);
	    volume_fluxes_2_corner[node_m1p1] = fabs(vol_flux_2_m1p1 + err_vol_1_lag);
	    volume_fluxes_2_corner[node_m1p1] = std::max(volume_fluxes_2_corner[node_m1p1],0.0);
	    if (cell_ooo == (ny - 1) * nx) {
	      volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_m1p1];
	      volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1p1];
	      volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1p1];
	      volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_m1p1];
	      volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1p1];
	      volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1p1];
	    }
	  }
	  if (min_vol_flux_1 == vol_flux_1_p1p1) {
	    volume_fluxes_1_corner[node_p1p1] = fabs(vol_flux_1_p1p1 - err_vol_1_lag);
	    volume_fluxes_1_corner[node_p1p1] = std::max(volume_fluxes_1_corner[node_p1p1],0.0);
	    volume_fluxes_2_corner[node_p1p1] = fabs(vol_flux_2_p1p1 + err_vol_1_lag);
	    volume_fluxes_2_corner[node_p1p1] = std::max(volume_fluxes_2_corner[node_p1p1],0.0);
	    if (cell_ooo == ny * nx - 1) {
	      volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1p1];
	      volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_p1p1];
	      volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1p1];
	      volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1p1];
	      volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_p1p1];
	      volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1p1];
	    }
	  }

	}

	//printf("AFTER : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n\n\n",volume_fluxes_1[prev_face_x]*1.0e6,volume_fluxes_1[next_face_x]*1.0e6,volume_fluxes_1_y[prev_face_y]*1.0e6,volume_fluxes_1_y[next_face_y]*1.0e6,volume_fluxes_1_corner[node_m1m1]*1.0e6,volume_fluxes_1_corner[node_m1p1]*1.0e6,volume_fluxes_1_corner[node_p1m1]*1.0e6,volume_fluxes_1_corner[node_p1p1]*1.0e6);
      }

    }
  }
}



void CheckDirectionalVolumeFluxesXYCorner(//in
					  index_t nx, 
					  index_t ny, 
					  RealType dx,
					  RealType dy,
					  const RealType* RESTRICT directional_lagrangian_volume,
					  const RealType* RESTRICT cell_volumic_fraction,
					  const index_t* RESTRICT sign_x_corner_fluxes,
					  const index_t* RESTRICT sign_y_corner_fluxes,
					  // out
					  RealType* RESTRICT volume_fluxes_1,
					  RealType* RESTRICT volume_fluxes_1_y,
					  RealType* RESTRICT volume_fluxes_1_corner,
					  RealType* RESTRICT volume_fluxes_2,
					  RealType* RESTRICT volume_fluxes_2_y,
					  RealType* RESTRICT volume_fluxes_2_corner) {
				   
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);
      
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);
      
      const RealType vol_fraction_11 = cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType vol_1_lag = vol_fraction_1 * directional_lagrangian_volume[cell_ooo];
      const RealType vol_2_lag = vol_fraction_2 * directional_lagrangian_volume[cell_ooo];

      const RealType prev_vol_flux_1_x = volume_fluxes_1[prev_face_x];
      const RealType next_vol_flux_1_x = volume_fluxes_1[next_face_x];
      const RealType prev_vol_flux_2_x = volume_fluxes_2[prev_face_x];
      const RealType next_vol_flux_2_x = volume_fluxes_2[next_face_x];  
      const RealType prev_vol_flux_1_y = volume_fluxes_1_y[prev_face_y];
      const RealType next_vol_flux_1_y = volume_fluxes_1_y[next_face_y];
      const RealType prev_vol_flux_2_y = volume_fluxes_2_y[prev_face_y];
      const RealType next_vol_flux_2_y = volume_fluxes_2_y[next_face_y];

      const RealType vol_flux_1_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_1_corner[node_m1m1];
      const RealType vol_flux_1_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_1_corner[node_m1p1];
      const RealType vol_flux_1_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_1_corner[node_p1m1];
      const RealType vol_flux_1_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_1_corner[node_p1p1];
      const RealType vol_flux_2_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_2_corner[node_m1m1];
      const RealType vol_flux_2_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_2_corner[node_m1p1];
      const RealType vol_flux_2_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_2_corner[node_p1m1];
      const RealType vol_flux_2_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_2_corner[node_p1p1];

      const RealType vol_proj_1 = vol_1_lag + prev_vol_flux_1_x - next_vol_flux_1_x + prev_vol_flux_1_y - next_vol_flux_1_y + vol_flux_1_m1m1 + vol_flux_1_m1p1 + vol_flux_1_p1m1 + vol_flux_1_p1p1;
      const RealType vol_proj_2 = vol_2_lag + prev_vol_flux_2_x - next_vol_flux_2_x + prev_vol_flux_2_y - next_vol_flux_2_y + vol_flux_2_m1m1 + vol_flux_2_m1p1 + vol_flux_2_p1m1 + vol_flux_2_p1p1;

      if ( (vol_proj_1 != 0) && (vol_proj_1 / (dx * dy) < epsilon_check)) {

	/*
	printf("\n\n\n ----------------- ++++++++++++++++++ LET'S CHECK THE DVOL FOR MAT 1 !!  ----------------- ++++++++++++++++++\n");
	printf("In cell %d,   vol proj 1 = %f e-12,   vol lag 1 = %f \n\n\n",cell_ooo,vol_proj_1*1.0e12,vol_1_lag*1.0e12);
	printf("BEFORE : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n",volume_fluxes_1[prev_face_x]*1.0e12,volume_fluxes_1[next_face_x]*1.0e12,volume_fluxes_1_y[prev_face_y]*1.0e12,volume_fluxes_1_y[next_face_y]*1.0e12,volume_fluxes_1_corner[node_m1m1]*1.0e12,volume_fluxes_1_corner[node_m1p1]*1.0e12,volume_fluxes_1_corner[node_p1m1]*1.0e12,volume_fluxes_1_corner[node_p1p1]*1.0e12);
	*/

	RealType min_vol_flux_1 = std::min(prev_vol_flux_1_x, std::min(- next_vol_flux_1_x, std::min(prev_vol_flux_1_y, std::min(- next_vol_flux_1_y, std::min(vol_flux_1_m1m1, std::min(vol_flux_1_m1p1, std::min(vol_flux_1_p1m1, vol_flux_1_p1p1))))))); 

	if (min_vol_flux_1 == 0.0) {
	  min_vol_flux_1 = std::max(prev_vol_flux_1_x, std::max(- next_vol_flux_1_x, std::max(prev_vol_flux_1_y, std::max(- next_vol_flux_1_y, std::max(vol_flux_1_m1m1, std::max(vol_flux_1_m1p1, std::max(vol_flux_1_p1m1, vol_flux_1_p1p1)))))));
	} 

	if (min_vol_flux_1 == 0.0) {
	  min_vol_flux_1 = 1.0e20;
	}

	if (min_vol_flux_1 == prev_vol_flux_1_x) {
	  volume_fluxes_1[prev_face_x] -= vol_proj_1;
	  volume_fluxes_2[prev_face_x] += vol_proj_1;
	  if (ix == 0) {
	    volume_fluxes_1[prev_face_x + nx] = volume_fluxes_1[prev_face_x];
	    volume_fluxes_2[prev_face_x + nx] = volume_fluxes_2[prev_face_x];
	  }
	}
	if (min_vol_flux_1 == - next_vol_flux_1_x) {
	  volume_fluxes_1[next_face_x] += vol_proj_1;
	  volume_fluxes_2[next_face_x] -= vol_proj_1;
	  if (ix == nx) {
	    volume_fluxes_1[next_face_x - nx] = volume_fluxes_1[next_face_x];
	    volume_fluxes_2[next_face_x - nx] = volume_fluxes_2[next_face_x];
	  }
	}
	if (min_vol_flux_1 == prev_vol_flux_1_y) {
	  volume_fluxes_1_y[prev_face_y] -= vol_proj_1;
	  volume_fluxes_2_y[prev_face_y] += vol_proj_1;
	  if (iy == 0) {
	    volume_fluxes_1_y[nx * ny + prev_face_y] = volume_fluxes_1_y[prev_face_y];
	    volume_fluxes_2_y[nx * ny + prev_face_y] = volume_fluxes_2_y[prev_face_y];
	  }
	}
	if (min_vol_flux_1 == - next_vol_flux_1_y) {
	  volume_fluxes_1_y[next_face_y] += vol_proj_1;
	  volume_fluxes_2_y[next_face_y] -= vol_proj_1;
	  if (iy == ny) {
	    volume_fluxes_1_y[next_face_y - nx * ny] = volume_fluxes_1_y[next_face_y];
	    volume_fluxes_2_y[next_face_y - nx * ny] = volume_fluxes_2_y[next_face_y];
	  }
	}
	if (min_vol_flux_1 == vol_flux_1_m1m1) {
	  volume_fluxes_1_corner[node_m1m1] = fabs(vol_flux_1_m1m1 - vol_proj_1);
	  volume_fluxes_1_corner[node_m1m1] = std::max(volume_fluxes_1_corner[node_m1m1],0.0);
	  volume_fluxes_2_corner[node_m1m1] = fabs(vol_flux_2_m1m1 + vol_proj_1);
	  volume_fluxes_2_corner[node_m1m1] = std::max(volume_fluxes_2_corner[node_m1m1],0.0);
	  if (cell_ooo == 0) {
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_m1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1m1];
	  }
	}
	if (min_vol_flux_1 == vol_flux_1_p1m1) {
	  volume_fluxes_1_corner[node_p1m1] = fabs(vol_flux_1_p1m1 - vol_proj_1);
	  volume_fluxes_1_corner[node_p1m1] = std::max(volume_fluxes_1_corner[node_p1m1],0.0);
	  volume_fluxes_2_corner[node_p1m1] = fabs(vol_flux_2_p1m1 + vol_proj_1);
	  volume_fluxes_2_corner[node_p1m1] = std::max(volume_fluxes_2_corner[node_p1m1],0.0);
	  if (cell_ooo == nx - 1) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_p1m1];
	  }
	}
	if (min_vol_flux_1 == vol_flux_1_m1p1) {
	  volume_fluxes_1_corner[node_m1p1] = fabs(vol_flux_1_m1p1 - vol_proj_1);
	  volume_fluxes_1_corner[node_m1p1] = std::max(volume_fluxes_1_corner[node_m1p1],0.0);
	  volume_fluxes_2_corner[node_m1p1] = fabs(vol_flux_2_m1p1 + vol_proj_1);
	  volume_fluxes_2_corner[node_m1p1] = std::max(volume_fluxes_2_corner[node_m1p1],0.0);
	  if (cell_ooo == (ny - 1) * nx) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_m1p1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1p1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1p1];
	  }
	}
	if (min_vol_flux_1 == vol_flux_1_p1p1) {
	  volume_fluxes_1_corner[node_p1p1] = fabs(vol_flux_1_p1p1 - vol_proj_1);
	  volume_fluxes_1_corner[node_p1p1] = std::max(volume_fluxes_1_corner[node_p1p1],0.0);
	  volume_fluxes_2_corner[node_p1p1] = fabs(vol_flux_2_p1p1 + vol_proj_1);
	  volume_fluxes_2_corner[node_p1p1] = std::max(volume_fluxes_2_corner[node_p1p1],0.0);
	  if (cell_ooo == ny * nx - 1) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1p1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_p1p1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1p1];
	  }
	}
	//printf("AFTER : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n\n\n",volume_fluxes_1[prev_face_x]*1.0e12,volume_fluxes_1[next_face_x]*1.0e12,volume_fluxes_1_y[prev_face_y]*1.0e12,volume_fluxes_1_y[next_face_y]*1.0e12,volume_fluxes_1_corner[node_m1m1]*1.0e12,volume_fluxes_1_corner[node_m1p1]*1.0e12,volume_fluxes_1_corner[node_p1m1]*1.0e12,volume_fluxes_1_corner[node_p1p1]*1.0e12);
     }

      if ( (vol_proj_2 != 0) && (vol_proj_2 / (dx * dy) < epsilon_check) ) {

	/*
	printf("\n\n\n ----------------- ++++++++++++++++++ LET'S CHECK THE DVOL FOR MAT 2 !!  ----------------- ++++++++++++++++++\n");
	printf("In cell %d,   vol proj 2 = %f e-12,  vol lag 2 = %f \n",cell_ooo,vol_proj_2*1.0e12,vol_2_lag*1.0e12);
	printf("BEFORE : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n",volume_fluxes_2[prev_face_x]*1.0e12,volume_fluxes_2[next_face_x]*1.0e12,volume_fluxes_2_y[prev_face_y]*1.0e12,volume_fluxes_2_y[next_face_y]*1.0e12,volume_fluxes_2_corner[node_m1m1]*1.0e12,volume_fluxes_2_corner[node_m1p1]*1.0e12,volume_fluxes_2_corner[node_p1m1]*1.0e12,volume_fluxes_2_corner[node_p1p1]*1.0e12);
	*/

	RealType min_vol_flux_2 = std::min(prev_vol_flux_2_x, std::min(- next_vol_flux_2_x, std::min(prev_vol_flux_2_y, std::min(- next_vol_flux_2_y, std::min(vol_flux_2_m1m1, std::min(vol_flux_2_m1p1, std::min(vol_flux_2_p1m1, vol_flux_2_p1p1))))))); 

	if (min_vol_flux_2 == 0.0) {
	  min_vol_flux_2 = std::max(prev_vol_flux_2_x, std::max(- next_vol_flux_2_x, std::max(prev_vol_flux_2_y, std::max(- next_vol_flux_2_y, std::max(vol_flux_2_m1m1, std::max(vol_flux_2_m1p1, std::max(vol_flux_2_p1m1, vol_flux_2_p1p1)))))));
	}

	if (min_vol_flux_2 == 0.0) {
	  min_vol_flux_2 = 1.0e20;
	}

	if (min_vol_flux_2 == prev_vol_flux_2_x) {
	  volume_fluxes_2[prev_face_x] -= vol_proj_2;
	  volume_fluxes_1[prev_face_x] += vol_proj_2;
	  if (ix == 0) {
	    volume_fluxes_1[prev_face_x + nx] = volume_fluxes_1[prev_face_x];
	    volume_fluxes_2[prev_face_x + nx] = volume_fluxes_2[prev_face_x];
	  }
	}
	if (min_vol_flux_2 == - next_vol_flux_2_x) {
	  volume_fluxes_2[next_face_x] += vol_proj_2;
	  volume_fluxes_1[next_face_x] -= vol_proj_2;
	  if (ix == nx) {
	    volume_fluxes_1[next_face_x - nx] = volume_fluxes_1[next_face_x];
	    volume_fluxes_2[next_face_x - nx] = volume_fluxes_2[next_face_x];
	  }
	}
	if (min_vol_flux_2 == prev_vol_flux_2_y) {
	  volume_fluxes_2_y[prev_face_y] -= vol_proj_2;
	  volume_fluxes_1_y[prev_face_y] += vol_proj_2;
	  if (iy == 0) {
	    volume_fluxes_1_y[nx * ny + prev_face_y] = volume_fluxes_1_y[prev_face_y];
	    volume_fluxes_2_y[nx * ny + prev_face_y] = volume_fluxes_2_y[prev_face_y];
	  }
	}
	if (min_vol_flux_2 == - next_vol_flux_2_y) {
	  volume_fluxes_2_y[next_face_y] += vol_proj_2;
	  volume_fluxes_1_y[next_face_y] -= vol_proj_2;
	  if (iy == ny) {
	    volume_fluxes_1_y[next_face_y - nx * ny] = volume_fluxes_1_y[next_face_y];
	    volume_fluxes_2_y[next_face_y - nx * ny] = volume_fluxes_2_y[next_face_y];
	  }
	}
	if (min_vol_flux_2 == vol_flux_2_m1m1) {
	  volume_fluxes_2_corner[node_m1m1] = fabs(vol_flux_2_m1m1 - vol_proj_2);
	  volume_fluxes_2_corner[node_m1m1] = std::max(volume_fluxes_2_corner[node_m1m1],0.0);
	  volume_fluxes_1_corner[node_m1m1] = fabs(vol_flux_1_m1m1 + vol_proj_2);
	  volume_fluxes_1_corner[node_m1m1] = std::max(volume_fluxes_1_corner[node_m1m1],0.0);
	  if (cell_ooo == 0) {
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1m1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_m1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1m1];
	  }
	}
	if (min_vol_flux_2 == vol_flux_2_p1m1) {
	  volume_fluxes_2_corner[node_p1m1] = fabs(vol_flux_2_p1m1 - vol_proj_2);
	  volume_fluxes_2_corner[node_p1m1] = std::max(volume_fluxes_2_corner[node_p1m1],0.0);
	  volume_fluxes_1_corner[node_p1m1] = fabs(vol_flux_1_p1m1 + vol_proj_2);
	  volume_fluxes_1_corner[node_p1m1] = std::max(volume_fluxes_1_corner[node_p1m1],0.0);
	  if (cell_ooo == nx - 1) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_p1m1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1m1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_p1m1];
	  }
	}
	if (min_vol_flux_2 == vol_flux_2_m1p1) {
	  volume_fluxes_2_corner[node_m1p1] = fabs(vol_flux_2_m1p1 - vol_proj_2);
	  volume_fluxes_2_corner[node_m1p1] = std::max(volume_fluxes_2_corner[node_m1p1],0.0);
	  volume_fluxes_1_corner[node_m1p1] = fabs(vol_flux_1_m1p1 + vol_proj_2);
	  volume_fluxes_1_corner[node_m1p1] = std::max(volume_fluxes_1_corner[node_m1p1],0.0);
	  if (cell_ooo == (ny - 1) * nx) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_1_corner[(nx + 1) * ny + nx] = volume_fluxes_1_corner[node_m1p1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_m1p1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_m1p1];
	    volume_fluxes_2_corner[(nx + 1) * ny + nx] = volume_fluxes_2_corner[node_m1p1];
	  }
	}
	if (min_vol_flux_2 == vol_flux_2_p1p1) {
	  volume_fluxes_2_corner[node_p1p1] = fabs(vol_flux_2_p1p1 - vol_proj_2);
	  volume_fluxes_2_corner[node_p1p1] = std::max(volume_fluxes_2_corner[node_p1p1],0.0);
	  volume_fluxes_1_corner[node_p1p1] = fabs(vol_flux_1_p1p1 + vol_proj_2);
	  volume_fluxes_1_corner[node_p1p1] = std::max(volume_fluxes_1_corner[node_p1p1],0.0);
	  if (cell_ooo == ny * nx - 1) {
	    volume_fluxes_1_corner[0] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_1_corner[nx] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_1_corner[(nx + 1) * ny] = volume_fluxes_1_corner[node_p1p1];
	    volume_fluxes_2_corner[0] = volume_fluxes_2_corner[node_p1p1];
	    volume_fluxes_2_corner[nx] = volume_fluxes_2_corner[node_p1p1];
	    volume_fluxes_2_corner[(nx + 1) * ny] = volume_fluxes_2_corner[node_p1p1];
	  }
	}
	//printf("AFTER : prev_x = %f, next_x = %f, prev_y = %f, next_y = %f, m1m1 = %f, m1p1 = %f, p1m1 = %f, p1p1 = %f\n\n\n",volume_fluxes_2[prev_face_x]*1.0e12,volume_fluxes_2[next_face_x]*1.0e12,volume_fluxes_2_y[prev_face_y]*1.0e12,volume_fluxes_2_y[next_face_y]*1.0e12,volume_fluxes_2_corner[node_m1m1]*1.0e12,volume_fluxes_2_corner[node_m1p1]*1.0e12,volume_fluxes_2_corner[node_p1m1]*1.0e12,volume_fluxes_2_corner[node_p1p1]*1.0e12);
      }      

    }
  }
} // end CheckDirectionalVolumeFluxesXYCorner


void CheckDirectionalMassFluxesXYCorner(//in
					  index_t nx, 
					  index_t ny, 
					  RealType dx,
					  RealType dy,
					  const RealType* RESTRICT directional_lagrangian_density_1,
					  const RealType* RESTRICT directional_lagrangian_density_2,
					  const RealType* RESTRICT in_cell_mass_1,
					  const RealType* RESTRICT in_cell_mass_2,
					  const RealType* RESTRICT in_cell_volumic_fraction,
					  const index_t* RESTRICT sign_x_corner_fluxes,
					  const index_t* RESTRICT sign_y_corner_fluxes,
					  // out
					  RealType* RESTRICT mass_flux_1_x,
					  RealType* RESTRICT mass_flux_1_y,
					  RealType* RESTRICT mass_1_corner_fluxes,
					  RealType* RESTRICT mass_flux_2_x,
					  RealType* RESTRICT mass_flux_2_y,
					  RealType* RESTRICT mass_2_corner_fluxes) {
				   
  #pragma omp parallel for
  for (index_t iy = 0; iy < ny; ++iy) {
    for (index_t ix = 0; ix < nx; ++ix) {

      const index_t cell_ooo = nx * iy + ix;
      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);
      
      const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

      const RealType vol_fraction_11 = in_cell_volumic_fraction[cell_ooo];
      const RealType vol_fraction_22 = 1.0 - vol_fraction_11;
      const RealType vol_fraction_1 = vol_fraction_11 / (vol_fraction_11 + vol_fraction_22);
      const RealType vol_fraction_2 = vol_fraction_22 / (vol_fraction_11 + vol_fraction_22);

      const RealType mass_1_lag = in_cell_mass_1[cell_ooo];
      const RealType mass_2_lag = in_cell_mass_2[cell_ooo];
      const RealType mass_lag = mass_1_lag + mass_2_lag;
      const RealType rho_1_lag = directional_lagrangian_density_1[cell_ooo];
      const RealType rho_2_lag = directional_lagrangian_density_2[cell_ooo];
      const RealType rho_lag = vol_fraction_1 * rho_1_lag + vol_fraction_2 * rho_2_lag;

      const RealType prev_mass_flux_1_x = mass_flux_1_x[prev_face_x];
      const RealType next_mass_flux_1_x = mass_flux_1_x[next_face_x];
      const RealType prev_mass_flux_2_x = mass_flux_2_x[prev_face_x];
      const RealType next_mass_flux_2_x = mass_flux_2_x[next_face_x];  
      const RealType prev_mass_flux_1_y = mass_flux_1_y[prev_face_y];
      const RealType next_mass_flux_1_y = mass_flux_1_y[next_face_y];
      const RealType prev_mass_flux_2_y = mass_flux_2_y[prev_face_y];
      const RealType next_mass_flux_2_y = mass_flux_2_y[next_face_y];

      const RealType mass_flux_1_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * mass_1_corner_fluxes[node_m1m1];
      const RealType mass_flux_1_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * mass_1_corner_fluxes[node_m1p1];
      const RealType mass_flux_1_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * mass_1_corner_fluxes[node_p1m1];
      const RealType mass_flux_1_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * mass_1_corner_fluxes[node_p1p1];
      const RealType mass_flux_2_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * mass_2_corner_fluxes[node_m1m1];
      const RealType mass_flux_2_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * mass_2_corner_fluxes[node_m1p1];
      const RealType mass_flux_2_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * mass_2_corner_fluxes[node_p1m1];
      const RealType mass_flux_2_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * mass_2_corner_fluxes[node_p1p1];

      const RealType mass_proj_1 = mass_1_lag + prev_mass_flux_1_x - next_mass_flux_1_x + prev_mass_flux_1_y - next_mass_flux_1_y + mass_flux_1_m1m1 + mass_flux_1_m1p1 + mass_flux_1_p1m1 + mass_flux_1_p1p1;
      const RealType mass_proj_2 = mass_2_lag + prev_mass_flux_2_x - next_mass_flux_2_x + prev_mass_flux_2_y - next_mass_flux_2_y + mass_flux_2_m1m1 + mass_flux_2_m1p1 + mass_flux_2_p1m1 + mass_flux_2_p1p1;

      if ( (mass_proj_1 != 0) && (mass_proj_1 / mass_lag < (rho_1_lag / rho_lag) * epsilon_check) ) {

	RealType min_mass_flux_1 = std::min(prev_mass_flux_1_x, std::min(- next_mass_flux_1_x, std::min(prev_mass_flux_1_y, std::min(- next_mass_flux_1_y, std::min(mass_flux_1_m1m1, std::min(mass_flux_1_m1p1, std::min(mass_flux_1_p1m1, mass_flux_1_p1p1))))))); 

	if (min_mass_flux_1 == 0.0) {
	  min_mass_flux_1 = std::max(prev_mass_flux_1_x, std::max(- next_mass_flux_1_x, std::max(prev_mass_flux_1_y, std::max(- next_mass_flux_1_y, std::max(mass_flux_1_m1m1, std::max(mass_flux_1_m1p1, std::max(mass_flux_1_p1m1, mass_flux_1_p1p1)))))));
	} 

	if (min_mass_flux_1 == 0.0) {
	  min_mass_flux_1 = 1.0e20;
	}

	if (min_mass_flux_1 == prev_mass_flux_1_x) {
	  mass_flux_1_x[prev_face_x] -= mass_proj_1;
	  mass_flux_2_x[prev_face_x] += mass_proj_1;
	  if (ix == 0) {
	    mass_flux_1_x[prev_face_x + nx] = mass_flux_1_x[prev_face_x];
	    mass_flux_2_x[prev_face_x + nx] = mass_flux_2_x[prev_face_x];
	  }
	}
	if (min_mass_flux_1 == - next_mass_flux_1_x) {
	  mass_flux_1_x[next_face_x] += mass_proj_1;
	  mass_flux_2_x[next_face_x] -= mass_proj_1;
	  if (ix == nx) {
	    mass_flux_1_x[next_face_x - nx] = mass_flux_1_x[next_face_x];
	    mass_flux_2_x[next_face_x - nx] = mass_flux_2_x[next_face_x];
	  }
	}
	if (min_mass_flux_1 == prev_mass_flux_1_y) {
	  mass_flux_1_y[prev_face_y] -= mass_proj_1;
	  mass_flux_2_y[prev_face_y] += mass_proj_1;
	  if (iy == 0) {
	    mass_flux_1_y[nx * ny + prev_face_y] = mass_flux_1_y[prev_face_y];
	    mass_flux_2_y[nx * ny + prev_face_y] = mass_flux_2_y[prev_face_y];
	  }
	}
	if (min_mass_flux_1 == - next_mass_flux_1_y) {
	  mass_flux_1_y[next_face_y] += mass_proj_1;
	  mass_flux_2_y[next_face_y] -= mass_proj_1;
	  if (iy == ny) {
	    mass_flux_1_y[next_face_y - nx * ny] = mass_flux_1_y[next_face_y];
	    mass_flux_2_y[next_face_y - nx * ny] = mass_flux_2_y[next_face_y];
	  }
	}
	if (min_mass_flux_1 == mass_flux_1_m1m1) {
	  mass_1_corner_fluxes[node_m1m1] = fabs(mass_flux_1_m1m1 - mass_proj_1);
	  mass_1_corner_fluxes[node_m1m1] = std::max(mass_1_corner_fluxes[node_m1m1],0.0);
	  mass_2_corner_fluxes[node_m1m1] = fabs(mass_flux_2_m1m1 + mass_proj_1);
	  mass_2_corner_fluxes[node_m1m1] = std::max(mass_2_corner_fluxes[node_m1m1],0.0);
	  if (cell_ooo == 0) {
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_m1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_m1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_m1m1];
	  }
	}
	if (min_mass_flux_1 == mass_flux_1_p1m1) {
	  mass_1_corner_fluxes[node_p1m1] = fabs(mass_flux_1_p1m1 - mass_proj_1);
	  mass_1_corner_fluxes[node_p1m1] = std::max(mass_1_corner_fluxes[node_p1m1],0.0);
	  mass_2_corner_fluxes[node_p1m1] = fabs(mass_flux_2_p1m1 + mass_proj_1);
	  mass_2_corner_fluxes[node_p1m1] = std::max(mass_2_corner_fluxes[node_p1m1],0.0);
	  if (cell_ooo == nx - 1) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_p1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_p1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_p1m1];
	  }
	}
	if (min_mass_flux_1 == mass_flux_1_m1p1) {
	  mass_1_corner_fluxes[node_m1p1] = fabs(mass_flux_1_m1p1 - mass_proj_1);
	  mass_1_corner_fluxes[node_m1p1] = std::max(mass_1_corner_fluxes[node_m1p1],0.0);
	  mass_2_corner_fluxes[node_m1p1] = fabs(mass_flux_2_m1p1 + mass_proj_1);
	  mass_2_corner_fluxes[node_m1p1] = std::max(mass_2_corner_fluxes[node_m1p1],0.0);
	  if (cell_ooo == (ny - 1) * nx) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_m1p1];
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_m1p1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_m1p1];
	  }
	}
	if (min_mass_flux_1 == mass_flux_1_p1p1) {
	  mass_1_corner_fluxes[node_p1p1] = fabs(mass_flux_1_p1p1 - mass_proj_1);
	  mass_1_corner_fluxes[node_p1p1] = std::max(mass_1_corner_fluxes[node_p1p1],0.0);
	  mass_2_corner_fluxes[node_p1p1] = fabs(mass_flux_2_p1p1 + mass_proj_1);
	  mass_2_corner_fluxes[node_p1p1] = std::max(mass_2_corner_fluxes[node_p1p1],0.0);
	  if (cell_ooo == ny * nx - 1) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_p1p1];
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_p1p1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_p1p1];
	  }
	}
     }

      if ( (mass_proj_2 != 0) && (mass_proj_2 / mass_lag < (rho_2_lag / rho_lag) * epsilon_check) ) {

	RealType min_mass_flux_2 = std::min(prev_mass_flux_2_x, std::min(- next_mass_flux_2_x, std::min(prev_mass_flux_2_y, std::min(- next_mass_flux_2_y, std::min(mass_flux_2_m1m1, std::min(mass_flux_2_m1p1, std::min(mass_flux_2_p1m1, mass_flux_2_p1p1))))))); 

	if (min_mass_flux_2 == 0.0) {
	  min_mass_flux_2 = std::max(prev_mass_flux_2_x, std::max(- next_mass_flux_2_x, std::max(prev_mass_flux_2_y, std::max(- next_mass_flux_2_y, std::max(mass_flux_2_m1m1, std::max(mass_flux_2_m1p1, std::max(mass_flux_2_p1m1, mass_flux_2_p1p1)))))));
	}

	if (min_mass_flux_2 == 0.0) {
	  min_mass_flux_2 = 1.0e20;
	}

	if (min_mass_flux_2 == prev_mass_flux_2_x) {
	  mass_flux_2_x[prev_face_x] -= mass_proj_2;
	  mass_flux_1_x[prev_face_x] += mass_proj_2;
	  if (ix == 0) {
	    mass_flux_1_x[prev_face_x + nx] = mass_flux_1_x[prev_face_x];
	    mass_flux_2_x[prev_face_x + nx] = mass_flux_2_x[prev_face_x];
	  }
	}
	if (min_mass_flux_2 == - next_mass_flux_2_x) {
	  mass_flux_2_x[next_face_x] += mass_proj_2;
	  mass_flux_1_x[next_face_x] -= mass_proj_2;
	  if (ix == nx) {
	    mass_flux_1_x[next_face_x - nx] = mass_flux_1_x[next_face_x];
	    mass_flux_2_x[next_face_x - nx] = mass_flux_2_x[next_face_x];
	  }
	}
	if (min_mass_flux_2 == prev_mass_flux_2_y) {
	  mass_flux_2_y[prev_face_y] -= mass_proj_2;
	  mass_flux_1_y[prev_face_y] += mass_proj_2;
	  if (iy == 0) {
	    mass_flux_1_y[nx * ny + prev_face_y] = mass_flux_1_y[prev_face_y];
	    mass_flux_2_y[nx * ny + prev_face_y] = mass_flux_2_y[prev_face_y];
	  }
	}
	if (min_mass_flux_2 == - next_mass_flux_2_y) {
	  mass_flux_2_y[next_face_y] += mass_proj_2;
	  mass_flux_1_y[next_face_y] -= mass_proj_2;
	  if (iy == ny) {
	    mass_flux_1_y[next_face_y - nx * ny] = mass_flux_1_y[next_face_y];
	    mass_flux_2_y[next_face_y - nx * ny] = mass_flux_2_y[next_face_y];
	  }
	}
	if (min_mass_flux_2 == mass_flux_2_m1m1) {
	  mass_2_corner_fluxes[node_m1m1] = fabs(mass_flux_2_m1m1 - mass_proj_2);
	  mass_2_corner_fluxes[node_m1m1] = std::max(mass_2_corner_fluxes[node_m1m1],0.0);
	  mass_1_corner_fluxes[node_m1m1] = fabs(mass_flux_1_m1m1 + mass_proj_2);
	  mass_1_corner_fluxes[node_m1m1] = std::max(mass_1_corner_fluxes[node_m1m1],0.0);
	  if (cell_ooo == 0) {
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_m1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_m1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_m1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_m1m1];
	  }
	}
	if (min_mass_flux_2 == mass_flux_2_p1m1) {
	  mass_2_corner_fluxes[node_p1m1] = fabs(mass_flux_2_p1m1 - mass_proj_2);
	  mass_2_corner_fluxes[node_p1m1] = std::max(mass_2_corner_fluxes[node_p1m1],0.0);
	  mass_1_corner_fluxes[node_p1m1] = fabs(mass_flux_1_p1m1 + mass_proj_2);
	  mass_1_corner_fluxes[node_p1m1] = std::max(mass_1_corner_fluxes[node_p1m1],0.0);
	  if (cell_ooo == nx - 1) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_p1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_p1m1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_p1m1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_p1m1];
	  }
	}
	if (min_mass_flux_2 == mass_flux_2_m1p1) {
	  mass_2_corner_fluxes[node_m1p1] = fabs(mass_flux_2_m1p1 - mass_proj_2);
	  mass_2_corner_fluxes[node_m1p1] = std::max(mass_2_corner_fluxes[node_m1p1],0.0);
	  mass_1_corner_fluxes[node_m1p1] = fabs(mass_flux_1_m1p1 + mass_proj_2);
	  mass_1_corner_fluxes[node_m1p1] = std::max(mass_1_corner_fluxes[node_m1p1],0.0);
	  if (cell_ooo == (ny - 1) * nx) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_m1p1];
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_m1p1];
	    mass_1_corner_fluxes[(nx + 1) * ny + nx] = mass_1_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_m1p1];
	    mass_2_corner_fluxes[(nx + 1) * ny + nx] = mass_2_corner_fluxes[node_m1p1];
	  }
	}
	if (min_mass_flux_2 == mass_flux_2_p1p1) {
	  mass_2_corner_fluxes[node_p1p1] = fabs(mass_flux_2_p1p1 - mass_proj_2);
	  mass_2_corner_fluxes[node_p1p1] = std::max(mass_2_corner_fluxes[node_p1p1],0.0);
	  mass_1_corner_fluxes[node_p1p1] = fabs(mass_flux_1_p1p1 + mass_proj_2);
	  mass_1_corner_fluxes[node_p1p1] = std::max(mass_1_corner_fluxes[node_p1p1],0.0);
	  if (cell_ooo == ny * nx - 1) {
	    mass_1_corner_fluxes[0] = mass_1_corner_fluxes[node_p1p1];
	    mass_1_corner_fluxes[nx] = mass_1_corner_fluxes[node_p1p1];
	    mass_1_corner_fluxes[(nx + 1) * ny] = mass_1_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[0] = mass_2_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[nx] = mass_2_corner_fluxes[node_p1p1];
	    mass_2_corner_fluxes[(nx + 1) * ny] = mass_2_corner_fluxes[node_p1p1];
	  }
	}
      }      

    }
  }
} // end CheckDirectionalMassFluxesXYCorner



void ProjectDensityMultimatCornerFluxes(// in
					index_t nx, 
					index_t ny, 
					RealType dx,
					RealType dy,
					const RealType* RESTRICT directional_lagrangian_volume,
					const RealType* RESTRICT in_cell_volumic_fraction,
					const RealType* RESTRICT volume_fluxes_1_x,			     
					const RealType* RESTRICT volume_fluxes_1_y,			     
					const RealType* RESTRICT volume_fluxes_1_corner,			     
					const RealType* RESTRICT volume_fluxes_2_x,
					const RealType* RESTRICT volume_fluxes_2_y,
					const RealType* RESTRICT volume_fluxes_2_corner,
					const RealType* RESTRICT in_rho_1,
					const RealType* RESTRICT in_rho_2,
					const RealType* RESTRICT in_cell_mass,
					const RealType* RESTRICT mass_flux_x,
					const RealType* RESTRICT mass_flux_y,
					const RealType* RESTRICT mass_corner_fluxes,
					const index_t* RESTRICT sign_x_corner_fluxes,
					const index_t* RESTRICT sign_y_corner_fluxes,
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

      const index_t prev_face_x = CellFaceM1O(cell_ooo, iy, nx);
      const index_t next_face_x = CellFaceP1O(cell_ooo, iy, nx);
      const index_t prev_face_y = CellFaceOM1(cell_ooo, iy, nx);
      const index_t next_face_y = CellFaceOP1(cell_ooo, iy, nx);

      const index_t node_m1m1 = CellNodeM1M1(cell_ooo, iy, nx);
      const index_t node_m1p1 = CellNodeM1P1(cell_ooo, iy, nx);
      const index_t node_p1m1 = CellNodeP1M1(cell_ooo, iy, nx);
      const index_t node_p1p1 = CellNodeP1P1(cell_ooo, iy, nx);

      const RealType mass_flux_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * mass_corner_fluxes[node_m1m1];
      const RealType mass_flux_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * mass_corner_fluxes[node_m1p1];
      const RealType mass_flux_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * mass_corner_fluxes[node_p1m1];
      const RealType mass_flux_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * mass_corner_fluxes[node_p1p1];

      const RealType out_mass_tot = in_cell_mass[cell_ooo] + mass_flux_x[prev_face_x] - mass_flux_x[next_face_x] + mass_flux_y[prev_face_y] - mass_flux_y[next_face_y] + mass_flux_m1m1 + mass_flux_m1p1 + mass_flux_p1m1 + mass_flux_p1p1;
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
      const RealType volume_flux_1_prev_x = volume_fluxes_1_x[prev_face_x];
      const RealType volume_flux_1_next_x = volume_fluxes_1_x[next_face_x];
      const RealType volume_flux_2_prev_x = volume_fluxes_2_x[prev_face_x];
      const RealType volume_flux_2_next_x = volume_fluxes_2_x[next_face_x];
      const RealType volume_flux_1_prev_y = volume_fluxes_1_y[prev_face_y];
      const RealType volume_flux_1_next_y = volume_fluxes_1_y[next_face_y];
      const RealType volume_flux_2_prev_y = volume_fluxes_2_y[prev_face_y];
      const RealType volume_flux_2_next_y = volume_fluxes_2_y[next_face_y];
      const RealType vol_flux_1_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_1_corner[node_m1m1];
      const RealType vol_flux_1_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_1_corner[node_m1p1];
      const RealType vol_flux_1_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_1_corner[node_p1m1];
      const RealType vol_flux_1_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_1_corner[node_p1p1];
      const RealType vol_flux_2_m1m1 = 0.5 * (1 + sign_x_corner_fluxes[node_m1m1] * sign_y_corner_fluxes[node_m1m1]) * sign(sign_x_corner_fluxes[node_m1m1]) * volume_fluxes_2_corner[node_m1m1];
      const RealType vol_flux_2_m1p1 = 0.5 * (1 - sign_x_corner_fluxes[node_m1p1] * sign_y_corner_fluxes[node_m1p1]) * sign(sign_x_corner_fluxes[node_m1p1]) * volume_fluxes_2_corner[node_m1p1];
      const RealType vol_flux_2_p1m1 = - 0.5 * (1 - sign_x_corner_fluxes[node_p1m1] * sign_y_corner_fluxes[node_p1m1]) * sign(sign_x_corner_fluxes[node_p1m1]) * volume_fluxes_2_corner[node_p1m1];
      const RealType vol_flux_2_p1p1 = - 0.5 * (1 + sign_x_corner_fluxes[node_p1p1] * sign_y_corner_fluxes[node_p1p1]) * sign(sign_x_corner_fluxes[node_p1p1]) * volume_fluxes_2_corner[node_p1p1];

      const RealType rho_1_lag = in_rho_1[cell_ooo] * (dx * dy) / lagrangian_volume;
      const RealType rho_2_lag = in_rho_2[cell_ooo] * (dx * dy) / lagrangian_volume;
      const RealType rho_lag = vol_fraction_1 * rho_1_lag + vol_fraction_2 * rho_2_lag;

      RealType out_volume_11 = vol_fraction_1 * lagrangian_volume + volume_flux_1_prev_x - volume_flux_1_next_x + volume_flux_1_prev_y - volume_flux_1_next_y + vol_flux_1_m1m1 + vol_flux_1_m1p1 + vol_flux_1_p1m1 + vol_flux_1_p1p1;
      RealType out_volume_22 = vol_fraction_2 * lagrangian_volume + volume_flux_2_prev_x - volume_flux_2_next_x + volume_flux_2_prev_y - volume_flux_2_next_y + vol_flux_2_m1m1 + vol_flux_2_m1p1 + vol_flux_2_p1m1 + vol_flux_2_p1p1;
      RealType out_vol_fraction_1 =  out_volume_11 / (out_volume_11 + out_volume_22);
      RealType out_vol_fraction_2 =  out_volume_22 / (out_volume_11 + out_volume_22);      

      RealType out_volume_1 = out_vol_fraction_1 * dx * dy;
      RealType out_volume_2 = out_vol_fraction_2 * dx * dy;

      RealType out_density_1, out_density_2;
 
      if ((out_vol_fraction_1 < epsilon_proj) || (out_c1 < (rho_1_lag / rho_lag) * epsilon_proj)) {
	//if ((out_vol_fraction_1 < 1.0e-10) || (out_c1 == 0.0)) {
	
	out_vol_fraction_1 = 0.0;
	out_cell_mass_1[cell_ooo] = 0.0;
	out_cell_mass_2[cell_ooo] = out_mass_tot;
	out_density_1 = 0.0;
	out_density_2 = out_cell_mass_2[cell_ooo] / (dx * dy);
	
      } else if ((out_vol_fraction_2 < epsilon_proj) || (out_c2 < (rho_2_lag / rho_lag) * epsilon_proj)) {
	//} else if ((out_vol_fraction_2 < 1.0e-10) || (out_c2 == 0.0)) {
	
	out_vol_fraction_1 = 1.0;
	out_cell_mass_1[cell_ooo] = out_mass_tot;
	out_cell_mass_2[cell_ooo] = 0.0;
	out_density_1 = out_cell_mass_1[cell_ooo] / (dx * dy);
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

} // end ProjectDensityMultimatCornerFluxes




  // ========================
  //    Functions for Lagrangian phase  
  // ========================


RealType TimeStepSGPCMultimat(int nx,
			      int ny,
			      const RealType dx,
			      const RealType dy,
			      const RealType CFL,
			      const RealType gamma_1,
			      const RealType gamma_2,
			      const RealType pi_1,
			      const RealType pi_2,
			      const RealType* RESTRICT density_1,
			      const RealType* RESTRICT density_2,
			      const RealType* RESTRICT pressure_1,
			      const RealType* RESTRICT pressure_2,
			      const RealType* RESTRICT in_cell_volumic_fraction,
			      const RealType* RESTRICT in_velocity_x,
			      const RealType* RESTRICT in_velocity_y) {
 
  RealType fluid_velocity = 0.0;
  RealType max_velocity = 0.0;
  RealType total_velocity = 0.0;
  RealType speed_of_sound_1 = 0.0;
  RealType speed_of_sound_2 = 0.0;
  
  #pragma omp parallel 
  {
    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny + 1; ++iy) {
      for (int ix = 0; ix < nx + 1; ++ix) {
	fluid_velocity = std::max(fluid_velocity,fabs(in_velocity_x[iy * (nx + 1) + ix]));
      }
    }


    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny + 1; ++iy) {
      for (int ix = 0; ix < nx + 1; ++ix) {
	fluid_velocity =  std::max(fluid_velocity,fabs(in_velocity_y[iy * (nx + 1) + ix]));
      }
    }

    
    //    #pragma omp for reduction(max:max_velocity) nowait
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
	const RealType	p_1_ooo = pressure_1[iy * nx + ix];
	const RealType	p_2_ooo = pressure_2[iy * nx + ix];
	const RealType rho_1_ooo =  density_1[iy * nx + ix];
	const RealType rho_2_ooo =  density_2[iy * nx + ix];
	if (in_cell_volumic_fraction[iy * nx + ix] > 0.01) { 
	  speed_of_sound_1 = SpeedOfSound(gamma_1, rho_1_ooo, p_1_ooo, pi_1);
	}
	if (in_cell_volumic_fraction[iy * nx + ix] < 0.99) { 
	  speed_of_sound_2 = SpeedOfSound(gamma_2, rho_2_ooo, p_2_ooo, pi_2);
	}

        total_velocity = fluid_velocity + std::max(speed_of_sound_1,speed_of_sound_2);
	max_velocity =  std::max(max_velocity,total_velocity);

      }
    }
  }

  RealType dt = CFL * std::min(dx,dy)/max_velocity;
  return dt;
}



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

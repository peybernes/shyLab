#include "geom.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

double Determinant3By3(RealType x1, RealType y1, RealType z1,
		       RealType x2, RealType y2, RealType z2,
		       RealType x3, RealType y3, RealType z3) {
  
  const double result =
    (x1 * y2 * z3) + (y1 * z2 * x3) + (z1 * x2 * y3) -
    (z1 * y2 * x3) - (x1 * z2 * y3) - (y1 * x2 * z3);

  return result;

}


double TetraedronVolume(RealType xa, RealType ya, RealType za,
			RealType xb, RealType yb, RealType zb,
			RealType xc, RealType yc, RealType zc,
			RealType xd, RealType yd, RealType zd) {
  
  const double determinant =
    Determinant3By3(xa - xd, ya - yd, za - zd,
    		    xb - xd, yb - yd, zb - zd,
    		    xc - xd, yc - yd, zc - zd);

  return fabs(determinant) / 6.0;

}

RealType Orient2d(RealType xa, RealType ya, RealType za,
		  RealType xb, RealType yb, RealType zb,
		  RealType xc, RealType yc, RealType zc) {

  return ((xa - xc) * (yb - yc) - (ya - yc) * (xb -xc));

}

// Is the point p inside the triangle abc ? the z component is discarded
bool IsInside(RealType xp, RealType yp, RealType zp,
	      RealType xa, RealType ya, RealType za,
	      RealType xb, RealType yb, RealType zb,
	      RealType xc, RealType yc, RealType zc) {


  // First, check if the point is in the bounding box of the triangle
  // vertices (else, the algorithm is not nearly robust enough)
  const RealType xmin = std::min(std::min(xa, xb), xc);
  const RealType xmax = std::max(std::max(xa, xb), xc);
  const RealType ymin = std::min(std::min(ya, yb), yc);
  const RealType ymax = std::max(std::max(ya, yb), yc);
  const RealType zmin = std::min(std::min(za, zb), zc);
  const RealType zmax = std::max(std::max(za, zb), zc);

  if ((xp < xmin) || (xp > xmax) ||
      (yp < ymin) || (yp > ymax) ||
      (zp < zmin) || (zp > zmax))
    return false;
  
  // Case where the point is in the bounding box. Check if the
  // Triangle vertices are clockwise or counter-clockwise.
  RealType insider = 1.0;
  
  if (Orient2d(xa, ya, za, 
	       xb, yb, zb, 
	       xc, yc, zc) > 0) {

    // counter clockwise.
    insider = Orient2d(xa, ya, za, 
		       xp, yp, zp, 
		       xc, yc, zc);
    insider *= Orient2d(xa, ya, za, 
			xb, yb, zb, 
			xp, yp, zp);
    insider *= Orient2d(xb, yb, zb, 
			xc, yc, zc, 
			xp, yp, zp);
  }
  
  else {
    // clockwise
    insider = Orient2d(xa, ya, za, 
		       xp, yp, zp, 
		       xb, yb, zb);
    insider *= Orient2d(xa, ya, za, 
			xc, yc, zc, 
			xp, yp, zp);
    insider *= Orient2d(xc, yc, zc, 
			xb, yb, zb, 
			xp, yp, zp);

  }

  return (insider > 0.0);
}

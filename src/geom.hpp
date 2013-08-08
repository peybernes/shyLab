#ifndef GEOM_HPP
#define GEOM_HPP

#include <cassert>
#include <iostream>

double TetraedronVolume(RealType xa, RealType ya, RealType za,
			RealType xb, RealType yb, RealType zb,
			RealType xc, RealType yc, RealType zc,
			RealType xd, RealType yd, RealType zd);

double Determinant3By3(RealType x1, RealType y1, RealType z1,
		       RealType x2, RealType y2, RealType z2,
		       RealType x3, RealType y3, RealType z3);
  
// are pa, pb pc in clockwise or counterclockwise order ?
// (Not robust). The z component is discarded
RealType Orient2d(RealType xa, RealType ya, RealType za,
		  RealType xb, RealType yb, RealType zb,
		  RealType xc, RealType yc, RealType zc);

// Is the point p inside the triangle abc ? the z component is discarded
bool IsInside(RealType xp, RealType yp, RealType zp,
	      RealType xa, RealType ya, RealType za,
	      RealType xb, RealType yb, RealType zb,
	      RealType xc, RealType yc, RealType zc);

#endif // GEOM_HPP

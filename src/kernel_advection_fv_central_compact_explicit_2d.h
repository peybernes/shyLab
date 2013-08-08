#ifndef KERNEL_ADVECTION_FD_CENTRAL_COMPACT_EXPLICIT_2D_H
#define KERNEL_ADVECTION_FD_CENTRAL_COMPACT_EXPLICIT_2D_H

#include <assert.h>
#include <math.h>

void AdvectionFdCentralCompact2dNodal(int nx, 
				      int ny, 
				      RealType dt,
				      RealType dx,
				      RealType dy,
				      const RealType* in_value_cell, 
				      const RealType* in_value_node, 
				      const RealType* in_velocity_x_cell,
				      const RealType* in_velocity_y_cell,
				      const RealType* in_velocity_x_node,
				      const RealType* in_velocity_y_node,
				      RealType* out_value_node);

void AdvectionFdCentralCompact2dCentered(int nx, 
					 int ny, 
					 RealType dt,
					 RealType dx,
					 RealType dy,
					 const RealType* in_value_cell, 
					 const RealType* in_value_node, 
					 const RealType* in_velocity_x_cell,
					 const RealType* in_velocity_y_cell,
					 const RealType* in_velocity_x_node,
					 const RealType* in_velocity_y_node,
					 RealType* out_value_cell);



#endif // KERNEL_ADVECTION_FD_CENTRAL_COMPACT_EXPLICIT_2D_H

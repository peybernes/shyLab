#ifndef KERNEL_LAGRANGE_2D_H
#define KERNEL_LAGRANGE_2D_H

typedef int index_t;

#ifdef __cplusplus
extern "C" {
#endif

RealType TimeStep(int nx,
		  int ny,
		  const RealType dx,
		  const RealType dy,
		  const RealType CFL,
		  const RealType gamma,
		  const RealType pi,
		  const RealType* density,
		  const RealType* pressure,
		  const RealType* in_velocity_x,
		  const RealType* in_velocity_y);


void LagrangePressurePredicted(int nx,
			       int ny,
			       RealType dt,
		 	       RealType dx,
		       	       RealType dy,
		       	       RealType gamma,
		       	       RealType pi,
	       		       const RealType* in_mass,
			       const RealType* in_enery,	 
			       const RealType* in_velocity_x,
		       	       const RealType* in_velocity_y,
			       RealType* out_pressure,
	       	      	       RealType* out_predicted_pressure,
       			       RealType* out_pseudo_pressure);

void LagrangePressurePredictedOptimised(int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					RealType gamma,
					RealType pi,
					const RealType* in_mass,
					const RealType* in_enery,	 
					const RealType* in_velocity_x,
					const RealType* in_velocity_y,
					RealType* out_pressure,
					RealType* out_predicted_pressure,
					RealType* out_pseudo_pressure);

void LagrangeVelocityPredicted(int nx,
			    int ny,
			    RealType dt,
			    RealType dx,
			    RealType dy,
			    const RealType* in_mass,
			    const RealType* in_pressure,
			    const RealType* in_pseudo_pressure,
			    const RealType* in_velocity_x,
			    const RealType* in_velocity_y,
			    RealType* out_velocity_x,
			    RealType* out_velocity_y);

void LagrangeCorrection(int nx,
			int ny,
			RealType dt,
			RealType dx,
			RealType dy,
			const RealType* in_mass,
			const RealType* in_enery,
			const RealType* in_pressure,
			const RealType* in_pseudo_pressure,
			const RealType* in_velocity_x,
			const RealType* in_velocity_y,
			RealType* out_enery);

void LagrangeCorrectionOptimised(int nx,
				 int ny,
				 RealType dt,
				 RealType dx,
				 RealType dy,
				 const RealType* in_mass,
				 const RealType* in_enery,
				 const RealType* in_pressure,
				 const RealType* in_pseudo_pressure,
				 const RealType* in_velocity_x,
				 const RealType* in_velocity_y,
				 RealType* out_enery);


void LagrangeVelocityCorrection(int nx,
				int ny,
				const RealType*  in_velocity_x,
				const RealType*  in_velocity_y,
				const RealType*  predicted_velocity_x,
				const RealType*  predicted_velocity_y,
				RealType*  lagrangian_velocity_x,
				RealType*  lagrangian_velocity_y);


void PeriodicBoundaryCopy(int nx,
			  int ny,
			  RealType*  in_velocity_x,
			  RealType*  in_velocity_y);

void RtBoundaryCopy(int nx,
		    int ny,
		    RealType* in_velocity_x,
		    RealType* in_velocity_y);


void PeriodicBoundaryVelocityPrediction(int nx,
					int ny,
					RealType dt,
					RealType dx,
					RealType dy,
					const RealType*  in_mass,
					const RealType*  in_pressure,
					const RealType*  in_pseudo_pressure,
					const RealType*  in_velocity_x,
					const RealType*  in_velocity_y,
					RealType*  out_velocity_x,
					RealType*  out_velocity_y);

void RtBoundaryVelocityPrediction(int nx,
				  int ny,
				  RealType dt,
				  RealType dx,
				  RealType dy,
				  const RealType* in_mass,
				  const RealType* in_pressure,
				  const RealType* in_pseudo_pressure,
				  const RealType* in_velocity_x,
				  const RealType* in_velocity_y,
				  RealType* out_velocity_x,
				  RealType* out_velocity_y);

void WallBoundaryVelocityPrediction(int nx,
				  int ny,
				  RealType dt,
				  RealType dx,
				  RealType dy,
				  const RealType* in_mass,
				  const RealType* in_pressure,
				  const RealType* in_pseudo_pressure,
				  const RealType* in_velocity_x,
				  const RealType* in_velocity_y,
				  RealType* out_velocity_x,
				  RealType* out_velocity_y);

void CopyVelocity(int nx,
		  int ny,
		  const RealType* in_u,
		  RealType* out_u);

#ifdef __cplusplus
}
#endif


#endif // KERNEL_LAGRANGE_2D_H

#ifndef NUMERICAL_CORE_H
#define NUMERICAL_CORE_H

#ifdef __cplusplus
extern "C" {
#endif 
  
  void InterpolateCellToFace(int id_face_begin,
			     int id_face_end,
			     const int* __restrict face_back_cell_id,
			     const int* __restrict face_front_cell_id,
			     const RealType* __restrict in_rho, 
			     const RealType* __restrict in_u,
			     const RealType* __restrict in_v, 
			     const RealType* __restrict in_w,
			     const RealType* __restrict in_e,
			     const RealType* __restrict normal_x, 
			     const RealType* __restrict normal_y,
			     const RealType* __restrict normal_z,
 			     RealType* __restrict rhoL_face,  
			     RealType* __restrict uL_face,
			     RealType* __restrict vL_face,  
			     RealType* __restrict wL_face,
			     RealType* __restrict eL_face,
			     RealType* __restrict rhoR_face,  
			     RealType* __restrict uR_face,
			     RealType* __restrict vR_face,  
			     RealType* __restrict wR_face,
			     RealType* __restrict eR_face);

  void NumericalFluxCellToFace(int id_face_begin,
			       int id_face_end,
			       const int* __restrict face_back_cell_id,
			       const int* __restrict face_front_cell_id,
			       const RealType* __restrict face_areas,
			       const RealType* __restrict in_rho,
			       const RealType* __restrict in_u,
			       const RealType* __restrict in_v,
			       const RealType* __restrict in_w,
			       const RealType* __restrict in_e,
			       const RealType* __restrict normal_x, 
			       const RealType* __restrict normal_y,
			       const RealType* __restrict normal_z,
			       RealType* __restrict flux_rho, 
			       RealType* __restrict flux_rho_u, 
			       RealType* __restrict flux_rho_v,
			       RealType* __restrict flux_rho_w,
			       RealType* __restrict flux_rho_e,
			       RealType* __restrict face_max_eigenvalue);

  void NumericalFluxFaceToFace(int id_face_begin,
			       int id_face_end,
			       const RealType* __restrict face_areas,
			       const RealType* __restrict rhoL_face,
			       const RealType* __restrict uL_face,
			       const RealType* __restrict vL_face,
			       const RealType* __restrict wL_face,
			       const RealType* __restrict eL_face, 
			       const RealType* __restrict hR_face, 
			       const RealType* __restrict uR_face, 
			       const RealType* __restrict vR_face,
			       const RealType* __restrict wR_face,
			       const RealType* __restrict eR_face,
			       const RealType* __restrict normal_x, 
			       const RealType* __restrict normal_y,
			       const RealType* __restrict normal_z,
			       RealType* __restrict flux_rho, 
			       RealType* __restrict flux_rho_u, 
			       RealType* __restrict flux_rho_v,
			       RealType* __restrict flux_rho_w,
			       RealType* __restrict flux_rho_e,
			       RealType* __restrict face_max_eigen);

  void ComputeTimestep(int nb_cells,
		       RealType CFL, 
		       RealType max_timestep,
		       RealType* __restrict cell_timestep,
		       RealType* __restrict out_timestep_ptr);

  void CellResiduals(int id_cell_begin, 
		     int id_cell_end,
		     int nb_faces, 
		     const int* __restrict cell_signed_face1_id,
		     const int* __restrict cell_signed_face2_id,
		     const int* __restrict cell_signed_face3_id,
		     const int* __restrict cell_signed_face4_id,
		     const RealType* __restrict flux_rho,
		     const RealType* __restrict flux_rho_u, 
		     const RealType* __restrict flux_rho_v,
		     const RealType* __restrict flux_rho_w,
		     const RealType* __restrict flux_rho_e,
		     const RealType* __restrict face_max_eigenvalue,
		     const RealType* __restrict normal_x,
		     const RealType* __restrict normal_y,
		     const RealType* __restrict normal_z,
		     const RealType* __restrict cell_volume,
		     RealType* __restrict out_rho, 
		     RealType* __restrict out_rho_u, 
		     RealType* __restrict out_rho_v,
		     RealType* __restrict out_rho_w,
		     RealType* __restrict out_rho_e,
		     RealType* __restrict cell_timestep);

  void EvolveValuesRK1(int nb_cells,
		       RealType EPS, 
		       RealType dt,
		       const RealType* __restrict in_rho,
		       const RealType* __restrict in_u,
		       const RealType* __restrict in_v,
		       const RealType* __restrict in_w,
		       const RealType* __restrict in_e,
		       RealType* __restrict out_rho,
		       RealType* __restrict out_rho_u,
		       RealType* __restrict out_rho_v,
		       RealType* __restrict out_rho_w,
		       RealType* __restrict out_rho_e);

  void CopyValues(int nb_cells,
		  const RealType* __restrict src_rho,
		  const RealType* __restrict src_u,
		  const RealType* __restrict src_v,
		  const RealType* __restrict src_w,
		  const RealType* __restrict src_e,
		  RealType* __restrict dst_rho,
		  RealType* __restrict dst_u,
		  RealType* __restrict dst_v,
		  RealType* __restrict dst_w,
		  RealType* __restrict dst_e);

  void HMPPAllocate(int nb_boundary_faces, int nb_faces, int nb_cells);

  void HMPPRelease();

#ifdef __cplusplus
}
#endif 

#endif // NUMERICAL_CORE_H

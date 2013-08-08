#ifndef MESH_HPP
#define MESH_HPP

#include <iosfwd>
#include <map>
#include <debug/vector>

#include "cell.hpp"
#include "grid.hpp"
#include "mesh_elements.hpp"
#include "mpi_boundary_condition.hpp"
#include "serialize_fwd.hpp"

class VariableStore;

class Mesh : public Grid {
public:
  /// Default constructor.
  Mesh();
  /// Construct mesh from mesh file.
  Mesh(int nb_processes, int process_rank, const std::string& filename);
  /// Read mesh from file (format : gmsh binary).
  void ReadGmsh(int nb_processes, int process_rank, std::istream* is_ptr);
  /// Renumber cells to reduce cache misses. Use the reverse Cuthill 
  // Mc Kee greedy heuristic.
  void ReverseCuthillMcKeeRenumbering();
  /// Write mesh bandwith as an image.
  void WriteMeshBandwith(std::ofstream &);
  /// Compute face-->cells and cells-->faces connectivities.
  void ComputeFaceCellConnectivity();
  void ComputeCellFaceConnectivity();
  virtual int dimension() { return -1; }
  /// Compute geometric quantities.
  virtual void ComputeGeometricQuantities(VariableStore*, 
					  VariableStore*, 
					  VariableStore*);
  /// Locate the cell containing a given location.
  int TriangleIndex(RealType xp, RealType yp, RealType zp);
  virtual void WriteHeaderVTKXmlAscii(std::ofstream* os_ptr);
  virtual void WriteFooterVTKXmlAscii(std::ofstream* os_ptr);
  /// Write mesh in VTK XML ascii format.
  virtual void WriteVTKXmlAscii(std::ofstream* os_ptr);
  /// Write mesh facets in VTK XML ascii format.
  void WriteFacesVTKXmlAscii(std::ofstream* os_ptr);
  /// Write mesh in VTK legacy binary format.
  void WriteVTKBinary(std::ostream &os, bool swap);
  /// Load mesh.
  void Load(ptree&);
  /// Save mesh.
  void Save(ptree&);
  /// Accessor to input filename.
  void ComputeFacesMap();
  void ComputeVerticesMap();
  void ComputeCellsMap();
  /// Uses faces_map.
  void ReorderFaces(VariableStore*);
  /// Uses cells_map.
  void ReorderCells(VariableStore*);
  std::string filename();
  virtual inline int nb_vertices() {return m_nb_vertices;}
  virtual inline int nb_vertices() const {return m_nb_vertices;}
  virtual inline int nb_cells() {return m_nb_cells;}
  virtual int nb_cells() const {return m_nb_cells;}
  inline int nb_faces() {return m_nb_faces;}
  inline int nb_faces() const {return m_nb_faces;}
  inline int nb_interior_faces() {return m_nb_interior_faces;}
  inline int nb_interior_faces() const {return m_nb_interior_faces;}
  inline int nb_boundary_faces() {return m_nb_boundary_faces;}
  inline int nb_boundary_faces() const {return m_nb_boundary_faces;}
  inline int const* face_front_cell_id() const {
    return (m_face_front_cell_id.size() == 0 ? NULL:
	    &(m_face_front_cell_id.at(0)));}
  inline int const* face_back_cell_id() const {
    return (m_face_back_cell_id.size() == 0 ? NULL:
	    &(m_face_back_cell_id.at(0)));}
  inline int const* cell_signed_face1_id() const {
    return (m_cell_signed_face1_id.size() == 0 ? NULL:
	    &(m_cell_signed_face1_id.at(0)));}
  inline int const* cell_signed_face2_id() const {
    return (m_cell_signed_face2_id.size() == 0 ? NULL:
	    &(m_cell_signed_face2_id.at(0)));}
  inline int const* cell_signed_face3_id() const {
    return (m_cell_signed_face3_id.size() == 0 ? NULL:
	    &(m_cell_signed_face3_id.at(0)));}
  inline int const* cell_signed_face4_id() const {
    return (m_cell_signed_face4_id.size() == 0 ? NULL:
	    &(m_cell_signed_face4_id.at(0)));}
  inline std::map<int,int> const& unique_to_local_cell_ids() const {
    return m_unique_to_local_cell_ids;}
  inline std::vector<MpiBoundaryCondition> const&
  mpi_boundary_conditions() const {
    return m_mpi_boundary_conditions;
  }
  inline std::vector<MpiBoundaryCondition>&
  mpi_boundary_conditions() {
    return m_mpi_boundary_conditions;
  }
private:
  std::vector<int> m_vertices_map;
  std::vector<int> m_faces_map;
  std::vector<int> m_cells_map;
  std::vector<RealType> m_vertices_x;
  std::vector<RealType> m_vertices_y;
  std::vector<RealType> m_vertices_z;
  int m_nb_cells;
  int m_nb_vertices;
  int m_nb_faces;
  int m_nb_interior_faces;
  int m_nb_boundary_faces;
  std::string m_filename;
  std::vector<Cell> m_cells;
  std::vector< std::vector<int> > m_cell_cell_connectivity;
  std::vector< std::vector<int> > m_face_vertices_id;
  std::vector<int> m_cell_signed_face1_id;
  std::vector<int> m_cell_signed_face2_id;
  std::vector<int> m_cell_signed_face3_id;
  std::vector<int> m_cell_signed_face4_id;
  std::vector<int> m_face_front_cell_id;
  std::vector<int> m_face_back_cell_id;
  std::map<int,int> m_unique_to_local_cell_ids;
  std::vector<MpiBoundaryCondition> m_mpi_boundary_conditions;
};

#endif // MESH_HPP

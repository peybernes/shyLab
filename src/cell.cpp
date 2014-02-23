#include "cell.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>

#include "mesh_elements.hpp"

const int num_buckets = 1024;

Cell::Cell() {
}

int Cell::Dimension() {

  int result = -1;

  switch (m_cell_type) {

  case SEGMENT:
    result = 1;
    break;

  case TRIANGLE:
    result = 2;
    break;

  case QUAD:
    result = 2;
    break;

  case TETRA:
    result = 3;
    break;

  case HEXA:
    result = 3;
    break;

  case PRISM:
    result = 3;
    break;

  case PYRAMID:
    result = 3;
    break;

  default:
    assert(0);

  }

  return result;
}

Cell::Cell(MeshElement cell_type, const std::vector<int>& vertices):
  m_cell_type(cell_type), m_vertices(vertices) {}

Cell::Cell (MeshElement cell_type, const std::vector<int>& vertices, 
	    const std::vector<int>& partition_ids):
  m_cell_type(cell_type), m_vertices(vertices), m_partition_ids(partition_ids) {}

// Necessary to use FaceEntry as a key in a std::map. Not necessary
// to use FaceEntry as a key in a std::tr1::unordered_map or
// boost::unordered_map
bool Cell::operator< (const Cell &other) const {
  // CAUTION : We assume that the indices are sorted... (implemented
  // in any constructor).
  return std::lexicographical_compare(this->m_vertices.begin(), 
				      this->m_vertices.end(), 
				      other.m_vertices.begin(), 
				      other.m_vertices.end());
}

// Hash function to make FaceEntry work with boost::hash (and thus
// with boost::unordered_map).
std::size_t hash_value(Cell const& input) {

  // This simple hash is actually more efficient than the ones below
  // with boost::hash_combine, because there are very few
  // collisions: there are few edges containing one given vertice.
  // boost::hash<int> int_hash;

  // Semms to work better than with input.m_vertices[0]...
  //return int_hash(input.m_vertices[1]);
  return input.m_vertices[1];
  //return input.m_vertices[0];
  
  // // From Christer Ericson "Realtime collision detection".
  // const unsigned int magic1 = 0x8da6b343;
  // const unsigned int magic2 = 0xd8163841;
  
  // return (input.m_vertices[0] * magic1  + 
  // 	    input.m_vertices[1] * magic2) % num_buckets;
  
  // std::size_t seed = 0;
  // // boost::hash_combine(seed, input.m_vertices[0]);
  // // boost::hash_combine(seed, input.m_vertices[1]);
  // return seed;
}

std::vector<Cell> Cell::ComputeFaceEntries() {

  int taille_face = 0;
  int nb_faces = 0;

  switch (m_cell_type) {

  case TRIANGLE:
    taille_face = taille_faces_tri3;
    nb_faces = nb_faces_tri3;
    break;

  case QUAD:
    taille_face = taille_faces_quad4;
    nb_faces = nb_faces_quad4;
    break;

  case TETRA:
    taille_face = taille_faces_tetra4;
    nb_faces = nb_faces_tetra4;
    break;

  case HEXA:
    taille_face = taille_faces_hexa8;
    nb_faces = nb_faces_hexa8;
    break;

  case PRISM:
    taille_face = taille_faces_prism6;
    nb_faces = nb_faces_prism6;
    break;

  case PYRAMID:
    taille_face = taille_faces_pyr5;
    nb_faces = nb_faces_pyr5;
    break;

  default:
    assert(0);

  }
  
  std::vector<Cell> face_entries;

  for (int id_face = 0; id_face < nb_faces; ++id_face) {

    std::vector<int> face_vertices;

    MeshElement face_type = OTHER;

    switch (m_cell_type) {

    case TRIANGLE:
      face_type = SEGMENT;
      break;

    case QUAD:
      face_type = SEGMENT;
      break;

    case TETRA:
      face_type = TRIANGLE;
      break;

    case HEXA:
      face_type = QUAD;
      break;

    case PRISM:
      if ((id_face == 0) || (id_face == 2))
	face_type = QUAD;
      else
	face_type = TRIANGLE;
      break;

    case PYRAMID:
      if (id_face == 0)
	face_type = QUAD;
      else
	face_type = TRIANGLE;
      break;

    default:
      assert(0);

    }

    for (int i = 0; i < taille_face; ++i) {

      int id_vertice = DUMMY;

      switch (m_cell_type) {

      case TRIANGLE:
	id_vertice = indices_tri3[id_face][i];
	break;

      case QUAD:
	id_vertice = indices_quad4[id_face][i];
	break;

      case TETRA:
	id_vertice = indices_tetra4[id_face][i];
	break;

      case HEXA:
	id_vertice = indices_hexa8[id_face][i];
	break;

      case PRISM:
	id_vertice = indices_prism6[id_face][i];
	break;

      case PYRAMID:
	id_vertice = indices_pyr5[id_face][i];
	break;

      default:
	assert(0);

      }

      if (id_vertice != DUMMY)
	face_vertices.push_back(m_vertices.at(id_vertice));

    }
  
    face_entries.push_back(Cell(face_type, face_vertices));
       
  }

  return face_entries;
}

bool Cell::operator== (const Cell& other) const {

  return (m_vertices == other.m_vertices);

}

std::ostream& operator << (std::ostream& os, Cell cell) {
  
  os << "Vertices (C and VTK numbering) :";

  for (int i = 0; i < 3; ++i)
    os << " " << cell.vertices().at(i);

  os << "\n";
  
  std::vector<Cell> face_entries = cell.ComputeFaceEntries();

  os << "faces: ";
  
  for (int i = 0; i < 3; ++i)
    os << face_entries.at(i);

  os << "\n";
  
  return os;
}

std::vector<int> Cell::vertices() const {
  return m_vertices;
}

std::vector<int> Cell::partition_ids() const {
  return m_partition_ids;
}


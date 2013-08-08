#include "mesh.hpp"

#include <stdint.h>

#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>

#include <boost/unordered_map.hpp>

#include "array_io.hpp"
#include "byte_swap.hpp"
#include "cell_variable_metadata.hpp"
#include "face_variable_metadata.hpp"
#include "geom.hpp"
#include "mesh_elements.hpp"
#include "serialize.hpp"
#include "token_parser.hpp"
#include "variable_store.hpp"
#include "vertice_variable_metadata.hpp"

typedef int32_t vtkint_t;
typedef int32_t vtktype_t;

#define STRINGIFY(s) #s

Mesh::Mesh():
  m_filename("") {}

Mesh::Mesh(int nb_processes, int process_rank, const std::string& filename):
  m_filename(filename) {
  
  if (!m_filename.empty()) {
    // Read mesh from file.
    std::ifstream ifs(m_filename.c_str());
    ReadGmsh(nb_processes, process_rank, &ifs);
    ifs.close();
  }
}

///
void Mesh::ReadGmsh(int nb_processes, int process_rank, std::istream* input_stream_ptr) {

  int filetype = 0;             // 0 --> ascii, 1 --> binary
  int datasize = sizeof(double); // the gmsh standard for now
  
  if (!*input_stream_ptr) {
    std::cerr << "Error while reading mesh from stream. Aborting.\n";
    std::abort();
  }

  m_mpi_boundary_conditions.assign(nb_processes, MpiBoundaryCondition());

  std::string line = "";                // string holding the current line

  while (getline(*input_stream_ptr ,line)) { // Iterate through stream, line by line

    if (line.substr(0, 2) == std::string("//")) // ignore comments
      {}

    // Read file header
    else if (line.substr(0, 11) == std::string("$MeshFormat")) {

      std::cerr << "Parsing header..." << std::endl;
      getline(*input_stream_ptr, line);

      std::istringstream isstream(line);
      std::string token;
      // get file version -- for now, do nothing
      getline(isstream, token, ' ');

      // get file type (ascii or binary)
      getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &filetype)){
        std::cerr << "  bad input for gmsh filetype (expecting an integer):"
		  << token << std::endl;
	std::abort();
      }

      if (!(filetype <= 1)) {
        std::cerr << "  invalid value for gmsh filetype (expecting 0 or 1): "
		  << filetype << std::endl;
	std::abort();
      }

      if (filetype == 0)
        std::cerr << "  gmsh filetype: ascii" << std::endl;

      else
        std::cerr << "  gmsh filetype: binary" << std::endl;

      std::cerr << "done.\n" << std::endl;
    }


    // Read mesh nodes
    else if (line.substr(0, 6) == std::string("$Nodes")) {
    
      std::cerr << "Parsing nodes..." << std::endl;
      getline(*input_stream_ptr, line);

      // Get node numbers
      int nb_nodes;

      if (ParseToken<int>(line, &nb_nodes))
	std::cerr << "  number of nodes: " << nb_nodes << std::endl;

      else {
	std::cerr << "  bad input (number of nodes): "
		  << "  expected an integer, got \'" 
		  << line << "\'" << std::endl;
	std::abort();
      }

      m_nb_vertices = nb_nodes;

      m_vertices_x.assign(m_nb_vertices, 0.0);
      m_vertices_y.assign(m_nb_vertices, 0.0);
      m_vertices_z.assign(m_nb_vertices, 0.0);

      m_vertices_map.assign(m_nb_vertices, -1);

      // Get node list

      // ascii
      if (filetype == 0) {
	
	for (int i = 0; i < m_nb_vertices; ++i) {

	  getline(*input_stream_ptr, line);
	  
	  std::istringstream isstream(line);
	  std::string token;

	  // Get node index.
	  getline(isstream, token, ' ');

	  int node_index = -1;

	  if (!ParseToken<int>(token, &node_index)) {
	    
	    std::cerr << " bad input (node index): "
		      << " expected an integer, got \'"
		      << token << "\'" << std::endl;
	    std::abort();

	  }
	  
	  // Get node x coordinate.
	  getline(isstream, token, ' ');

	  RealType node_x = 0.0;

	  if (!ParseToken<RealType>(token, &node_x)) {
	    
	    std::cerr << " bad input (node x coordinate): "
		      << " expected a floating point value, got \'"
		      << token << "\'" << std::endl;
	    std::abort();
	    
	  }

	  // Get node y coordinate.
	  getline(isstream, token, ' ');

	  RealType node_y = 0.0;

	  if (!ParseToken<RealType>(token, &node_y)) {
	    
	    std::cerr << " bad input (node y coordinate): "
		      << " expected a floating point value, got \'"
		      << token << "\'" << std::endl;
	    std::abort();
	    
	  }

	  // Get node z coordinate.
	  getline(isstream, token, ' ');
	  
	  RealType node_z = 0.0;
	  
	  if (!ParseToken<RealType>(token, &node_z)) {
	    
	    std::cerr << " bad input (node z coordinate): "
		      << " expected a floating point value, got \'"
		      << token << "\'" << std::endl;
	    std::abort();
	    
	  }
	  
	  m_vertices_map.at(i) = node_index;
	  m_vertices_x.at(i) = node_x;
	  m_vertices_y.at(i) = node_y;
	  m_vertices_z.at(i) = node_z;

	}
	
      }      

      // binary
      else { 

	for (int i = 0; i < m_nb_vertices; ++i) {

	  int node_index = -1;
	  input_stream_ptr->read((char*)&node_index, 4);

	  double x = 0.0;
	  input_stream_ptr->read((char*)&x, 8);

	  double y = 0.0;
	  input_stream_ptr->read((char*)&y, 8);

	  double z = 0.0;
	  input_stream_ptr->read((char*)&z, 8);

	  m_vertices_map.at(i) = node_index;
	  m_vertices_x.at(i) = x;
	  m_vertices_y.at(i) = y;
	  m_vertices_z.at(i) = z;

	}

        // Sanity check: is the next character a newline ?
        assert(input_stream_ptr->peek() == '\n');
      }

      std::cerr << "done.\n" << std::endl;

    }

    // Read mesh elements
    else if (line.substr(0, 9) == std::string("$Elements")) {
      
      std::cerr << "Parsing elements..." << std::endl;

      // Get number of elements of the global mesh.
      int nb_elements_global;
      getline(*input_stream_ptr, line);

      if (ParseToken<int>(line, &nb_elements_global)) {
       
	if (process_rank == 0)
	  std::cerr << "  global number of elements (including boundaries): "
		    << nb_elements_global << std::endl;

      }

      else {

	if (process_rank == 0)
          std::cerr << "  bad input (number of elements): "
		    << "  expected an integer, got \'" 
		    << line << "\'" << std::endl;
	
	std::abort();

      }
      
      if (filetype == 0) { // ascii
	
	for (int i = 0; i < nb_elements_global; ++i) {
	  
	  getline(*input_stream_ptr, line);
	  
	  std::istringstream isstream(line);
	  std::string token;

	  // Get element index.
	  getline(isstream, token, ' ');

	  int cell_index = -1;

	  if (!ParseToken<int>(token, &cell_index)) {
	    
	    std::cerr << " bad input (cell index): "
		      << " expected an integer, got \'"
		      << token << "\'" << std::endl;
	    std::abort();

	  }

	  // Get element type.
	  getline(isstream, token, ' ');

	  int element_type = -1;

	  if (!ParseToken<int>(token, &element_type)) {
	    
	    std::cerr << " bad input (element type): "
		      << " expected an integer, got \'"
		      << token << "\'" << std::endl;
	    std::abort();

	  }

	  MeshElement current_cell_type = OTHER;

	  switch (element_type) {

	  case GMSH_POINT:
	    break;
	    
	  case GMSH_SEGMENT:
	    break;
	    
	  case GMSH_TRIANGLE:
	    current_cell_type = TRIANGLE;
	    break;

	  case GMSH_QUAD:
	    current_cell_type = QUAD;
	    break;

	  case GMSH_TETRA:
	    current_cell_type = TETRA;
	    break;

	  case GMSH_HEXA:
	    current_cell_type = HEXA;
	    break;

	  case GMSH_PRISM:
	    current_cell_type = PRISM;
	    break;

	  case GMSH_PYRAMID:
	    current_cell_type = PYRAMID;
	    break;

	  default:
	    assert(0);

	  }

	  // Get number of tags.
	  getline(isstream, token, ' ');

	  int nb_tags = -1;

	  if (!ParseToken<int>(token, &nb_tags)) {
	    
	    std::cerr << " bad input (number of tags): "
		      << " expected an integer, got \'"
		      << token << "\'" << std::endl;
	    std::abort();

	  }

	  // Get tags (we use gmsh default convention : first tag is
	  // physical id, second is elementary id, third, is number of
	  // partitions, then the list of partition ids the cell
	  // belongs to (if the cell is a ghost cell, the partition id
	  // is negative).

	  int physical_id = -1;
	  int elementary_id = -1;

	  int nb_partitions = 1;
	  std::vector<int> partition_ids;

	  int nb_vertices = 0;
	  std::vector<int> vertice_ids;
	  
	  if (nb_tags >= 1) {

	    // Get physical id.
	    getline(isstream, token, ' ');

	    if (!ParseToken<int>(token, &physical_id)) {
	    
	      std::cerr << " bad input (physical id): "
			<< " expected an integer, got \'"
			<< token << "\'" << std::endl;
	      std::abort();

	    }

	  }

	  if (nb_tags >= 2) {

	    // Get elementary id.
	    getline(isstream, token, ' ');

	    if (!ParseToken<int>(token, &elementary_id)) {
	    
	      std::cerr << " bad input (elementary id): "
			<< " expected an integer, got \'"
			<< token << "\'" << std::endl;
	      std::abort();

	    }

	  }

	  if (nb_tags >= 3) {

	    // Get number of partitions.
	    getline(isstream, token, ' ');

	    if (!ParseToken<int>(token, &nb_partitions)) {
	    
	      std::cerr << " bad input (elementary id): "
			<< " expected an integer, got \'"
			<< token << "\'" << std::endl;
	      std::abort();

	    }

	    assert(nb_tags == 3 + nb_partitions);

	    for (int j = 0; j < nb_partitions; ++j) {

	      int partition_id = 0;
	      // Get partition id.
	      getline(isstream, token, ' ');

	      if (!ParseToken<int>(token, &partition_id)) {
	    
		std::cerr << " bad input (partition id): "
			  << " expected an integer, got \'"
			  << token << "\'" << std::endl;
		std::abort();

	      }

	      partition_ids.push_back(partition_id);

	    }

	  }

	  // Get vertice ids.

	  // Get number of vertices depending on element type.
	  switch (element_type) {

	  case GMSH_SEGMENT:
	    nb_vertices = 2;
	    break;
	    
	  case GMSH_TRIANGLE:
	    nb_vertices = 3;
	    break;

	  case GMSH_QUAD:
	    nb_vertices = 4;
	    break;

	  case GMSH_TETRA:
	    nb_vertices = 4;
	    break;

	  case GMSH_HEXA:
	    nb_vertices = 8;
	    break;
	    
	  case GMSH_PRISM:
	    nb_vertices = 6;
	    break;

	  case GMSH_PYRAMID:
	    nb_vertices = 5;
	    break;

	  case GMSH_POINT:
	    nb_vertices = 1;
	    break;

	  default:
	    assert(0);

	  }

	  for (int j = 0; j < nb_vertices; ++j) {

	    int vertice_id = -1;
	    // Get vertice id.
	    getline(isstream, token, ' ');

	    if (!ParseToken<int>(token, &vertice_id)) {
	    
	      std::cerr << " bad input (vertice id): "
			<< " expected an integer, got \'"
			<< token << "\'" << std::endl;
	      std::abort();

	    }

	    vertice_ids.push_back(vertice_id);

	  }

	  // We only fill in maximum dimension elements.

	  bool element_is_cell = false;

	  element_is_cell =
	    ((ambient_dimension == 2) &&
	     ((element_type == GMSH_TRIANGLE) || 
	      (element_type == GMSH_QUAD))) 
	    ||
	    ((ambient_dimension == 3) &&
	     ((element_type == GMSH_TETRA) || 
	      (element_type == GMSH_HEXA)  ||
	      (element_type == GMSH_PRISM) ||
	      (element_type == GMSH_PYRAMID)));
	    
	  if (element_is_cell) {

	    // Scan for ghost cells.
	    int ghost_cell_counter = 0;

	    for (int j = 0; j < nb_partitions; ++j) {

	      if (partition_ids.at(j) < 0)
		++ghost_cell_counter;
	      
	    }

	    const bool parallel_mode = false;

	    bool belongs_to_current_mpi_process = false;

	    if (parallel_mode == true) {

	    // If a cell belongs to several partition, it must be 
	    // a ghost cell in every partition except 1 (owner
	    // computes rule).
	    assert((nb_partitions == 1) || 
		   ((nb_partitions > 1) && 
		    (ghost_cell_counter == nb_partitions - 1)));

	    // Update MPI ghost node info for elements belonging to
	    // several partitions.

	    for (int j = 0; j < nb_partitions; ++j) {
	      
	      const int partition_id = partition_ids.at(j);
	      
	      // Gmsh numbers partition Fortran style (e.g. from 1 to n).
	      // Moreover, negative number means a ghost node.
	      if (process_rank == std::abs(partition_id) - 1) {

		belongs_to_current_mpi_process = true;
		
		// Case of a cell belonging to several partitions.
		if (nb_partitions > 1) {

		  // The cell is a true cell on this mpi process.
		  if (partition_id > 0) {
		    
		    // Current cell values must be send to all other
		    // partition sharing it.
		    for (int k = 0; k < nb_partitions; ++k) {
		      
		      if (k != j) {

			const int other_partition_id = 
			  partition_ids.at(k);
			assert(other_partition_id < 0);
		
			m_mpi_boundary_conditions.
			  at(std::abs(other_partition_id) - 1).
			  send_cell_ids.push_back(cell_index);
	
		      }

		    }
		    
		  }
		  
		  // The cell is a ghost cell on this mpi process.
		  else if (partition_id < 0) {

		    // We determine the partition where the cell is a
		    // true cell, and mark it as the partition from
		    // which we must recieve data.
		    for (int k = 0; k < nb_partitions; ++k) {

		      const int other_partition_id = partition_ids.at(k);

		      if (other_partition_id > 0) {
			m_mpi_boundary_conditions.
			  at(std::abs(other_partition_id) - 1).
			  recieve_cell_ids.push_back(cell_index);
		      }

		    }

		  }

		  else
		    assert(0);

		}
		
	      }

	    }
	    }
	    
	    if (belongs_to_current_mpi_process == true || parallel_mode == false) {

	      m_cells.push_back(Cell(current_cell_type, vertice_ids, partition_ids));
	      m_unique_to_local_cell_ids[cell_index] = m_cells.size() - 1;

	    }
	    
	    }

	  }
	

      }
 
      else { // binary

	// Currently not up to date.
	assert(0);
	
	int nb_elements_read = 0;

	// In gmsh binary file format, elements are grouped by their types.
	while (nb_elements_read < nb_elements_global) {

	  // First we read header.

	  // Type id.
	  int elements_type = -1;
	  input_stream_ptr->read((char*)&elements_type, sizeof(int));

	  // Number of element of a given type.
	  int nb_elements_in_group = 0;
	  input_stream_ptr->read((char*)&nb_elements_in_group, sizeof(int));

	  // Number of tags per element.
	  int nb_tags = 0;
	  input_stream_ptr->read((char*)&nb_tags, sizeof(int));

	  nb_elements_read += nb_elements_in_group;

	  int number_of_vertices = 0;
	    
	  switch (elements_type) {
	    
	  case 1:
	    number_of_vertices = 2;
	    break;

	  case 2:
	    number_of_vertices = 3;
	    break;

	  case 15:
	    number_of_vertices = 1;
	    break;

	  default:
	    assert(0);
	    break;

	  }

	  for (int i = 0; i < nb_elements_in_group; ++i) {

	    int node_index = -1;
	    input_stream_ptr->read((char*)&node_index, sizeof(int));	    

	    int tags[nb_tags];
	    input_stream_ptr->read((char*)&tags, sizeof(tags));	    
	    
	    // By default tags[0] is physical region id, tags[1] is
	    // elementary region id, tags[2] is number of partitions,
	    // tags[3] is the first partition number, tags[4] is the
	    // opposite of other partition number (case of a ghost
	    // cell).

	    const int nb_partitions = 1;

	    std::vector<int> partition_ids;
	    
	    // Partition info is provided by the mesh file.
	    if (nb_tags == 3) {

	      partition_ids.push_back(tags[2]);

	    }

	    // Partition info is provided by the mesh file.
	    if (nb_tags > 3) {

	      const int nb_partitions = tags[2];

	      for (int id_partition = 0; id_partition < nb_partitions;
		   ++id_partition) {

		// partition ids are counted from 1 in gmsh (Fortran style).
		partition_ids.push_back(std::abs(tags[2 + id_partition])- 1);

	      }
	    }

	    // Partition info is not provided (thus the mesh is
	    // sequential).
	    else {

	      partition_ids.push_back(1);
	      
	    }

	    int vertices[number_of_vertices];
	    input_stream_ptr->read((char*)&vertices, sizeof(vertices));	    
	    
	    if (elements_type == 2) {

	      std::vector<int> vertice_ids;
	      vertice_ids.push_back(vertices[0]);
	      vertice_ids.push_back(vertices[1]);
	      vertice_ids.push_back(vertices[2]);

	      bool belongs_to_current_mpi_process = false;

	      for (int i = 0; i < nb_partitions; ++i) {
		
	      	const int partition_id = partition_ids.at(i);

	      	// Gmsh numbers partition Fortran style (e.g. from 1 to n).
	      	// Moreover, negative number means a ghost node.
	      	if (process_rank + 1 == std::abs(partition_id))
	      	  belongs_to_current_mpi_process = true;

	      }

	      if (belongs_to_current_mpi_process == true) {

		//m_cells.push_back(Cell(vertice_ids, partition_ids));

	      }

	    }

	  }
	  
	}

	// Sanity check : is the next character a new line ?
        assert(input_stream_ptr->peek() == '\n');

      }

      m_nb_cells = m_cells.size();

      std::cerr << "  number of cells : " << m_nb_cells << "\n";

      std::cerr << "done.\n" << std::endl;
    }

    else // ignore other lines
      {}

    }
    

}

/// Reverse Cuthill McKee renumbering algorithm.
///
/// This algorithm reorders the mesh cells so the mesh bandwith is
/// reduced. It is a variant of a breadth first search traversal. It
/// returns a mesh for which neighbouring cells tend to be close in
/// indexing.
void Mesh::ReverseCuthillMcKeeRenumbering() {

  m_cells_map = std::vector<int>(m_nb_cells, -1);
  
  // Heuristic choice for the seed cell. For now, we choose the first
  // boundary cell we meet.

  int seed_cell_id = -1;

  for (int i = 0; i < m_nb_faces; ++i) {
    
    if (m_face_back_cell_id.at(i) == -1) {
      seed_cell_id = m_face_front_cell_id.at(i);
      break;

    }
    
    if (i == m_nb_faces - 1)
      assert(0);
  }

  assert((0 < seed_cell_id) && (seed_cell_id < m_nb_cells));

  // Traverse the mesh graph in a modified breadth first search
  // fashion, starting from cell seed_cell_id. m_cells_map holds the
  // renumbering permutation.

  int counter = 0;

  std::queue<int> next_candidate_cells;

  next_candidate_cells.push(seed_cell_id);
  m_cells_map.at(seed_cell_id) = seed_cell_id;
  
  ++counter;

  while (!(next_candidate_cells.empty())) {

    // first cell id in the queue.
    const int current_cell_id = next_candidate_cells.front();
    next_candidate_cells.pop();

    // Get the current cell neighbor ids.
    for (int i = 0; i < 3; ++i) {
      
      const int neighbor_cell_id = 
	m_cell_cell_connectivity.at(i).at(current_cell_id);
      
      if ((neighbor_cell_id != -1) && 
	  (m_cells_map.at(neighbor_cell_id) == -1)) {
	// This cell id is unvisited.
	m_cells_map.at(neighbor_cell_id) = counter;
	next_candidate_cells.push(neighbor_cell_id);
	++counter;
      }
    }
  }
}

void Mesh::WriteMeshBandwith(std::ofstream &os) {

  // write mesh connectivity in PPM format.
  if (m_nb_cells > 5000)
    std::cerr << "Mesh size too big, mesh bandwith not written.\n";

  else {
    os << "P2\n";
    os << "# File generated by Shy: mesh bandwith\n";
    os << m_nb_cells << "\n";
    os << m_nb_cells << "\n";
    os << 255 << "\n";

    for (int j = 0; j < m_nb_cells; ++j) {
      for (int i = 0; i < m_nb_cells; ++i) {

	std::vector<int> neighbors(3);

	for (int id_neighbor = 0; id_neighbor < 3; ++id_neighbor)
	  neighbors[id_neighbor] = 
	    m_cell_cell_connectivity.at(id_neighbor).at(i);

        if (neighbors[0] == j || 
	    neighbors[1] == j ||
	    neighbors[2] == j)
          os << 255 << " ";
        else
          os << 0 << " ";
      }
      os << "\n";
    }
  }

}

/// Compute face<->cell and cell->cell connectivity. Uses TR1
/// unordered_map (= hash table) container.
void Mesh::ComputeFaceCellConnectivity() {

  m_face_front_cell_id.clear();
  m_face_back_cell_id.clear();

  m_cell_cell_connectivity.clear();
  
  typedef 
    boost::unordered_map<Cell, std::vector<int> > connectivity_hashtable_t;
  
  // Will be discarded at the end of function. Connectivity indices
  // will be put in m_face_cell_connectivity_{front,back}_
  connectivity_hashtable_t face_cell_connectivity;

  std::cerr << "Computing mesh connectivity..." << std::endl;
  
  // Walk throught mesh cells, and store their faces in a hash table.
  for (int j = 0; j < m_nb_cells; ++j) {

    Cell cell = m_cells.at(j);

    std::vector<Cell> cell_faces = cell.ComputeFaceEntries();

    for (int i = 0; i < 3; ++i)
      // A O(log N) cost is hidden here (hash table insertion)
      face_cell_connectivity[cell_faces.at(i)].push_back(j);
    
  }

  m_nb_faces = face_cell_connectivity.size();

  // Allocation of face/cell connectivity arrays.
  m_face_front_cell_id.reserve(m_nb_faces);
  m_face_back_cell_id.reserve(m_nb_faces);

  // Allocation of cell/cell connectivity arrays.
  for (int i = 0; i < 3; ++i)
    m_cell_cell_connectivity.push_back(std::vector<int>(m_nb_cells, -1));

  connectivity_hashtable_t::iterator face_iterator;
  for (face_iterator = face_cell_connectivity.begin(); 
       face_iterator != face_cell_connectivity.end(); ++face_iterator) {

    const Cell key = face_iterator->first;

    // The face belongs to one cell : it is a boundary face.
    if (face_cell_connectivity[key].size() == 1)
      face_cell_connectivity[key].push_back(-1);

    const int cell_id_1 = face_cell_connectivity[key].at(0);
    const int cell_id_2 = face_cell_connectivity[key].at(1);

    const int front_cell_id = std::max(cell_id_1, cell_id_2);
    const int back_cell_id = std::min(cell_id_1, cell_id_2);

    // Face cells.
    m_face_front_cell_id.push_back(front_cell_id);
    m_face_back_cell_id.push_back(back_cell_id);
    
    // Cell cell.
    if ((cell_id_1 != -1) && (cell_id_2 != -1)) {

      for (int i = 0; i < 3; ++i) {

      	if (m_cell_cell_connectivity.at(i).at(cell_id_1) == -1) {
      	  m_cell_cell_connectivity.at(i).at(cell_id_1) = cell_id_2;
      	  break;
      	}

      	else {
      	  if (i == 2)
      	    assert(0);
      	}

      }

      for (int i = 0; i < 3; ++i) {

      	if (m_cell_cell_connectivity.at(i).at(cell_id_2) == -1) {
      	  m_cell_cell_connectivity.at(i).at(cell_id_2) = cell_id_1;
      	  break;
      	}

      	else {
      	  if (i == 2)
      	    assert(0);
      	}

      }
    }

    // Face vertices.
    m_face_vertices_id.push_back(key.vertices());

    if (face_cell_connectivity[key].size() == 2) {
      // The face belongs to cells cell_id_1 and cell_id_2, so cell_id_1 and
      // cell_id_2 are neighbors.
    }

    // Should never happen.
    else
      assert(0);
  }

  std::cerr << "  * number of faces : " << face_cell_connectivity.size()
	    << "\n";
  std::cerr << "  * load factor : " 
	    << face_cell_connectivity.load_factor() << "\n";
  std::cerr << "...done" << std::endl;
}

/// Compute cell-->face connectivity from face-->cell connectivity.
void Mesh::ComputeCellFaceConnectivity() {

  /// Reset.
  m_cell_signed_face1_id.assign(m_nb_cells, -1);
  m_cell_signed_face2_id.assign(m_nb_cells, -1);
  m_cell_signed_face3_id.assign(m_nb_cells, -1);
  m_cell_signed_face4_id.assign(m_nb_cells, -1);

  // Ugly hack, because 0 == -0.

  for (int i = 0; i < m_nb_faces; ++i) {

    const int cell_id1 = m_face_back_cell_id[i];
    const int cell_id2 = m_face_front_cell_id[i];

    if (cell_id1 != -1) {
      
      if (m_cell_signed_face1_id[cell_id1] == -1)
	m_cell_signed_face1_id[cell_id1] = i + 1;

      else if (m_cell_signed_face2_id[cell_id1] == -1)
	m_cell_signed_face2_id[cell_id1] = i + 1;

      else if (m_cell_signed_face3_id[cell_id1] == -1)
	m_cell_signed_face3_id[cell_id1] = i + 1;

      else if (m_cell_signed_face4_id[cell_id1] == -1)
	m_cell_signed_face4_id[cell_id1] = i + 1;

      else
	assert(0);

    }

    if (cell_id2 != -1) {

      if (m_cell_signed_face1_id[cell_id2] == -1)
	m_cell_signed_face1_id[cell_id2] = -i - 1;
    
      else if (m_cell_signed_face2_id[cell_id2] == -1)
	m_cell_signed_face2_id[cell_id2] = -i - 1;
      
      else if (m_cell_signed_face3_id[cell_id2] == -1)
	m_cell_signed_face3_id[cell_id2] = -i - 1;

      else if (m_cell_signed_face4_id[cell_id1] == -1)
	m_cell_signed_face4_id[cell_id1] = - i - 1;
      
      else
	assert(0);

    }

  }
}

void Mesh::WriteFacesVTKXmlAscii(std::ofstream* os_ptr) {
  
  // Write vertices coordinates.

  *os_ptr << "<Points>\n"
	  << "<DataArray "
	  << "type=\"Float32\" "
	  << "NumberOfComponents=\"3\" "
	  << "format=\"ascii\">\n";

  for (int i = 0; i < m_nb_vertices; ++i) {
    
    *os_ptr << m_vertices_x.at(i) << " " 
	    << m_vertices_y.at(i) << " " 
	    << m_vertices_z.at(i) << "\n";

  }
  
  *os_ptr << "</DataArray>\n"
	  << "</Points>\n";
  
  // Write cells connectivities (1., 2., 3.).
  int current_offset = 0;
  
  *os_ptr << "<Cells>\n";
  
  // 1. Nodal connectivities.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"connectivity\" "
	  << "format=\"ascii\">\n";

  for (int i = 0; i < m_nb_cells; ++i) {

    const int nb_cell_vertices = m_cells.at(i).vertices().size();

    for (int j = 0; j < nb_cell_vertices; ++j)
      *os_ptr << m_cells.at(i).vertices().at(j) - 1 << " ";

    *os_ptr << "\n";
    
  }

  *os_ptr << "</DataArray>\n";

  // 2. Offsets.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"offsets\" "
	  << "format=\"ascii\">\n";
  
  for (int i = 0; i < m_nb_cells; ++i) {

    const int nb_cell_vertices = m_cells.at(i).vertices().size();

    current_offset += nb_cell_vertices;

    *os_ptr <<  current_offset << "\n";

  }

  *os_ptr << "</DataArray>\n";

  // 3. Cell types.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"types\" "
	  << "format=\"ascii\">\n";
  
  for (int i = 0; i < m_nb_cells; ++i) {

    MeshElement cell_type = m_cells.at(i).cell_type();

    int vtk_cell_type = -1;

    switch (cell_type) {

    case TRIANGLE:
      vtk_cell_type = VTK_TRIANGLE;
      break;

    case QUAD:
      vtk_cell_type = VTK_QUAD;
      break;

    case TETRA:
      vtk_cell_type = VTK_TETRA;
      break;

    case HEXA:
      vtk_cell_type = VTK_HEXAHEDRON;
      break;

    case PRISM:
      vtk_cell_type = VTK_WEDGE;
      break;

    case PYRAMID:
      vtk_cell_type = VTK_PYRAMID;
      break;

    default:
      assert(0);

    }

    *os_ptr << vtk_cell_type << "\n";

  }

  *os_ptr << "</DataArray>\n";

  *os_ptr << "</Cells>\n";

}

void Mesh::WriteHeaderVTKXmlAscii(std::ofstream* os_ptr) {

    *os_ptr << "<?xml version=\"1.0\"?>\n"
	    << "<VTKFile type=\"UnstructuredGrid\" "
	    << "version=\"0.1\" "
	    << "byte_order=\"LittleEndian\">\n";

    *os_ptr << "<UnstructuredGrid>\n"
	    << "<Piece NumberOfPoints=\"" << m_nb_vertices << "\" "
	    << "NumberOfCells=\"" << m_nb_cells << "\">\n";

}

void Mesh::WriteFooterVTKXmlAscii(std::ofstream* os_ptr) {

    *os_ptr << "</Piece>\n"
	    << "</UnstructuredGrid>\n";

    *os_ptr << "</VTKFile>\n";

}

void Mesh::WriteVTKXmlAscii(std::ofstream* os_ptr) {
  
  // Write vertices coordinates.

  *os_ptr << "<Points>\n"
	  << "<DataArray "
	  << "type=\"Float32\" "
	  << "NumberOfComponents=\"3\" "
	  << "format=\"ascii\">\n";

  for (int i = 0; i < m_nb_vertices; ++i) {
    
    *os_ptr << m_vertices_x.at(i) << " " 
	    << m_vertices_y.at(i) << " " 
	    << m_vertices_z.at(i) << "\n";

  }
  
  *os_ptr << "</DataArray>\n"
	  << "</Points>\n";
  
  // Write cells connectivities (1., 2., 3.).
  int current_offset = 0;
  
  *os_ptr << "<Cells>\n";
  
  // 1. Nodal connectivities.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"connectivity\" "
	  << "format=\"ascii\">\n";

  for (int i = 0; i < m_nb_cells; ++i) {

    const int nb_cell_vertices = m_cells.at(i).vertices().size();

    for (int j = 0; j < nb_cell_vertices; ++j)
      *os_ptr << m_cells.at(i).vertices().at(j) - 1 << " ";

    *os_ptr << "\n";
    
  }

  *os_ptr << "</DataArray>\n";

  // 2. Offsets.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"offsets\" "
	  << "format=\"ascii\">\n";
  
  for (int i = 0; i < m_nb_cells; ++i) {

    const int nb_cell_vertices = m_cells.at(i).vertices().size();

    current_offset += nb_cell_vertices;

    *os_ptr <<  current_offset << "\n";

  }

  *os_ptr << "</DataArray>\n";

  // 3. Cell types.
  *os_ptr << "<DataArray "
	  << "type=\"Int32\" "
	  << "Name=\"types\" "
	  << "format=\"ascii\">\n";
  
  for (int i = 0; i < m_nb_cells; ++i) {

    MeshElement cell_type = m_cells.at(i).cell_type();

    int vtk_cell_type = -1;

    switch (cell_type) {

    case TRIANGLE:
      vtk_cell_type = VTK_TRIANGLE;
      break;

    case QUAD:
      vtk_cell_type = VTK_QUAD;
      break;

    case TETRA:
      vtk_cell_type = VTK_TETRA;
      break;

    case HEXA:
      vtk_cell_type = VTK_HEXAHEDRON;
      break;

    case PRISM:
      vtk_cell_type = VTK_WEDGE;
      break;

    case PYRAMID:
      vtk_cell_type = VTK_PYRAMID;
      break;

    default:
      assert(0);

    }

    *os_ptr << vtk_cell_type << "\n";

  }

  *os_ptr << "</DataArray>\n";

  *os_ptr << "</Cells>\n";

}

void Mesh::WriteVTKBinary(std::ostream &os, bool swap) {

  if (!os) {
    std::cerr << "Error while opening file for writing. Aborting.\n";
    std::abort();
  }

  // Header.
  os << "# vtk DataFile Version 3.0\nOutput from Shy.\n";
  os << "BINARY\nDATASET UNSTRUCTURED_GRID\n\n";

  // Mesh vertices coords
  os << "POINTS " << m_nb_vertices << " " << STRINGIFY(RealType) << "\n";
  
  for (int i = 0; i < m_nb_vertices; ++i) {

    RealType x = m_vertices_x.at(i);

    if (swap == true)
      x = ByteSwap<RealType>(x);

    os.write((char*)&x, sizeof(RealType));

    RealType y = m_vertices_y.at(i);

    if (swap == true)
      y = ByteSwap<RealType>(y);

    os.write((char*)&y, sizeof(RealType));

    RealType z = static_cast<RealType>(0.0);

    if (swap == true)
      z = ByteSwap<RealType>(z);

    os.write((char*)&z, sizeof(RealType));

  }

  os << "\n";

  // write cells
  os << "CELLS " << m_nb_cells << " " << (3 + 1) * m_nb_cells << "\n";

  for (int i = 0; i < m_nb_cells; ++i) {

    int nb_cell_vertices = 3;
    
    int nb_vertices_dump = nb_cell_vertices;

    if (swap == true)
      nb_vertices_dump = ByteSwap<int>(nb_cell_vertices);

    os.write((char*)&nb_vertices_dump, sizeof(vtkint_t));

    for (int j = 0; j < nb_cell_vertices; ++j) {
      
      int vertex_indice = m_cells.at(i).vertices()[j] - 1;

      vertex_indice = static_cast<vtkint_t>(vertex_indice);

      if (swap == true)
	vertex_indice = ByteSwap<int>(vertex_indice);

      os.write((char*)&vertex_indice, sizeof(vtkint_t));
    }
  }
  
  os << "\n";

  // write cell types (5 for triangles)
  os << "CELL_TYPES " << m_nb_cells << "\n";

  for (int i = 0; i < m_nb_cells; ++i) {
    
    int cell_type = 5;
    cell_type = static_cast<vtktype_t>(cell_type);

    if (swap == true)
      cell_type = ByteSwap<vtktype_t>(cell_type);

    os.write((char*)&cell_type, sizeof(vtktype_t));

  }

  os << "\n\n";
}


int Mesh::TriangleIndex(RealType xp, RealType yp, RealType zp) {

  int result = -1;

  int i = 0;
  while ((result < 0) && (i < m_nb_cells)) {

    const RealType xa = m_vertices_x[m_cells.at(i).vertices()[0]];
    const RealType xb = m_vertices_x[m_cells.at(i).vertices()[1]];
    const RealType xc = m_vertices_x[m_cells.at(i).vertices()[2]];

    const RealType ya = m_vertices_y[m_cells.at(i).vertices()[0]];
    const RealType yb = m_vertices_y[m_cells.at(i).vertices()[1]];
    const RealType yc = m_vertices_y[m_cells.at(i).vertices()[2]];

    const RealType za = m_vertices_z[m_cells.at(i).vertices()[0]];
    const RealType zb = m_vertices_z[m_cells.at(i).vertices()[1]];
    const RealType zc = m_vertices_z[m_cells.at(i).vertices()[2]];
    
    if (IsInside(xp, yp, zp, xa, ya, za, xb, yb, zb, xc, yc, zc))
      result = i;

    ++i;     // go to next triangle
  }

  // This is possible for a mesh with degenerate triangles
  if (result < 0) {
    std::cerr << "Warning: could not find a triangle containing the point of coordinates "
              << "(" << xp << " " << yp << " " << zp << ")" << std::endl;
  }

  return result;
}

// void Mesh::ComputeGradientInterpolator() {

//   std::cerr << "Computing least square matrices for gradient ";
//   std::cerr << "reconstruction... ";

//   for ( int i = 0; i < m_nb_cells; ++i ) {

//     // We construct the least square interpolation matrix as
//     // M = sum( Di * (Di)^t ), where Di = Xi - X0

//     LSQMatrix lsq = LSQMatrix::Zero();

//     for ( int j = 0; j < ncell; ++j )

//       if ( Cells.at(i).neighbors()[j] != -1 ) {
//         Vector DeltaPoints;
// 	DeltaPoints.x() = 
// 	  CellCenters.x( Cells.at(i).neighbors()[j] ) -
//           CellCenters.x( i );
// 	DeltaPoints.y() = 
// 	  CellCenters.y( Cells.at(i).neighbors()[j] ) -
//           CellCenters.y( i );

//         DeltaPoints /= DeltaPoints.norm();
//         lsq += (DeltaPoints * DeltaPoints.transpose()).block<2,2>(0,0);
//       }

//     if ( !Cells.at(i).degenerate ) {
//       LSQMatrix inverse = LSQMatrix::Zero();
//       lsq.marked<Eigen::SelfAdjoint>().computeInverse( &inverse );
//       GradientInterpolator.push_back( inverse );
//     }

//     // In this degenerate case, the (boundary) cell has less neighbors than the
//     // ambient space dimension, and we don't reconstruct the gradient.
//     // TODO -- reconstruct the gradient using the triangles in the 2-ring
//     //else
//     GradientInterpolator.push_back( LSQMatrix::Zero() );

//   }

//   std::cerr << "done.\n";

// }

void Mesh::ComputeGeometricQuantities(VariableStore* cell_variables, 
				      VariableStore* face_variables,
				      VariableStore* vertice_variables) {

  assert(m_nb_cells == cell_variables->nb_elements());
  assert(m_nb_faces == face_variables->nb_elements());
  assert(m_nb_vertices == vertice_variables->nb_elements());

  RealType* face_centers_x = (*face_variables)(FACE_CENTERS_X);
  RealType* face_centers_y = (*face_variables)(FACE_CENTERS_Y);
  RealType* face_centers_z = (*face_variables)(FACE_CENTERS_Z);
  RealType* face_volumes = (*face_variables)(FACE_VOLUMES);

  RealType* normals_x = (*face_variables)(NORMALS_X);
  RealType* normals_y = (*face_variables)(NORMALS_Y);
  RealType* normals_z = (*face_variables)(NORMALS_Z);

  RealType* vertices_x = (*vertice_variables)(VERTICES_X);
  RealType* vertices_y = (*vertice_variables)(VERTICES_Y);
  RealType* vertices_z = (*vertice_variables)(VERTICES_Z);

  RealType* cell_centers_x = (*cell_variables)(CELL_CENTERS_X);
  RealType* cell_centers_y = (*cell_variables)(CELL_CENTERS_Y);
  RealType* cell_centers_z = (*cell_variables)(CELL_CENTERS_Z);
  RealType* cell_volumes = (*cell_variables)(CELL_VOLUMES);

  std::cerr << "Computing mesh geometric quantities... ";

  // Transfer vertice coordinates from mesh to variable.
#pragma omp parallel for
  for (int i = 0; i < m_nb_vertices; ++i) {

    vertices_x[i] = m_vertices_x.at(i);
    vertices_y[i] = m_vertices_y.at(i);
    vertices_z[i] = m_vertices_z.at(i);

  }
  
#pragma omp parallel for
  for (int i = 0; i < m_nb_cells; ++i) {

    RealType center_x = 0.0;
    RealType center_y = 0.0;
    RealType center_z = 0.0;

    const int nb_cell_vertices = m_cells.at(i).vertices().size();

    for (int j = 0; j < nb_cell_vertices; ++j) {
      
      const int cell_vertice_id = m_cells.at(i).vertices().at(j) - 1;

      assert(0 <= cell_vertice_id);
      assert(cell_vertice_id < m_nb_vertices);

      center_x += vertices_x[cell_vertice_id];
      center_y += vertices_y[cell_vertice_id];
      center_z += vertices_z[cell_vertice_id];

    }

    center_x /= nb_cell_vertices;
    center_y /= nb_cell_vertices;
    center_z /= nb_cell_vertices;

    cell_centers_x[i] = center_x;
    cell_centers_y[i] = center_y;
    cell_centers_z[i] = center_z;

    // Cell volumes

    const int cell_vertice_id1 = m_cells.at(i).vertices()[0] - 1;

    assert((0 <= cell_vertice_id1) && 
	   (cell_vertice_id1 < m_nb_vertices));

    const RealType x1 = vertices_x[cell_vertice_id1];
    const RealType y1 = vertices_y[cell_vertice_id1];
    const RealType z1 = vertices_z[cell_vertice_id1];

    const int cell_vertice_id2 = m_cells.at(i).vertices()[1] - 1;

    assert((0 <= cell_vertice_id2) && 
	   (cell_vertice_id2 < m_nb_vertices));

    const RealType x2 = vertices_x[cell_vertice_id2];
    const RealType y2 = vertices_y[cell_vertice_id2];
    const RealType z2 = vertices_z[cell_vertice_id2];

    const int cell_vertice_id3 = m_cells.at(i).vertices()[2] - 1;

    assert((0 <= cell_vertice_id3) && 
	   (cell_vertice_id3 < m_nb_vertices));

    const RealType x3 = vertices_x[cell_vertice_id3];
    const RealType y3 = vertices_y[cell_vertice_id3];
    const RealType z3 = vertices_z[cell_vertice_id3];
    
    RealType volume = -1.0;

    if (m_cells.at(i).cell_type() == TRIANGLE)
      volume = 0.5 * std::abs(Orient2d(x1, y1, z1, 
				       x2, y2, z2,
				       x3, y3, z3));
    
    else if (m_cells.at(i).cell_type() == TETRA) {

      const int cell_vertice_id4 = m_cells.at(i).vertices()[3] - 1;

      assert((0 <= cell_vertice_id4) && 
	     (cell_vertice_id4 < m_nb_vertices));

      const RealType x4 = vertices_x[cell_vertice_id4];
      const RealType y4 = vertices_y[cell_vertice_id4];
      const RealType z4 = vertices_z[cell_vertice_id4];
      
      volume = TetraedronVolume(x1, y1, z1, 
				x2, y2, z2, 
				x3, y3, z3,
				x4, y4, z4);

    }

    else 
      assert(0);

    cell_volumes[i] = volume;

  }

#pragma omp parallel for
  for (int i = 0; i < m_nb_faces; ++i) {
    
    RealType face_center_x = 0.0;
    RealType face_center_y = 0.0;
    RealType face_center_z = 0.0;

    const int nb_face_vertices = m_face_vertices_id.at(i).size();

    for (int j = 0; j < nb_face_vertices; ++j) {

      const int vertice_id = m_face_vertices_id.at(i).at(j) - 1;
      assert((0 <= vertice_id) && (vertice_id < m_nb_vertices));

      face_center_x += m_vertices_x[vertice_id];
      face_center_y += m_vertices_y[vertice_id];
      face_center_z += m_vertices_z[vertice_id];

    }

    face_center_x /= nb_face_vertices;
    face_center_y /= nb_face_vertices;
    face_center_z /= nb_face_vertices;
    
    face_centers_x[i] = face_center_x;
    face_centers_y[i] = face_center_y;
    face_centers_z[i] = face_center_z;

    RealType normal_x = 0.0;
    RealType normal_y = 0.0;
    RealType normal_z = 0.0;

    for (int j = 0; j < nb_face_vertices; ++j) {

      const int id = m_face_vertices_id.at(i).at(j) - 1;
      assert((0 <= id) && (id < m_nb_vertices));

      const int next_id = 
	m_face_vertices_id.at(i).at((j + 1) % nb_face_vertices) - 1;
      assert((0 <= id) && (id < m_nb_vertices));

      normal_x += 
	((m_vertices_y.at(id) - m_vertices_y.at(next_id)) * 
	 (m_vertices_z.at(id) + m_vertices_z.at(next_id)));

      normal_y += 
	((m_vertices_z.at(id) - m_vertices_z.at(next_id)) * 
	 (m_vertices_x.at(id) + m_vertices_x.at(next_id)));

      normal_z +=
	((m_vertices_x.at(id) - m_vertices_x.at(next_id)) * 
	 (m_vertices_y.at(id) + m_vertices_y.at(next_id)));

    }

    const RealType norm = 
      (normal_x * normal_x) + (normal_y * normal_y) + (normal_z * normal_z);

    normal_x /= norm;
    normal_y /= norm;
    normal_z /= norm;

    const int front_cell_id = m_face_front_cell_id.at(i);

    const RealType orienter =
      (face_centers_x[i] - cell_centers_x[front_cell_id]) * normal_x +
      (face_centers_y[i] - cell_centers_y[front_cell_id]) * normal_y +
      (face_centers_z[i] - cell_centers_z[front_cell_id]) * normal_z;

    if (orienter <= 0.0) {

      normal_x *= -1;
      normal_y *= -1;
      normal_z *= -1;

    }

    normals_x[i] = normal_x;
    normals_y[i] = normal_y;
    normals_y[i] = normal_z;

    // Face area
    const int vertice1 = m_face_vertices_id.at(i)[0] - 1;
    assert((0 <= vertice1) && (vertice1 < m_nb_vertices));
    
    const RealType x1 = vertices_x[vertice1];
    const RealType y1 = vertices_y[vertice1];
    const RealType z1 = vertices_z[vertice1];
    
    const int vertice2 = m_face_vertices_id.at(i)[1] - 1;
    assert((0 <= vertice2) && (vertice2 < m_nb_vertices));

    const RealType x2 = vertices_x[vertice2];
    const RealType y2 = vertices_y[vertice2];
    const RealType z2 = vertices_z[vertice2];

    RealType face_area = -1.0;

    if (m_face_vertices_id.at(i).size() == TRIANGLE) {

      face_area =
	std::sqrt((x2 - x1) * (x2 - x1) + 
		  (y2 - y1) * (y2 - y1) +
		  (z2 - z1) * (z2 - z1));

    }

    else if (m_face_vertices_id.at(i).size() == TETRA) {

      int vertice3 = m_face_vertices_id.at(i)[1] - 1;
      assert((0 <= vertice3) && (vertice3 < m_nb_vertices));

      const RealType x3 = vertices_x[vertice3];
      const RealType y3 = vertices_y[vertice3];
      const RealType z3 = vertices_z[vertice3];

      face_area =
	0.5 * fabs(Orient2d(x1, y1, z1, x2, y2, z2, x3, y3, z3));

    }



    face_volumes[i] = face_area;

  }

  std::cerr << "done.\n\n";
}

std::string Mesh::filename() {
  return m_filename;
}

void Mesh::Load(ptree &pt) {
  
  m_filename = pt.get<std::string>("filename", "");

}

void Mesh::Save(ptree &pt) {

  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree mesh_ptree;

  const int checkpoint = 0;

  if (!checkpoint) {

    pt.put("filename", m_filename);

  }

  else {

    //pt.put_child("Points", mesh_ptree);
    
    

  }
  
}

/// Computes a reordering pattern for mesh facets so that boundary
/// facets come first.
void Mesh::ComputeFacesMap() {

  m_faces_map = std::vector<int>(m_nb_faces, -1);

  int cnt_boundary_faces = 0;
  int cnt_interior_faces = 0;

  for (int i = 0; i < m_nb_faces; ++i) {

    assert(0 <= m_face_front_cell_id[i]);

    // Boundary face. Put it in the beginning.
    if (m_face_back_cell_id[i] == -1) {

      m_faces_map[i] = cnt_boundary_faces;
      cnt_boundary_faces += 1;
      
    }

    // Interior face. Put it at the end.
    else {
      m_faces_map[i] = m_nb_faces - 1 - cnt_interior_faces;
      cnt_interior_faces += 1;
    }

  }

  m_nb_boundary_faces = cnt_boundary_faces;
  m_nb_interior_faces = cnt_interior_faces;

  std::cerr << "boundary faces " << m_nb_boundary_faces << 
    " interior faces " << m_nb_interior_faces << "\n";
}

void Mesh::ReorderCells(VariableStore* var) {

  assert(m_nb_cells == var->nb_elements());

  RealType* cell_centers_x = (*var)(CELL_CENTERS_X);
  RealType* cell_centers_y = (*var)(CELL_CENTERS_Y);
  RealType* cell_volumes = (*var)(CELL_VOLUMES);

  for (int i = 0; i < m_nb_cells; ++i) {

    const int new_cell_id = m_cells_map.at(i);
    
    // swap for all mesh cell values.

    std::swap(cell_centers_x[i], cell_centers_x[new_cell_id]);
    std::swap(cell_centers_y[i], cell_centers_y[new_cell_id]);
    
    for (int id_neighbor = 0; id_neighbor < 3; ++id_neighbor)
      std::swap(m_cell_cell_connectivity.at(id_neighbor).at(i),
		m_cell_cell_connectivity.at(id_neighbor).at(new_cell_id));

  }

}

void Mesh::ReorderFaces(VariableStore* var) {

  // N.B. Swapping only won't cut it...

  RealType* normals_x = (*var)(NORMALS_X);
  RealType* normals_y = (*var)(NORMALS_Y);

  RealType* face_centers_x = (*var)(FACE_CENTERS_X);
  RealType* face_centers_y = (*var)(FACE_CENTERS_Y);
  RealType* face_volumes = (*var)(FACE_VOLUMES);

  std::vector<int> temp_face_front_cell_id(m_nb_faces, -1);
  std::vector<int> temp_face_back_cell_id(m_nb_faces, -1);

  std::vector<RealType> temp_face_centers_x(m_nb_faces, 0.0);
  std::vector<RealType> temp_face_centers_y(m_nb_faces, 0.0);

  std::vector<RealType> temp_face_volumes(m_nb_faces, 0.0);
  
  std::vector<RealType> temp_face_normals_x(m_nb_faces, 0.0);
  std::vector<RealType> temp_face_normals_y(m_nb_faces, 0.0);

  for (int i = 0; i < m_nb_faces; ++i) {

    const int new_face_id = m_faces_map.at(i);

    assert((0 <= new_face_id) && (new_face_id < m_nb_faces));

    temp_face_front_cell_id.at(new_face_id) =
      m_face_front_cell_id.at(i);

    temp_face_back_cell_id.at(new_face_id) = 
      m_face_back_cell_id.at(i);

    temp_face_centers_x.at(new_face_id) = 
      face_centers_x[i];

    temp_face_centers_y.at(new_face_id) = 
      face_centers_y[i];

    temp_face_volumes.at(new_face_id) = 
      face_volumes[i];

    temp_face_normals_x.at(new_face_id) = 
      normals_x[i];

    temp_face_normals_y.at(new_face_id) = 
      normals_y[i];

  }

  for (int i = 0; i < m_nb_faces; ++i) {

    m_face_front_cell_id.at(i) = temp_face_front_cell_id.at(i);
    m_face_back_cell_id.at(i) = temp_face_back_cell_id.at(i);
    
    face_centers_x[i] = temp_face_centers_x.at(i);
    face_centers_y[i] = temp_face_centers_y.at(i);

    face_volumes[i] = temp_face_volumes.at(i);

    normals_x[i] = temp_face_normals_x.at(i);
    normals_y[i] = temp_face_normals_y.at(i);

  }

}

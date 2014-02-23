#include "structured_grid.hpp"

#include <cassert>
#include <fstream>
#include <iostream>

#include "variable_metadata.hpp"
#include "variable_store.hpp"

static void SampleArray(double value_first, 
			double value_last, 
			int resolution,
			std::vector<double>* array_ptr) {

  assert(value_first <= value_last);
  assert(0 < resolution);

  array_ptr->resize(resolution);

  for (int i = 0; i < resolution; ++i)
    array_ptr->at(i) = value_first + 
      (double)i * (value_last - value_first) / (double)resolution;

}

StructuredGrid::StructuredGrid():
  m_nx(0), m_ny(0), m_nz(0),
  m_xmin(0.0), m_xmax(0.0),
  m_ymin(0.0), m_ymax(0.0),
  m_zmin(0.0), m_zmax(0.0), m_dimension(0) {};

StructuredGrid::StructuredGrid(double xmin, double xmax, int nx,
			       double ymin, double ymax, int ny):
  m_nx(nx), m_ny(ny), m_nz(0),
  m_xmin(xmin), m_xmax(xmax),
  m_ymin(ymin), m_ymax(ymax),
  m_zmin(0.0), m_zmax(0.0), 
  m_dimension(2) {

  SampleArray(m_xmin, m_xmax, m_nx + 1, &m_x_coordinates);
  SampleArray(m_ymin, m_ymax, m_ny + 1, &m_y_coordinates);
  SampleArray(m_zmin, m_zmax, m_nz + 1, &m_z_coordinates);

}

StructuredGrid::StructuredGrid(double xmin, double xmax, int nx,
			       double ymin, double ymax, int ny,
			       double zmin, double zmax, int nz):
  m_nx(nx), m_ny(ny), m_nz(nz),
  m_xmin(xmin), m_xmax(xmax),
  m_ymin(ymin), m_ymax(ymax),
  m_zmin(zmin), m_zmax(zmax), 
  m_dimension(3) {

  SampleArray(m_xmin, m_xmax, m_nx + 1, &m_x_coordinates);
  SampleArray(m_ymin, m_ymax, m_ny + 1, &m_y_coordinates);
  SampleArray(m_zmin, m_zmax, m_nz + 1, &m_z_coordinates);

};

void StructuredGrid::WriteHeaderVTKXmlAscii(std::ofstream* os_ptr) {

    *os_ptr << "<?xml version=\"1.0\"?>\n"
	    << "<VTKFile type=\"RectilinearGrid\" "
	    << "version=\"0.1\" "
	    << "byte_order=\"LittleEndian\">\n";

    *os_ptr << "<RectilinearGrid WholeExtent=\""
	    << 0 << " " << nx() << " "
	    << 0 << " " << ny() << " "
      	    << 0 << " " << nz()
	    << "\">\n";

    *os_ptr << "<Piece Extent=\""
	    << 0 << " " << nx() << " "
	    << 0 << " " << ny() << " "
	    << 0 << " " << nz() << "\""
	    << " Dimensions=\"" 
	    << nx() + 1 << " " << ny() + 1 << " " << nz() + 1 << "\""
	    << ">\n";


}

void StructuredGrid::WriteFooterVTKXmlAscii(std::ofstream* os_ptr) {

  *os_ptr << "</Piece>\n"
	  << "</RectilinearGrid>\n";
  
  *os_ptr << "</VTKFile>\n";

}


void StructuredGrid::WriteVTKXmlAscii(std::ofstream* os_ptr) {

  *os_ptr << "<Coordinates>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_x = m_x_coordinates.size();

  for (int i = 0; i < nb_points_x; ++i)
    *os_ptr << m_x_coordinates.at(i) << "\n";

  *os_ptr << "</DataArray>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_y = m_y_coordinates.size();

  for (int i = 0; i < nb_points_y; ++i)
    *os_ptr << m_y_coordinates.at(i) << "\n";

  *os_ptr << "</DataArray>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_z = m_z_coordinates.size();
  
  for (int i = 0; i < nb_points_z; ++i)
    *os_ptr << m_z_coordinates.at(i) << "\n";
  
  *os_ptr << "</DataArray>\n";

  *os_ptr << "</Coordinates>\n";

}

void StructuredGrid::ComputeGeometricQuantities(VariableStore* cell_variables, 
						VariableStore* face_variables,
						VariableStore* vertice_variables) {

  RealType* vertices_x = (*vertice_variables)(VERTICES_X);
  RealType* vertices_y = (*vertice_variables)(VERTICES_Y);
  RealType* vertices_z = (*vertice_variables)(VERTICES_Z);

  RealType* cell_centers_x = (*cell_variables)(CELL_CENTERS_X);
  RealType* cell_centers_y = (*cell_variables)(CELL_CENTERS_Y);
  RealType* cell_centers_z = (*cell_variables)(CELL_CENTERS_Z);
  RealType* cell_volumes = (*cell_variables)(CELL_VOLUMES);

  const int NX = nx();
  const int NY = ny();
  const int NZ = nz();
  
  if (m_dimension == 2) {

    for (int j = 0; j < NY + 1; ++j) {
      for (int i = 0; i < NX + 1; ++i) {

	const int local_id = (NX * j) + i;

	vertices_x[local_id] = m_x_coordinates[i];
	vertices_y[local_id] = m_y_coordinates[j];

      }
    }    

    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

	const int local_id = (NX * j) + i;
	
	cell_centers_x[local_id] = 0.5 * (m_x_coordinates[i + 1] + m_x_coordinates[i]);
	cell_centers_y[local_id] = 0.5 * (m_y_coordinates[j + 1] + m_y_coordinates[j]);

	cell_volumes[local_id] = 
	  (m_x_coordinates[i + 1] - m_x_coordinates[i]) * 
	  (m_y_coordinates[j + 1] - m_y_coordinates[j]);

      }
    }
      
  } else if (m_dimension == 3) {

    for (int k = 0; k < NZ + 1; ++k) {
      for (int j = 0; j < NY + 1; ++j) {
	for (int i = 0; i < NX + 1; ++i) {

	  const int local_id = (NX * j) + i;
	  
	  vertices_x[local_id] = m_x_coordinates[i];
	  vertices_y[local_id] = m_y_coordinates[j];
	  vertices_y[local_id] = m_z_coordinates[k];

	}
      }    
    }
    
    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
	for (int i = 0; i < NX; ++i) {

	  const int local_id = (NX * NY * k) + (NX * j) + i;
	
	  cell_centers_x[local_id] = 0.5 * (m_x_coordinates[i + 1] + m_x_coordinates[i]);
	  cell_centers_y[local_id] = 0.5 * (m_y_coordinates[j + 1] + m_y_coordinates[j]);
	  cell_centers_z[local_id] = 0.5 * (m_z_coordinates[k + 1] + m_z_coordinates[k]);

	  cell_volumes[local_id] = 
	    (m_x_coordinates[i + 1] - m_x_coordinates[i]) * 
	    (m_y_coordinates[j + 1] - m_y_coordinates[j]) * 
	    (m_z_coordinates[k + 1] - m_z_coordinates[k]);

	}
      }
    }
  }

  RealType* face_centers_x = (*face_variables)(FACE_CENTERS_X);
  RealType* face_centers_y = (*face_variables)(FACE_CENTERS_Y);
  RealType* face_centers_z = (*face_variables)(FACE_CENTERS_Z);

  if (m_dimension == 2) {

    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX + 1; ++i) {

	const int local_id = (NX * j) + i;
	
	face_centers_x[local_id] = m_x_coordinates[i];
	face_centers_y[local_id] = 0.5 * (m_y_coordinates[j + 1] + m_y_coordinates[j]);

      }
    }
      
  } else if (m_dimension == 3) {
    
    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
	for (int i = 0; i < NX; ++i) {

	  const int local_id = (NX * NY * k) + (NX * j) + i;

	  face_centers_x[local_id] = m_x_coordinates[i];
	  face_centers_y[local_id] = 0.5 * (m_y_coordinates[j + 1] + m_y_coordinates[j]);
	  face_centers_z[local_id] = 0.5 * (m_z_coordinates[j + 1] + m_z_coordinates[j]);
	
	}
      }
    }
  }



}

void StructuredGrid::Load(const ptree& pt) {

  m_nx = pt.get<int>("nx", 0);
  m_ny = pt.get<int>("ny", 0);

  if (pt.count("nz") == 0) {

    m_dimension = 2;
    m_nz = 0;

  } else {

    m_dimension = 3;
    m_nz = pt.get<int>("nz", 0);

  }

  m_xmin = pt.get<double>("xmin", 0.0);
  m_xmax = pt.get<double>("xmax", 0.0);
  m_ymin = pt.get<double>("ymin", 0.0);
  m_ymax = pt.get<double>("ymax", 0.0);
  m_zmin = pt.get<double>("zmin", 0.0);
  m_zmax = pt.get<double>("zmax", 0.0);

  SampleArray(m_xmin, m_xmax, m_nx + 1, &m_x_coordinates);
  SampleArray(m_ymin, m_ymax, m_ny + 1, &m_y_coordinates);
  SampleArray(m_ymin, m_ymax, m_nz + 1, &m_z_coordinates);
}

void Save(ptree& pt) {

  assert(0);

}

int StructuredGrid::nb_cells() { 
  return (m_dimension == 2 ? m_nx * m_ny : m_nx * m_ny * m_nz); }

int StructuredGrid::nb_vertices() { 
  return (m_dimension == 2 ? (m_nx + 1) * (m_ny + 1) : (m_nx + 1) * (m_ny + 1) * (m_nz + 1)); }

int StructuredGrid::nb_faces() { 
  return (m_dimension == 2 ? (m_nx + 1) * (m_ny + 1) : (m_nx + 1) * (m_ny + 1) * (m_nz + 1)); }

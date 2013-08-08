#ifndef STRUCTURED_GRID_HPP
#define STRUCTURED_GRID_HPP

#include <iosfwd>
#include <vector>

#include "grid.hpp"
#include "serialize.hpp"

class StructuredGrid : public Grid {
public:
  StructuredGrid();
  StructuredGrid(double xmin, double xmax, int nx,
		 double ymin, double ymax, int ny);
  StructuredGrid(double xmin, double xmax, int nx,
		 double ymin, double ymax, int ny,
		 double zmin, double zmax, int nz);
  virtual void WriteHeaderVTKXmlAscii(std::ofstream* os_ptr);
  virtual void WriteFooterVTKXmlAscii(std::ofstream* os_ptr);
  virtual void WriteVTKXmlAscii(std::ofstream* os_ptr);
  virtual void ComputeGeometricQuantities(VariableStore*, 
					  VariableStore*, 
					  VariableStore*);
  virtual int nb_cells();
  virtual int nb_vertices();
  virtual int nb_faces();
  double xmin() const { return m_xmin; }
  double xmax() const { return m_xmax; }
  double ymin() const { return m_ymin; }
  double ymax() const { return m_ymax; }
  double zmin() const { return m_zmin; }
  double zmax() const { return m_zmax; }
  int nx() const { return m_nx; }
  int ny() const { return m_ny; }
  int nz() const { return m_nz; }
  void Load(const ptree& pt);
  void Save(ptree&);
  virtual int dimension() { return m_dimension; }
private:
  std::vector<double> m_x_coordinates;
  std::vector<double> m_y_coordinates;
  std::vector<double> m_z_coordinates;
  int m_nx;
  int m_ny;
  int m_nz;
  double m_xmin;
  double m_xmax;
  double m_ymin;
  double m_ymax;
  double m_zmin;
  double m_zmax;
  int m_dimension;
};

#endif // STRUCTURED_GRID_HPP

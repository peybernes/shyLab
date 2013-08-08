#ifndef GRID_HPP
#define GRID_HPP

#include <iosfwd>

class VariableStore;

class Grid {
  virtual void WriteHeaderVTKXmlAscii(std::ofstream* os_ptr) = 0;
  virtual void WriteFooterVTKXmlAscii(std::ofstream* os_ptr) = 0;
  virtual void WriteVTKXmlAscii(std::ofstream* os_ptr) = 0;
  virtual void ComputeGeometricQuantities(VariableStore*, 
					  VariableStore*, 
					  VariableStore*) = 0;
  virtual int nb_cells() = 0;
  virtual int nb_vertices() = 0;
  virtual int nb_faces() = 0;
  virtual int dimension() = 0;
};

#endif // GRID_HPP

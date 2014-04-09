#ifndef GRID_HPP
#define GRID_HPP

#include <iosfwd>
#include <string>

class VariableStore;

class Grid {
public:
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
  inline std::string name() { return m_name;}
protected:
  std::string m_name;
private:
};

#endif // GRID_HPP

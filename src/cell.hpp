#ifndef CELL_HPP
#define CELL_HPP

#include <iosfwd>
#include <vector>

#include "mesh_elements.hpp"

/// Unstructured mesh object.
class Cell {
public:
  Cell();
  Cell(MeshElement cell_type, const std::vector<int>& vertices);
  Cell(MeshElement cell_type, const std::vector<int>& vertices, 
       const std::vector<int>& partition_ids);
  bool operator< (const Cell& other) const;
  bool operator== (const Cell& other) const;
  std::vector<Cell> ComputeFaceEntries();
  std::vector<int> vertices() const;
  std::vector<int> partition_ids() const;
  inline MeshElement cell_type() { return m_cell_type; }
  inline MeshElement cell_type() const { return m_cell_type; }
  friend std::size_t hash_value(const Cell& input);
  friend std::ostream& operator << (std::ostream&, Cell cell);
  int Dimension();
private:
  MeshElement m_cell_type;
  std::vector<int> m_vertices;
  std::vector<int> m_partition_ids;
};

#endif // CELL_HPP

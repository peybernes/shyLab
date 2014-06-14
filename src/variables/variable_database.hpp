#ifndef VARIABLE_DATABASE_HPP
#define VARIABLE_DATABASE_HPP

#include <map>
#include <string>

#include "variable_attribute.hpp"

typedef enum {NONE, VERTICE, FACET, CELL, PARTICLE} VariableSupport;

typedef enum {BOUNDARY, TOTAL} VariableRegion;

class VariableEntry {
public:
  VariableEntry();
  VariableEntry(int id, VariableSupport support,
		int attribute, VariableRegion region);
  inline bool written() { return m_attribute & WRITTEN; }
  inline bool written() const { return m_attribute & WRITTEN; }
  inline bool communicated() { return m_attribute & COMMUNICATED; }
  inline bool initializable() { return m_attribute & INITIALIZABLE; }
  inline bool checkpointed() { return m_attribute & PROTECTED; }
  inline int id() { return m_id; }
  inline int id() const { return m_id; }
  inline VariableSupport support() { return m_support; }
  inline VariableSupport support() const { return m_support; }
private:
  int m_id;
  VariableSupport m_support;
  int m_attribute;
  VariableRegion m_region;
};

typedef std::map<std::string, VariableEntry> VariableDatabase;

void AddVariableInDatabase(const std::string& variable_name, 
			   const VariableEntry& variable_entry, 
			   VariableDatabase* database_ptr);

void InitVariableDatabase(VariableDatabase* database_ptr);

int NumberOfVariableOfGivenSupport(VariableSupport support, const VariableDatabase& database);

int NumberOfCellVariables(const VariableDatabase& database);

int NumberOfVerticeVariables(const VariableDatabase& database);

int NumberOfFacetVariables(const VariableDatabase& database);


#endif // VARIABLE_DATABASE_HPP

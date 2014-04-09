#ifndef INIT_HPP
#define INIT_HPP

#include <string>

#include "event.hpp"
#include "serialize_fwd.hpp"
#include "variable_database.hpp"

// Row major layout : C/C++. Column major : Fortran.
typedef enum {ROW_MAJOR, COLUMN_MAJOR} LayoutType;

class InitVariable : public Event {
public:
  InitVariable(Simulation*, const Timetable&);
  virtual void Execute();
  virtual void Save(ptree&);
  virtual void Load(ptree&);
  std::string expression();
  const std::string expression() const;
  int id();
private:
  VariableEntry m_variable_entry;
  std::string m_expression;
  std::string m_stream_name;
  std::string m_var_name;
  LayoutType m_layout;
};

#endif // INIT_HPP

#ifndef INIT_HPP
#define INIT_HPP

#include <string>

#include "event.hpp"
#include "serialize_fwd.hpp"
#include "variable_database.hpp"

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
};

#endif // INIT_HPP

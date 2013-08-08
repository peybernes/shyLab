#ifndef MATH_PARSER_HPP
#define MATH_PARSER_HPP

#include <string>
#include <iosfwd>

#include <muParser.h>

#include "variable_database.hpp"

class NumericalParams;
class PhysicalParams;
class Simulation;

class MathParser : public mu::Parser {
public:
  MathParser();
  void RegisterSimulationVariables(const VariableDatabase& database);
  double SafeEval(); // Eval with error checking and reporting.
  void SetParameters(const PhysicalParams&, 
		     const NumericalParams&);
  void RegisterExpression(const std::string&);
  void UpdateReservedVariables(Simulation*, VariableSupport, int);
  void ListVariables(std::ostream&);
  void ListExpressions(std::ostream&);
  void RegisterVariables();
  std::map<std::string, mu::value_type> const& local_variables() const;
private:
  std::map<std::string, mu::value_type> m_local_variables;
};

#endif // MATH_PARSER_HPP

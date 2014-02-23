#include "math_parser.hpp"

#include <iostream>

#include <muParserDef.h>
#include <boost/foreach.hpp>

#include "variable_metadata.hpp"
#include "numerical_params.hpp"
#include "physical_params.hpp"
#include "simulation.hpp"

const int PARSER_VARIABLE_BUF_SIZE = 10000;

//! Factory function for creating new parser variables. This could as
// well be a function performing database queries.
static double* AddVariable(const mu::char_type *a_szName, 
			   void *a_pUserData) {
  // I don't want dynamic allocation here, so i used this static
  // buffer If you want dynamic allocation you must allocate all
  // variables dynamically in order to delete them later on. Or you
  // find other ways to keep track of variables that have been created
  // implicitely.
  static double afValBuf[PARSER_VARIABLE_BUF_SIZE];  
  static int iVal = 0;          

  mu::console() << _T("Math. parser : generating new variable \"") 
                << a_szName << "\"\n";

  afValBuf[iVal] = 0;

  if (iVal >= PARSER_VARIABLE_BUF_SIZE)
    throw mu::ParserError(_T("Variable buffer overflow."));

  return &afValBuf[iVal++];
}


MathParser::MathParser() {

  try {
    this->SetVarFactory(AddVariable, this);
  }

  catch(mu::Parser::exception_type &e) {
    mu::console() << _T("\nMathematical parser initialization error:\n");
    mu::console() << _T("Message:     ")   << e.GetMsg()   << _T("\n");
    mu::console() << _T("Expression:  \"") << e.GetExpr()  << _T("\"\n");
    mu::console() << _T("Token:       \"") << e.GetToken()    << _T("\"\n");
  }

  // Set local variables.

  m_local_variables["Pi"] = 3.1415926;

  m_local_variables["pi"] = 3.1415926;

  m_local_variables["PI"] = 3.1415926;

  // time.
  m_local_variables["t"] = 0.0;

  // Physical parameters.
  m_local_variables["g"] = 0.0;

  this->RegisterVariables();
}

void MathParser::RegisterSimulationVariables(const VariableDatabase& database) {

  typedef VariableDatabase::const_iterator IterType;

  for (IterType iter = database.begin();
       iter != database.end(); ++iter) {

    const std::string variable_name = iter->first;

    assert(database.count(variable_name) != 0);
    
    const VariableEntry variable_entry = iter->second;
    
    m_local_variables[variable_name] = 0.0;

  }  
  
  this->RegisterVariables();
}

void MathParser::SetParameters(const PhysicalParams &phy_params, 
			       const NumericalParams &num_params) {

  m_local_variables["g"] = phy_params.g;
  
}

void MathParser::RegisterVariables() {

  typedef std::pair<std::string, mu::value_type> var_entry_t;

  BOOST_FOREACH(var_entry_t var, m_local_variables) {

    DefineVar(var.first, &(m_local_variables[var.first]));

  }

  // Specific (abbreviated) names

}


double MathParser::SafeEval() {

  try {
    this->Eval();
  }
  
  catch(mu::Parser::exception_type &e) {
    mu::console() << _T("\nMath parsing error:\n");
    mu::console() << _T("-------------------\n");
    mu::console() << _T("Message:     ")   << e.GetMsg()   << _T("\n");
    mu::console() << _T("Expression:  \"") << e.GetExpr()  << _T("\"\n");
    mu::console() << _T("Token:       \"") << e.GetToken()    << _T("\"\n");

  }
}

void MathParser::RegisterExpression(const std::string &expr) {

  try {
    this->SetExpr(expr);
  }

  catch(mu::Parser::exception_type &e) {
    mu::console() << _T("\nMath parsing error:\n");
    mu::console() << _T("-------------------\n");
    mu::console() << _T("Message:     ")   << e.GetMsg()   << _T("\n");
    mu::console() << _T("Expression:  \"") << e.GetExpr()  << _T("\"\n");
    mu::console() << _T("Token:       \"") << e.GetToken()    << _T("\"\n");
    
 }

}

std::map<std::string, mu::value_type> const&
  MathParser::local_variables() const {
  return m_local_variables;
}

void MathParser::UpdateReservedVariables(Simulation* sim,
					 VariableSupport support, 
					 int i) {

  typedef VariableDatabase::const_iterator IterType;

  for (IterType iter = sim->variables_database.begin();
       iter != sim->variables_database.end(); ++iter) {

    const std::string variable_name = iter->first;

    assert(sim->variables_database.count(variable_name) != 0);
    
    const VariableEntry variable_entry = iter->second;

    if (variable_entry.support() == CELL) {
      
      m_local_variables[variable_name] = 
	sim->cell_variables(variable_entry.id())[i];

    } else if (variable_entry.support() == VERTICE) {

      m_local_variables[variable_name] = 
	sim->vertice_variables(variable_entry.id())[i];

    } else if (variable_entry.support() == FACET) {

      m_local_variables[variable_name] = 
	sim->face_variables(variable_entry.id())[i];

    } else {

      assert(0);

    }

  }

  m_local_variables["t"] = sim->clock.time();
}

void MathParser::ListVariables(std::ostream &os) {

  mu::varmap_type variables = this->GetVar();
  
  if (!variables.size())
    return;

  os << "Registered parser variables: (";

  mu::varmap_type::const_iterator item = variables.begin();
  for (; item!=variables.end(); ++item)
    mu::console() <<  _T(" ") << item->first;
  
  mu::console() << " )" << std::endl;

}

void MathParser::ListExpressions(std::ostream &os)
{
  mu::string_type sExpr = this->GetExpr();
  
  if (sExpr.length() == 0) {
    
    os << _T("Expression string is empty\n");
    return;
  }

  // Query the used variables (must be done after calc)
  mu::console() << _T("\nExpression variables:\n");
  mu::console() <<   _T("---------------------\n");
  mu::console() << _T("Expression: ") << this->GetExpr() << _T("\n");

  mu::varmap_type variables = this->GetUsedVar();

  if (!variables.size())
    mu::console() << _T("Expression does not contain variables\n");

  else {
    mu::console() << _T("Number: ") << (int)variables.size() << _T("\n");
    mu::varmap_type::const_iterator item = variables.begin();
    for (; item!=variables.end(); ++item)
      mu::console() << _T("Name: ") 
		    << item->first 
		    << _T("   Address: [0x") 
		    << item->second << _T("]\n");
  }
}

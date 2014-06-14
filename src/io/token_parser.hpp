#ifndef TOKEN_PARSER_HPP_
#define TOKEN_PARSER_HPP_

#include <string>
#include <sstream>

/// Convert a string into any type, returning 0 for invalid input
/// (N.B. as an alternative, we could use boost lexical_cast).
template <class T>
bool ParseToken(const std::string& s, T* t_ptr) {
  
  std::string garbage;
  std::istringstream iss(s);
  bool begins_good = !(iss >> *t_ptr).fail();
  iss.clear();
  iss.seekg(0, std::ios::beg);
  iss >> *t_ptr >> garbage;

  return begins_good && garbage.empty();

}

#endif // TOKEN_PARSER_HPP_

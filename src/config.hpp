#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include <limits>

const RealType INFTY = std::numeric_limits<RealType>::max();
const RealType MINUS_INFTY = std::numeric_limits<RealType>::min();

const int MAXINT = std::numeric_limits<int>::max();

/// A macro to disallow the copy constructor and operator= functions
/// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName)	\
  TypeName(const TypeName&);			\
  void operator=(const TypeName&)

#endif // CONFIG_HPP_

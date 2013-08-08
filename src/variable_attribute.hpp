#ifndef VARIABLE_ATTRIBUTE_HPP
#define VARIABLE_ATTRIBUTE_HPP

/// Attribute of a variable. See for instance
/// cell_variable_metadata.hpp for uses of this enum.
///
/// PROTECTED : the variable will be checkpointed.
///
/// WRITTEN : the variable will be dumped by #OutputSimulation (VTK
/// dump).
///
/// COMMUNICATED : the variable will be communicated using ghost cells
/// in MPI mode.
///
/// INITIALIZABLE : the variable can be initialized from the parameter
/// file.
enum VariableAttribute {
  PROTECTED = 1, // << 0
  WRITTEN = 1 << 1,
  COMMUNICATED = 1 << 2,
  INITIALIZABLE = 1 << 3
};

#endif // VARIABLE_ATTRIBUTE_HPP

#ifndef CELL_VARIABLE_METADATA_HPP
#define CELL_VARIABLE_METADATA_HPP

#include "variable_attribute.hpp"

enum {CELL_VOLUMES,
      CELL_CENTERS_X,
      CELL_CENTERS_Y,
      CELL_CENTERS_Z,
      NB_CELL_VALUES};

// Default variable attributes. Can be overwritten with the parameter
// file.
static const unsigned int cell_variable_attributes[NB_CELL_VALUES] =
  { !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_VOLUMES
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_X
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Y
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Z
  };

static const char* cell_variable_names[NB_CELL_VALUES] = {
  "cell_volumes",
  "cell_centers_x", 
  "cell_centers_y",
  "cell_centers_z",
};

enum {CELL_CENTERS,
      NB_CELL_VECTOR_VALUES};

const int cell_vector_variable_coordinates[NB_CELL_VECTOR_VALUES][3] =
  {{CELL_CENTERS_X, CELL_CENTERS_Y, CELL_CENTERS_Z}};

static const char* cell_vector_variable_names[NB_CELL_VECTOR_VALUES] =
  {"cell_centers"};

static const unsigned int cell_vector_variable_attributes[NB_CELL_VECTOR_VALUES] =
  {PROTECTED |  WRITTEN |  !COMMUNICATED};


#endif // CELL_VARIABLE_METADATA_HPP

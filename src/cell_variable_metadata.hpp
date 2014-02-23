#ifndef CELL_VARIABLE_METADATA_HPP
#define CELL_VARIABLE_METADATA_HPP

#include "variable_attribute.hpp"

enum {IN_RHO,
      IN_U,
      IN_V,
      IN_W,
      IN_E,
      IN_P,
      PREDICTED_RHO,
      PREDICTED_U,
      PREDICTED_V, 
      PREDICTED_W,
      PREDICTED_E,
      OUT_RHO,
      OUT_U,
      OUT_V, 
      OUT_W,
      OUT_E,
      OUT_P,
      OUT_TIMESTEP, 
      CELL_VOLUMES,
      CELL_CENTERS_X,
      CELL_CENTERS_Y,
      CELL_CENTERS_Z,
      RHO_REF,
      NB_CELL_VALUES};

// Default variable attributes. Can be overwritten with the parameter
// file.
static const unsigned int cell_variable_attributes[NB_CELL_VALUES] =
  { PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_RHO
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_U
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_V
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_W
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_E
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_P
    PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // PREDICTED_RHO
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_U
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_V
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_W
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_E
    PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // OUT_RHO
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_U
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_V
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_W
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_E
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_P
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // OUT_TIMESTEP
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_VOLUMES
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_X
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Y
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Z
    !PROTECTED |  WRITTEN | !COMMUNICATED | INITIALIZABLE // RHO_REF
  };

static const char* cell_variable_names[NB_CELL_VALUES] = {
  "in_rho",
  "in_u",
  "in_v",
  "in_w",
  "in_e",
  "in_p",
  "predicted_rho",
  "predicted_u",
  "predicted_v",
  "predicted_w",
  "predicted_e",
  "out_rho",
  "out_u",
  "out_v",
  "out_w",
  "out_e",
  "out_p",
  "local_timestep",
  "cell_volumes",
  "cell_centers_x", 
  "cell_centers_y",
  "cell_centers_z",
  "rho_ref"
};

enum {VELOCITY,
      NB_CELL_VECTOR_VALUES};

const int cell_vector_variable_coordinates[NB_CELL_VECTOR_VALUES][3] =
  {{OUT_U, OUT_V, OUT_W}};

static const char* cell_vector_variable_names[NB_CELL_VECTOR_VALUES] =
  {"Velocity"};

static const unsigned int cell_vector_variable_attributes[NB_CELL_VECTOR_VALUES] =
  {PROTECTED |  WRITTEN |  !COMMUNICATED};

#endif // CELL_VARIABLE_METADATA_HPP

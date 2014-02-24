#ifndef CELL_VARIABLE_METADATA_HPP
#define CELL_VARIABLE_METADATA_HPP

#include "variable_attribute.hpp"

enum {IN_RHO,
      IN_E,
      IN_P,
      IN_CELL_MASS,
      PREDICTED_RHO,
      PREDICTED_E,
      OUT_RHO,
      OUT_E,
      OUT_P,
      OUT_CELL_MASS,
      OUT_TIMESTEP, 
      CELL_VOLUMES,
      CELL_CENTERS_X,
      CELL_CENTERS_Y,
      CELL_CENTERS_Z,
      RHO_REF,
      P_REF,
      DIRECTIONAL_LAGRANGIAN_VOLUME,
      DIRECTIONAL_LAGRANGIAN_DENSITY,
      NB_CELL_VALUES};

// Default variable attributes. Can be overwritten with the parameter
// file.
static const unsigned int cell_variable_attributes[NB_CELL_VALUES] =
  { PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_RHO
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_E
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_P
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_CELL_MASS
    PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // PREDICTED_RHO
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_E
    PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // OUT_RHO
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_E
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_P
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // OUT_CELL_MASS
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // OUT_TIMESTEP
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_VOLUMES
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_X
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Y
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, // CELL_CENTERS_Z
    !PROTECTED |  WRITTEN | !COMMUNICATED | INITIALIZABLE, // RHO_REF
    !PROTECTED |  WRITTEN | !COMMUNICATED | INITIALIZABLE, // P_REF    
    !PROTECTED |  WRITTEN | !COMMUNICATED | INITIALIZABLE, // DIRECTIONAL_LAGRANGIAN_VOLUME
    !PROTECTED |  WRITTEN | !COMMUNICATED | INITIALIZABLE, // DIRECTIONAL_LAGRANGIAN_DENSITY    
  };

static const char* cell_variable_names[NB_CELL_VALUES] = {
  "in_rho",
  "in_e",
  "in_p",
  "in_cell_mass",
  "predicted_rho",
  "predicted_e",
  "out_rho",
  "out_e",
  "out_p",
  "out_cell_mass",
  "local_timestep",
  "cell_volumes",
  "cell_centers_x", 
  "cell_centers_y",
  "cell_centers_z",
  "rho_ref",
  "p_ref",
  "directional_lagrangian_volume",
  "directional_lagrangian_density"
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

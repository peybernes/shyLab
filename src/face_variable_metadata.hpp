#ifndef FACE_VARIABLE_METADATA_HPP
#define FACE_VARIABLE_METADATA_HPP

#include "variable_attribute.hpp"

enum {FLUX_RHO,
      FLUX_RHO_U,
      FLUX_RHO_V,
      FLUX_RHO_W,
      FLUX_RHO_E, 
      MAX_EIGEN, 
      FACE_VOLUMES,
      NORMALS_X,
      NORMALS_Y,
      NORMALS_Z,
      FACE_CENTERS_X,
      FACE_CENTERS_Y,
      FACE_CENTERS_Z,
      U_FACE,
      V_FACE,
      W_FACE,
      IN_RHO_FACE,
      OUT_RHO_FACE,
      NB_FACE_VALUES};

static const char* face_variable_names[NB_FACE_VALUES] = {
  "flux_rho",
  "flux_rho_u",
  "flux_rho_v",
  "flux_rho_w",
  "flux_rho_e",
  "max_eigen",
  "face_volumes",
  "normals_x",
  "normals_y",
  "normals_z", 
  "face_centers_x",
  "face_centers_y",
  "face_centers_z",
  "u_face",
  "v_face",
  "w_face",
  "in_rho_face",
  "out_rho_face"
};

// Default variable attributes. Can be overwritten with the parameter
// file.
static const unsigned int face_variable_attributes[NB_FACE_VALUES] =
  { PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, 
    PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE, 
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE,
    !PROTECTED |  WRITTEN | !COMMUNICATED | !INITIALIZABLE
  };


enum {NORMALS,
      NB_FACE_VECTOR_VALUES};

const int face_vector_variable_coordinates[NB_FACE_VECTOR_VALUES][3] =
  {{NORMALS_X, NORMALS_Y, NORMALS_Z}};

static const char* face_vector_variable_names[NB_FACE_VECTOR_VALUES] =
  {"Normals"};

#endif // FACE_VARIABLE_METADATA_HPP

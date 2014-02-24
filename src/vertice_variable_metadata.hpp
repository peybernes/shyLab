#ifndef VERTICE_VARIABLE_METADATA_HPP
#define VERTICE_VARIABLE_METADATA_HPP

enum {VERTICES_X,
      VERTICES_Y,
      VERTICES_Z,
      IN_U,
      IN_V,
      PREDICTED_U,
      PREDICTED_V, 
      OUT_U,
      OUT_V, 
      OUT_W,
      U_REF,
      NB_VERTICE_VALUES};

static const unsigned int vertice_variable_attributes[NB_VERTICE_VALUES] =
  { PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE,
    PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE,
    PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE,
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_U
    PROTECTED |  WRITTEN | !COMMUNICATED |  INITIALIZABLE, // IN_V
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_U
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // PREDICTED_V
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_U
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_V
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // OUT_W
    !PROTECTED |  WRITTEN |  COMMUNICATED | !INITIALIZABLE, // U_REF
  };

static const char* vertice_variable_names[NB_VERTICE_VALUES] =
  {"vertices_x",
   "vertices_y",
   "vertices_z",
   "in_u",
   "in_v",
   "predicted_u",
   "predicted_v",
   "out_u",
   "out_v",
   "out_w",
   "u_ref"
  };

enum {VELOCITY,
      NB_VERTICE_VECTOR_VALUES};

const int vertice_vector_variable_coordinates[NB_VERTICE_VECTOR_VALUES][3] =
  {{OUT_U, OUT_V, OUT_W}};

static const char* vertice_vector_variable_names[NB_VERTICE_VECTOR_VALUES] =
  {"Velocity"};

static const unsigned int vertice_vector_variable_attributes[NB_VERTICE_VECTOR_VALUES] =
  {PROTECTED |  WRITTEN |  !COMMUNICATED};

#endif // VERTICE_VARIABLE_METADATA_HPP

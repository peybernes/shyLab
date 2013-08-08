#ifndef VERTICE_VARIABLE_METADATA_HPP
#define VERTICE_VARIABLE_METADATA_HPP

enum {VERTICES_X,
      VERTICES_Y,
      VERTICES_Z,
      NB_VERTICE_VALUES};

static const unsigned int vertice_variable_attributes[NB_VERTICE_VALUES] =
  { PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE,
    PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE,
    PROTECTED |  !WRITTEN | !COMMUNICATED |  !INITIALIZABLE
  };

static const char* vertice_variable_names[NB_VERTICE_VALUES] =
  {"vertices_x",
   "vertices_y",
   "vertices_z"};

#endif // VERTICE_VARIABLE_METADATA_HPP


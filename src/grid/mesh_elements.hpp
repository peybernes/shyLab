#ifndef MESH_ELEMENTS_HPP
#define MESH_ELEMENTS_HPP

const int ambient_dimension = 3;

enum {
  GMSH_POINT = 15,
  GMSH_SEGMENT = 1,
  GMSH_TRIANGLE = 2,
  GMSH_QUAD = 3,
  GMSH_TETRA = 4,
  GMSH_HEXA = 5,
  GMSH_PRISM = 6,
  GMSH_PYRAMID = 7
};

enum {
  VTK_VERTEX = 1,
  VTK_LINE = 3,
  VTK_TRIANGLE = 5,
  VTK_QUAD = 9,
  VTK_TETRA = 10,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14
};

typedef enum {
  SEGMENT, TRIANGLE, QUAD, TETRA, HEXA, PRISM, PYRAMID, OTHER} MeshElement;

const int DUMMY = -1;

const int nb_faces_segment = 2;
const int taille_faces_segment = 1;

static int indices_segment[nb_faces_segment][taille_faces_segment] =
  {{0},
   {1}};

const int nb_faces_tri3 = 3;
const int taille_faces_tri3 = 2;

static int indices_tri3[nb_faces_tri3][taille_faces_tri3] =
  {{0, 1},
   {1, 2},
   {2, 0}};

const int nb_faces_quad4 = 4;
const int taille_faces_quad4 = 2;

static int indices_quad4[nb_faces_quad4][taille_faces_quad4] =
  {{0, 1},
   {1, 2},
   {2, 3},
   {3, 0}};

const int nb_faces_tetra4 = 4;
const int taille_faces_tetra4 = 3;

static int indices_tetra4[nb_faces_tetra4][taille_faces_tetra4] =
  {{0, 1, 3},
   {0, 1, 3},
   {1, 2, 3},
   {0, 2, 3}};

const int nb_faces_hexa8 = 6;
const int taille_faces_hexa8 = 4;

static int indices_hexa8[nb_faces_hexa8][taille_faces_hexa8] =
  {{0, 1, 2, 3},
   {4, 5, 6, 7},
   {0, 1, 5, 4},
   {3, 2, 6, 7},
   {1, 2, 6, 5},
   {0, 3, 7, 4}};

const int nb_faces_prism6 = 5;
const int taille_faces_prism6 = 4;

static int indices_prism6[nb_faces_prism6][taille_faces_prism6] =
  {{0, 1, 4, 3},
   {0, 1, 2, DUMMY},
   {1, 4, 5, 2},
   {3, 4, 5, DUMMY},
   {0, 1, 5, 2}};

const int nb_faces_pyr5 = 5;
const int taille_faces_pyr5 = 4;

static int indices_pyr5[nb_faces_pyr5][taille_faces_pyr5] =
  {{0, 1, 2, 3},
   {0, 1, 4, DUMMY},
   {1, 2, 4, DUMMY},
   {2, 3, 4, DUMMY},
   {3, 0, 4, DUMMY}};

#endif

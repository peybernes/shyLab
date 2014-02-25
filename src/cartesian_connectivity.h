typedef int index_t;

inline index_t CellCellP1O(index_t cell_oo, index_t nx) {

  return cell_oo + 1;

}

inline index_t CellCellM1O(index_t cell_oo, index_t nx) {

  return cell_oo - 1;

}

inline index_t CellCellOP1(index_t cell_oo, index_t nx) {

  return cell_oo + nx;

}

inline index_t CellCellOM1(index_t cell_oo, index_t nx) {

  return cell_oo - nx;

}

inline index_t CellFaceP1O(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy + 1;

}

inline index_t CellFaceM1O(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy;

}

inline index_t CellFaceOM1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo;

}

inline index_t CellFaceOP1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + nx;

}

inline index_t CellNodeM1M1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy;

}

inline index_t CellNodeP1M1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy + 1;

}

inline index_t CellNodeM1P1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy + nx + 1;

}

inline index_t CellNodeP1P1(index_t cell_oo, index_t iy, index_t nx) {

  return cell_oo + iy + nx + 2;

}

inline index_t FaceXCellP1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo - iy;

}

inline index_t FaceXCellM1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo - iy - 1;

}

inline index_t FaceYCellP1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo;

}

inline index_t FaceYCellM1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo - nx;

}

inline index_t FaceXNodeM1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo;

}

inline index_t FaceXNodeP1(index_t face_oo, index_t iy, index_t nx) {
  
  return face_oo + nx + 1;

}

inline index_t FaceYNodeM1(index_t face_oo, index_t iy, index_t nx) {

  return face_oo + iy;

}

inline index_t FaceYNodeP1(index_t face_oo, index_t iy, index_t nx) {
  
  return face_oo + iy + 1;

}

inline index_t NodeCellM1M1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo - (nx + 1 + iy);

}

inline index_t NodeCellP1M1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo -(nx + iy);

}

inline index_t NodeCellM1P1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo -(iy + 1);

}

inline index_t NodeCellP1P1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo - iy;

}

inline index_t NodeNodeM1O(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo - 1;

}

inline index_t NodeNodeP1O(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo + 1;

}
inline index_t NodeNodeOM1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo - (nx+ 1);

}

inline index_t NodeNodeOP1(index_t node_oo, index_t iy, index_t nx) {
  
  return node_oo + (nx + 1);

}

inline index_t NodeFaceXOOM1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - (nx + 1);

}

inline index_t NodeFaceXOOP1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo;

}

inline index_t NodeFaceXM1M1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - nx - 2;

}

inline index_t NodeFaceXM1P1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - 1;

}

inline index_t NodeFaceXP1M1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - nx;

}

inline index_t NodeFaceXP1P1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo + 1;

}

inline index_t NodeFaceYM1OO(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - iy - 1;

}

inline index_t NodeFaceYP1OO(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - iy;

}

inline index_t NodeFaceYM1M1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - nx - iy - 1;

}

inline index_t NodeFaceYM1P1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo + nx - iy - 1;

}

inline index_t NodeFaceYP1M1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo - nx - iy;

}

inline index_t NodeFaceYP1P1(index_t node_oo, index_t iy, index_t nx) {

  return node_oo + nx - iy;

}

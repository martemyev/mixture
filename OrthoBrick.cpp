#include "OrthoBrick.h"

OrthoBrick::OrthoBrick() {
  nVertices = 8;
  vertices = new Node3D[nVertices];
  cpointsFileName = "cpoints/cpointsBrickSurf.msh";
  templateName = "OrthoBrickTemplate";
}

OrthoBrick::~OrthoBrick() {
  delete[] vertices;
}

// initialization of the vertices
void OrthoBrick::verticesInit() {
  // NOTE! lengths[0] - length of side parallel to x-axis
  //       lengths[1] - length of side parallel to y-axis
  //       lengths[2] - length of side parallel to z-axis
  // at first, we create standard orthobrick with center at origin of coordinates
  double coord[][4] = { { -0.5 * lengths[0], -0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0], -0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        { -0.5 * lengths[0], -0.5 * lengths[1],  0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0], -0.5 * lengths[1],  0.5 * lengths[2], 1.0 }, \
                        { -0.5 * lengths[0],  0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0],  0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        { -0.5 * lengths[0],  0.5 * lengths[1],  0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0],  0.5 * lengths[1],  0.5 * lengths[2], 1.0 } };
  double coord_new[4];
  for (int i = 0; i < nVertices; i++) {
    toNewCoord(coord[i], coord_new); // from standard orthobrick at origin to orthobrick with rotation and translation
    vertices[i].init(coord_new, i);
  }
  
#ifdef DEBUG
  std::cout << "OrthoBrick : \n";
  for (int i = 0; i < nVertices; i++)
    std::cout << vertices[i].getX() << " " << vertices[i].getY() << " " << vertices[i].getZ() << std::endl;
#endif  
  
}

// calculate volume of the orthobrick
double OrthoBrick::volume() {
  return lengths[0] * lengths[1] * lengths[2];
}

// check the containing of the point
bool OrthoBrick::hasPoint(Node3D *point, double tolerance) {
  // at first we need to convert our system of coordinate
  // in such state, that orthobrick would with center at the origin of coordinates
  // and with edges parallel to coordinate axes.
  // to do that we need to use inverse transformation matrix
  double coord[] = { point->getX(), point->getY(), point->getZ(), 1.0 };
  double coord_origin[4];
  toOldCoord(coord, coord_origin); // transform to 'old' system of coordinates
  // now we can check the containing of the point with orthobrick in standard position
  // instead of orthobrick that has a rotation and transation in space.
  double px = coord_origin[0]; // coordinates of point
  double py = coord_origin[1];
  double pz = coord_origin[2];
  double x0 = -0.5 * lengths[0]; // limits of standard orthobrick
  double y0 = -0.5 * lengths[1];
  double z0 = -0.5 * lengths[2];
  double x1 = -x0;
  double y1 = -y0;
  double z1 = -z0;
  
  // if point is inside orthobrick, then return true
  if (px < x1 - tolerance && \
      px > x0 + tolerance && \
      py < y1 - tolerance && \
      py > y0 + tolerance && \
      pz < z1 - tolerance && \
      pz > z0 + tolerance)
    return true; // orthobrick contains the point

  return false; // orthobrick doesn't contain the point
}

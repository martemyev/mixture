#include "Ellipsoid.h"

Ellipsoid::Ellipsoid() {
  nVertices = 6; // 2 vertices for every direction
  vertices = new Node3D[nVertices];
  cpointsFileName = "cpoints/cpointsSphereSurf.msh";
  templateName = "EllipsoidTemplate";
}

Ellipsoid::~Ellipsoid() {
  delete[] vertices;
}

// initialization of the vertices
void Ellipsoid::verticesInit() {
  // NOTE! lengths[0] - major diameter (axis)
  //       lengths[1] - first minor diameter (axis)
  //       lengths[2] - second minor diameter (axis)
  // at first, we create standard ellipsoid with center at origin of coordinates
  double coord[][4] = { { -0.5 * lengths[0], 0.0, 0.0, 1.0 }, \
                        {  0.5 * lengths[0], 0.0, 0.0, 1.0 }, \
                        { 0.0, -0.5 * lengths[1], 0.0, 1.0 }, \
                        { 0.0,  0.5 * lengths[1], 0.0, 1.0 }, \
                        { 0.0, 0.0, -0.5 * lengths[2], 1.0 }, \
                        { 0.0, 0.0,  0.5 * lengths[2], 1.0 } };
  double coord_new[4];
  for (int i = 0; i < nVertices; i++) {
    toNewCoord(coord[i], coord_new); // from standard ellipsoid at origin to ellipsoid with rotation and translation
    vertices[i].init(coord_new, i);
  }
}

// calculate volume of the ellipsoid
double Ellipsoid::volume() {
  double a = 0.5 * lengths[0]; // major semi-axis
  double b = 0.5 * lengths[1]; // first minor semi-axis
  double c = 0.5 * lengths[2]; // second minor semi-axis
  double vol = 4.0 / 3.0 * PI * a * b * c;
  return vol;
}

// check if ellipsoid contains some point with accuracy in tolerance
bool Ellipsoid::hasPoint(Node3D *point, double tolerance) {
  // at first we need to convert our system of coordinate
  // in such state, that ellipsoid would with center at the origin of coordinates
  // and with axes parallel to Cartesian coordinates axes
  // to do that we need to use inverse transformation matrix
  double coord[] = { point->getX(), point->getY(), point->getZ(), 1.0 };
  double coord_origin[4];
  toOldCoord(coord, coord_origin); // transform to 'old' system of coordinates
  double a = 0.5 * lengths[0];
  double b = 0.5 * lengths[1];
  double c = 0.5 * lengths[2];
  double x = coord_origin[0];
  double y = coord_origin[1];
  double z = coord_origin[2];
  if (((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c)) < (1 - tolerance))
    return true; // ellipsoid contains the point
  
  return false; // ellipsoid doesn't contain the point
}

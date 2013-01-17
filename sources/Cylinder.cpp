#include "Cylinder.h"

Cylinder::Cylinder() {
  nVertices = 8; // 4 points for every base
  vertices = new Node3D[nVertices];
  cpointsFileName = "cpoints/cpointsCylinderSurf.msh";
  templateName = "CylinderTemplate";
}

Cylinder::~Cylinder() {
  delete[] vertices;
}

// initialization of the vertices
void Cylinder::verticesInit() {
  // NOTE! lengths[0] - major diameter of elliptic base
  //       lengths[1] - minor diameter of elliptic base
  //       lengths[2] - height of cylinder
  // at first, we create standard cylinder with center at origin of coordinates
  double coord[][4] = { { -0.5 * lengths[0], 0.0, -0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0], 0.0, -0.5 * lengths[2], 1.0 }, \
                        { 0.0, -0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        { 0.0,  0.5 * lengths[1], -0.5 * lengths[2], 1.0 }, \
                        { -0.5 * lengths[0], 0.0,  0.5 * lengths[2], 1.0 }, \
                        {  0.5 * lengths[0], 0.0,  0.5 * lengths[2], 1.0 }, \
                        { 0.0, -0.5 * lengths[1],  0.5 * lengths[2], 1.0 }, \
                        { 0.0,  0.5 * lengths[1],  0.5 * lengths[2], 1.0 } };
  double coord_new[4];
  for (int i = 0; i < nVertices; i++) {
    toNewCoord(coord[i], coord_new); // from standard cylinder at origin to cylinder with rotation and translation
    vertices[i].init(coord_new, i);
  }
}

// calculate volume of the cylinder
double Cylinder::volume() {
  double a = 0.5 * lengths[0]; // major semi-axis of base
  double b = 0.5 * lengths[1]; // minor semi-axis of base
  double vol = PI * a * b * lengths[2];
  return vol;
}

// check the containing of the point
bool Cylinder::hasPoint(Node3D *point, double tolerance) {
  // at first we need to convert our system of coordinate
  // in such state, that cylinder would with center at the origin of coordinates
  // and with height-axis parallel to z-coordinate axis.
  // to do that we need to use inverse transformation matrix
  double coord[] = { point->getX(), point->getY(), point->getZ(), 1.0 };
  double coord_origin[4];
  toOldCoord(coord, coord_origin); // transform to 'old' system of coordinates
  // now we can check the containing of the point with cylinder in standard position
  // instead of cylinder that has a rotation and transation in space.
  double x = coord_origin[0]; // coordinates of point
  double y = coord_origin[1];
  double z = coord_origin[2];
  double a = 0.5 * lengths[0]; // major semi-axis of the base
  double b = 0.5 * lengths[1]; // minor semi-axis of the base
  
  // at first, compare using height of cylinder
  if (fabs(z) < 0.5 * lengths[2] - tolerance) {
    // then, compare in elliptic cut of cylinder,
    // using this rule:
    // point (x, y) is inside ellipse with a and b semi-axes,
    // if equation holds: x^2 / a^2 + y^2 / b^2 < 1
    if ((x * x) / (a * a) + (y * y) / (b * b) < 1 - tolerance)
      return true; // cylinder contains the point
  }
  
  return false; // cylinder doesn't contain the point
}

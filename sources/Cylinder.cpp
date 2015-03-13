//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#include "Cylinder.h"
#include "Config.h"
#include "Require.h"
#include <fstream>

/**
 * Constructor of Cylinder class.
 * Cylinder is defined by 8 vertices (4 vertices on every base).
 * The base of cylinder could be ellipsoid (not only circle).
 */
Cylinder::Cylinder() {
  nVertices = 8; // 4 points for every base
  vertices = new Node3D[nVertices];
  cpointsFileName = (std::string)CPOINTS_DIR + "/cpointsCylinderSurf.msh";
  templateName = "CylinderTemplate";
}

/**
 * Destructor of Cylinder class.
 * Here we only delete vertices.
 */
Cylinder::~Cylinder() {
  delete[] vertices;
}

/**
 * Initialization of cylinder vertices.
 * lengths[0] - major diameter of elliptic base
 * lengths[1] - minor diameter of elliptic base
 * lengths[2] - height of cylinder
 */
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

/**
 * Calculate the volume of the cylinder.
 */
double Cylinder::volume() const {
  const double a = 0.5 * lengths[0]; // major semi-axis of base
  const double b = 0.5 * lengths[1]; // minor semi-axis of base
  return PI * a * b * lengths[2];
}

/**
 * Check the containing of the point.
 * \param[in] point - the point that we want to check
 * \param[in] tolerance - the accuracy of checking
 */
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
  const double x = coord_origin[0]; // coordinates of point
  const double y = coord_origin[1];
  const double z = coord_origin[2];
  const double a = 0.5 * lengths[0]; // major semi-axis of the base
  const double b = 0.5 * lengths[1]; // minor semi-axis of the base
  
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

/**
 * Create control points of unit cylinder: height is 1, base axes lengths are 1, center at origin
 */
void Cylinder::createUnitElementCPoints(double cl) {
  std::string procName = "Cylinder::createUnitElementCPoints";
  std::cout << "  " << procName << std::endl;
  std::string geoName = "cpointsCylSurf.geo";
  std::ofstream out(geoName.c_str());
  frequire(out, geoName, procName);
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n";
  out << "x0 = -0.5; y0 = 0; z0 = -0.5;\n";
  out << "x1 = 0.5; y1 = 0; z1 = -0.5;\n";
  out << "x2 = 0; y2 = -0.5; z2 = -0.5;\n";
  out << "x3 = 0; y3 = 0.5; z3 = -0.5;\n";
  out << "x4 = -0.5; y4 = 0; z4 = 0.5;\n";
  out << "x5 = 0.5; y5 = 0; z5 = 0.5;\n";
  out << "x6 = 0; y6 = -0.5; z6 = 0.5;\n";
  out << "x7 = 0; y7 = 0.5; z7 = 0.5;\n";
  out << "cl = " << cl << ";\n";
  out << "surfNumber = 1;\n";
  out << "volNumber = 0;\n";
  out << "Call " << templateName << ";\n";
  out << "Physical Surface(3101) = { cylinderSurfaces[] };\n";
  out.close();
  std::string createMesh = GMSH_BIN + GMSH_BUILD_CPOINTS_OPTIONS + geoName + " -o " + cpointsFileName;
  require(!system(createMesh.c_str()), "Mesh cannot be built!", procName);
}

#include "Ellipsoid.h"
#include "Config.h"
#include "Require.h"
#include <fstream>

/**
 * Constructor of Ellipsoid class.
 * Ellipsoid is defined by 6 vertices (2 vertices for every axis).
 */
Ellipsoid::Ellipsoid() {
  nVertices = 6; // 2 vertices for every direction
  vertices = new Node3D[nVertices];
  cpointsFileName = (std::string)CPOINTS_DIR + "/cpointsSphereSurf.msh";
  templateName = "EllipsoidTemplate";
}

/**
 * Destructor of Ellipsoid class.
 * Here we only delete vertices.
 */
Ellipsoid::~Ellipsoid() {
  delete[] vertices;
}

/**
 * Initialization of the vertices.
 * lengths[0] - major diameter (axis)
 * lengths[1] - first minor diameter (axis)
 * lengths[2] - second minor diameter (axis)
 */
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

/**
 * Calculate the volume of the ellipsoid.
 */
double Ellipsoid::volume() const {
  const double a = 0.5 * lengths[0]; // major semi-axis
  const double b = 0.5 * lengths[1]; // first minor semi-axis
  const double c = 0.5 * lengths[2]; // second minor semi-axis
  return 4.0 / 3.0 * PI * a * b * c;
}

/**
 * Check if the ellipsoid contains some point with accuracy in tolerance.
 * \param[in] point - the pointer to point that we want to check
 * \param[in] tolerance - the accuracy of checking
 */
bool Ellipsoid::hasPoint(Node3D *point, double tolerance) {
  // at first we need to convert our system of coordinate
  // in such state, that ellipsoid would with center at the origin of coordinates
  // and with axes parallel to Cartesian coordinates axes
  // to do that we need to use inverse transformation matrix
  double coord[] = { point->getX(), point->getY(), point->getZ(), 1.0 };
  double coord_origin[4];
  toOldCoord(coord, coord_origin); // transform to 'old' system of coordinates
  const double a = 0.5 * lengths[0];
  const double b = 0.5 * lengths[1];
  const double c = 0.5 * lengths[2];
  const double x = coord_origin[0];
  const double y = coord_origin[1];
  const double z = coord_origin[2];
  if (((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c)) < (1 - tolerance))
    return true; // ellipsoid contains the point
  
  return false; // ellipsoid doesn't contain the point
}

/**
 * Create control points of unit ellispoid: axes lengths are 1, center at origin
 */
void Ellipsoid::createUnitElementCPoints(double cl) {
  std::string procName = "Ellipsoid::createUnitElementCPoints";
  std::string geoName = "cpointsSphereSurf.geo";
  std::ofstream out(geoName.c_str());
  frequire(out, geoName, procName);
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n";
  out << "x0 = -0.5; y0 = 0; z0 = 0;\n";
  out << "x1 = 0.5; y1 = 0; z1 = 0;\n";
  out << "x2 = 0; y2 = -0.5; z2 = 0;\n";
  out << "x3 = 0; y3 = 0.5; z3 = 0;\n";
  out << "x4 = 0; y4 = 0; z4 = -0.5;\n";
  out << "x5 = 0; y5 = 0; z5 = 0.5;\n";
  out << "cl = " << cl << ";\n";
  out << "surfNumber = 1;\n";
  out << "volNumber = 0;\n";
  out << "Call " << templateName << ";\n";
  out << "Physical Surface(3101) = { ellSurfaces[] };\n";
  out.close();
  std::string createMesh = GMSH_BIN + GMSH_BUILD_CPOINTS_OPTIONS + geoName + " -o " + cpointsFileName;
  require(!system(createMesh.c_str()), "Mesh cannot be built!", procName);
}

/*
 * mixture - Copyright (c) 2012-2013 Mikhail Artemiev
 *
 * http://code.google.com/p/mixture
 *
 * This library is provided under the terms of MIT license.
 * See the LICENSE file for license information.
 */

#include "GeoTetrahedron.h"
#include "Config.h"
#include "Require.h"
#include <fstream>

/**
 * Constructor of GeoTetrahedron class.
 * Tetrahedron is defined by 4 vertices.
 */
GeoTetrahedron::GeoTetrahedron() {
  nVertices = 4;
  vertices = new Node3D[nVertices];
  cpointsFileName = (std::string)CPOINTS_DIR + "/cpointsTetSurf.msh";
  templateName = "TetrahedronTemplate";
}

/**
 * Destructor of GeoTetrahedron class.
 */
GeoTetrahedron::~GeoTetrahedron() {
  delete[] vertices;
}

/**
 * Initialization of the vertices
 */
void GeoTetrahedron::verticesInit() {
  // at first, we create standard tetrahedron with center at origin of coordinates
  double a = lengths[0]; // the length of tetrahedron edge
  double coord[][4] = { { 0.0, 0.0, sqrt(6.0) * a / 4.0, 1.0 }, \
                        { 0.0, sqrt(3.0) * a / 3.0, -sqrt(6.0) * a / 12.0, 1.0 }, \
                        { a / 2.0, -sqrt(3.0) * a / 6.0, -sqrt(6.0) * a / 12.0, 1.0 }, \
                        { -a / 2.0, -sqrt(3.0) * a / 6.0, -sqrt(6.0) * a / 12.0, 1.0 } };
  double coord_new[4];
  for (int i = 0; i < nVertices; i++) {
    toNewCoord(coord[i], coord_new); // from standard tetrahedron at origin to tetrahedron with rotation and translation
    vertices[i].init(coord_new, i);
  }
}

/**
 * Calculate the volume of the tetrahedron
 */
double GeoTetrahedron::volume() const {
  double x[4], y[4], z[4]; // Cartesian coordinates of the tetrahedron vertices
  for (int i = 0; i < nVertices; i++) {
    x[i] = vertices[i].getX();
    y[i] = vertices[i].getY();
    z[i] = vertices[i].getZ();
  }
  double detD = Det(x[0], x[1], x[2], x[3],
                    y[0], y[1], y[2], y[3],
                    z[0], z[1], z[2], z[3],
                    1.0, 1.0, 1.0, 1.0); // some key determinant
  return fabs(detD) / 6.0; // tetrahedron volume
}

/**
 * Check if tetrahedron contains some point with accuracy in tolerance
 * \param[in] point - point that we want to check
 * \param[in] tolerance - accuracy of checking
 */
bool GeoTetrahedron::hasPoint(Node3D *point, double tolerance) {
  double x[4], y[4], z[4]; // Cartesian coordinates of the tetrahedron vertices
  for (int i = 0; i < nVertices; i++) {
    x[i] = vertices[i].getX();
    y[i] = vertices[i].getY();
    z[i] = vertices[i].getZ();
  }
  double px = point->getX(); // Cartesian coordinates of the point under investigation
  double py = point->getY();
  double pz = point->getZ();

  double V = volume(); // the volume of the main tetrahedron
  double v[4]; // the volumes of the 'small' tetrahedra

  double S = Det(px, py, pz, 1.0,
                 x[1], y[1], z[1], 1.0,
                 x[2], y[2], z[2], 1.0,
                 x[3], y[3], z[3], 1.0);
  v[0] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          px, py, pz, 1.0,
          x[2], y[2], z[2], 1.0,
          x[3], y[3], z[3], 1.0);
  v[1] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          x[1], y[1], z[1], 1.0,
          px, py, pz, 1.0,
          x[3], y[3], z[3], 1.0);
  v[2] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          x[1], y[1], z[1], 1.0,
          x[2], y[2], z[2], 1.0,
          px, py, pz, 1.0);
  v[3] = fabs(S) / 6.0;
  
  if (fabs(V - v[0] - v[1] - v[2] - v[3]) / V < tolerance)
    return true; // tetrahedron contains the point, because the sum of small volumes coicides with the main volume

  return false;
}

/**
 * Create control points of unit tetrahedron: edges lengths are 1, center at origin
 */
void GeoTetrahedron::createUnitElementCPoints(double cl) {
  std::string procName = "GeoTetrahedron::createUnitElementCPoints";
  std::string geoName = "cpointsTetSurf.geo";
  std::ofstream out(geoName.c_str());
  frequire(out, geoName, procName);
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n";
  out << "x0 = 0; y0 = 0; z0 = 0.612372;\n";
  out << "x1 = 0; y1 = 0.57735; z1 = -0.204124;\n";
  out << "x2 = 0.5; y2 = -0.288675; z2 = -0.204124;\n";
  out << "x3 = -0.5; y3 = -0.288675; z3 = -0.204124;\n";
  out << "cl = " << cl << ";\n";
  out << "surfNumber = 1;\n";
  out << "volNumber = 0;\n";
  out << "Call " << templateName << ";\n";
  out << "Physical Surface(3101) = { tetSurfaces[] };\n";
  out.close();
  std::string createMesh = GMSH_BIN + GMSH_BUILD_CPOINTS_OPTIONS + geoName + " -o " + cpointsFileName;
  require(!system(createMesh.c_str()), "Mesh cannot be built!", procName);
}

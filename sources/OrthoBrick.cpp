//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#include "OrthoBrick.h"
#include "Config.h"
#include "Require.h"
#include <fstream>

/**
 * Constructor of OrthoBrick class.
 * OrthoBrick is defined by 8 vertices.
 */
OrthoBrick::OrthoBrick() {
  nVertices = 8;
  vertices = new Node3D[nVertices];
  cpointsFileName = (std::string)CPOINTS_DIR + "/cpointsBrickSurf.msh";
  templateName = "OrthoBrickTemplate";
}

/**
 * Destructor of OrthoBrick class.
 */
OrthoBrick::~OrthoBrick() {
  delete[] vertices;
}

/**
 * Initialization of the vertices
 * lengths[0] - length of side parallel to x-axis
 * lengths[1] - length of side parallel to y-axis
 * lengths[2] - length of side parallel to z-axis
 */
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

/**
 * Calculate the volume of the orthobrick
 */
double OrthoBrick::volume() const {
  return lengths[0] * lengths[1] * lengths[2];
}

/**
 * Check the containing of the point.
 * \param[in] point - point that we want to check
 * \param[in] tolerance - accuracy of checking
 */
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
  const double px = coord_origin[0]; // coordinates of point
  const double py = coord_origin[1];
  const double pz = coord_origin[2];
  const double x0 = -0.5 * lengths[0]; // limits of standard orthobrick
  const double y0 = -0.5 * lengths[1];
  const double z0 = -0.5 * lengths[2];
  const double x1 = -x0;
  const double y1 = -y0;
  const double z1 = -z0;
  
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

/**
 * Create control points of unit orthobrick: sides lengths are 1, center at origin
 */
void OrthoBrick::createUnitElementCPoints(double cl) {
  std::string procName = "OrthoBrick::createUnitElementCPoints";
  std::string geoName = "cpointsBrickSurf.geo";
  std::ofstream out(geoName.c_str());
  frequire(out, geoName, procName);
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n";
  out << "x0 = -0.5; y0 = -0.5; z0 = -0.5;\n";
  out << "x1 = 0.5; y1 = -0.5; z1 = -0.5;\n";
  out << "x2 = -0.5; y2 = -0.5; z2 = 0.5;\n";
  out << "x3 = 0.5; y3 = -0.5; z3 = 0.5;\n";
  out << "x4 = -0.5; y4 = 0.5; z4 = -0.5;\n";
  out << "x5 = 0.5; y5 = 0.5; z5 = -0.5;\n";
  out << "x6 = -0.5; y6 = 0.5; z6 = 0.5;\n";
  out << "x7 = 0.5; y7 = 0.5; z7 = 0.5;\n";
  out << "cl = " << cl << ";\n";
  out << "surfNumber = 1;\n";
  out << "volNumber = 0;\n";
  out << "Call " << templateName << ";\n";
  out << "Physical Surface(3101) = { orthoBrickSurfaces[] };\n";
  out.close();
  std::string createMesh = GMSH_BIN + GMSH_BUILD_CPOINTS_OPTIONS + geoName + " -o " + cpointsFileName;
  require(!system(createMesh.c_str()), "Mesh cannot be built!", procName);
}

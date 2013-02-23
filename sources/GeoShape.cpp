/*
 * mixture - Copyright (c) 2012-2013 Mikhail Artemiev
 *
 * http://code.google.com/p/mixture
 *
 * This library is provided under the terms of MIT license.
 * See the LICENSE file for license information.
 */

#include "GeoShape.h"
#include "Require.h"
#include "TriangularMesh.h"
#include "MeshTriangle3D.h"
#include <iostream>
#include <boost/filesystem.hpp>

GeoShape::~GeoShape() { }

/**
 * Initialization of geometric shape
 * \param[in] cen - the center of shape
 * \param[in] len - lengths of shape
 * \param[in] rvec - rotation vector
 * \param[in] an - angle of rotation
 * \param[in] chalen - characteristic length
 */
void GeoShape::init(const Node3D &cen, double len[], double rvec[], double an, double chalen) {
  center = cen;
  for (int i = 0; i < 3; i++) {
    lengths[i] = len[i];
    rvector[i] = rvec[i];
  }
  angle = an;
  cl = chalen * (lengths[0] + lengths[1] + lengths[2]) / 3.0;
   
  // normalize the rotation vector
  double norm = Norm(rvector, 3); // the norm of the rotation vector
  require(norm > 1e-8, "The norm of rotation vector is too small!", "GeoShape::init");
  for (int i = 0; i < 3; i++)
    rvector[i] /= norm;
    
  // initialization of the transformation matrix
  tMatrixInit();
  
  // initialization of the inverse transformation matrix
  tInvMatInit();
  
//#ifdef DEBUG
//  checkTMatrices();
//#endif

  // initialization of the vertices must be realized in derived class
  verticesInit();
  
  // initialization of the control points
  if (!boost::filesystem::exists(cpointsFileName))
    createUnitElementCPoints(0.1); // if there is no file with control points we must create it (parameter is a characteristic length of mesh)
  cpointsInit(cpointsFileName);
}

/**
 * Initialization of the transformation matrix
 */
void GeoShape::tMatrixInit() {
  /*
      / xx(1-c)+c   xy(1-c)-zs  xz(1-c)+ys  h \
      | yx(1-c)+zs  yy(1-c)+c   yz(1-c)-xs  g |
  T = | xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c   k |
      \ 0           0           0           1 /,
  where (x, y, z) is a normalized rotation vector, ||(x, y, z)|| = 1
  c = cos(angle), s = sin(angle),
  (h, g, k) is a vector of translation = coordinates of the shape center,
            if shape center was at origin of coordinates
  / a' \       / a \
  | b' |       | b |
  | c' | = T * | c |
  \ 1  /       \ 1 /,
  where (a, b, c) is a vector of coordinates of vertex of shape at the origin,
  (a', b', c') is a vector of coordinates of vertex of shape at corresponding center and rotation
  */
  double c = cos(angle * PI / 180.0);
  double s = sin(angle * PI / 180.0);
  
  tMatrix[0][0] = rvector[0] * rvector[0] * (1.0 - c) + c;
  tMatrix[0][1] = rvector[0] * rvector[1] * (1.0 - c) - rvector[2] * s;
  tMatrix[0][2] = rvector[0] * rvector[2] * (1.0 - c) + rvector[1] * s;
  tMatrix[0][3] = center.getX();
  
  tMatrix[1][0] = rvector[1] * rvector[0] * (1.0 - c) + rvector[2] * s;
  tMatrix[1][1] = rvector[1] * rvector[1] * (1.0 - c) + c;
  tMatrix[1][2] = rvector[1] * rvector[2] * (1.0 - c) - rvector[0] * s;
  tMatrix[1][3] = center.getY();
  
  tMatrix[2][0] = rvector[2] * rvector[0] * (1.0 - c) - rvector[1] * s;
  tMatrix[2][1] = rvector[2] * rvector[1] * (1.0 - c) + rvector[0] * s;
  tMatrix[2][2] = rvector[2] * rvector[2] * (1.0 - c) + c;
  tMatrix[2][3] = center.getZ();
  
  tMatrix[3][0] = tMatrix[3][1] = tMatrix[3][2] = 0.0;
  tMatrix[3][3] = 1.0;
}

/**
 * Initialization of the inverse transformation matrix
 */
void GeoShape::tInvMatInit() {
  /*
      /    |   \
      | R  | t |
  T = |    |   |
      |--------|
      \ 0  | 1 /
      
         /    |       \
         | RT | -RT*t |
  T^-1 = |    |       |
         |------------|
         \ 0  |   1   /
  */
  double c = cos(angle * PI / 180.0);
  double s = sin(angle * PI / 180.0);
  
  tInvMat[0][0] = rvector[0] * rvector[0] * (1.0 - c) + c;
  tInvMat[0][1] = rvector[1] * rvector[0] * (1.0 - c) + rvector[2] * s;
  tInvMat[0][2] = rvector[2] * rvector[0] * (1.0 - c) - rvector[1] * s;
  tInvMat[0][3] = -(tInvMat[0][0] * center.getX() + tInvMat[0][1] * center.getY() + tInvMat[0][2] * center.getZ());
  
  tInvMat[1][0] = rvector[0] * rvector[1] * (1.0 - c) - rvector[2] * s;
  tInvMat[1][1] = rvector[1] * rvector[1] * (1.0 - c) + c;
  tInvMat[1][2] = rvector[2] * rvector[1] * (1.0 - c) + rvector[0] * s;
  tInvMat[1][3] = -(tInvMat[1][0] * center.getX() + tInvMat[1][1] * center.getY() + tInvMat[1][2] * center.getZ());
  
  tInvMat[2][0] = rvector[0] * rvector[2] * (1.0 - c) + rvector[1] * s;
  tInvMat[2][1] = rvector[1] * rvector[2] * (1.0 - c) - rvector[0] * s;
  tInvMat[2][2] = rvector[2] * rvector[2] * (1.0 - c) + c;
  tInvMat[2][3] = -(tInvMat[2][0] * center.getX() + tInvMat[2][1] * center.getY() + tInvMat[2][2] * center.getZ());
  
  tInvMat[3][0] = tInvMat[3][1] = tInvMat[3][2] = 0.0;
  tInvMat[3][3] = 1.0;
}

/**
 * Check transformation and its inverse matrices
 */
void GeoShape::checkTMatrices() {
  std::cout << "GeoShape::checkTMatrices\n";
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double res = 0.0;
      for (int k = 0; k < 4; k++) {
        res += tMatrix[i][k] * tInvMat[k][j];
      }
      std::cout << (fabs(res) > 1e-8 ? res : 0) << " ";
    }
    std::cout << "\n";
  }
}

/**
 * Transformation matrix multiplies on vector of coordinates to get new coordinates
 * \param[in] from - array of coordinates that we need to transform
 * \param[out] to - array of resulting coordinates
 */
void GeoShape::toNewCoord(double from[], double to[]) {
  for (int i = 0; i < 4; i++) {
    to[i] = 0.0;
    for (int j = 0; j < 4; j++)
      to[i] += tMatrix[i][j] * from[j];
  }
}

/**
 * Inverse transformation matrix multiplies on vector of coordinates to get "old" coordinates
 * \param[in] from - array of coordinates that we need to transform
 * \param[out] to - array of resulting coordinates
 */
void GeoShape::toOldCoord(double from[], double to[]) {
  for (int i = 0; i < 4; i++) {
    to[i] = 0.0;
    for (int j = 0; j < 4; j++)
      to[i] += tInvMat[i][j] * from[j];
  }
}

/**
 * Get the center of shape
 */
Node3D* GeoShape::getCenter() { return &center; }

/**
 * Get the number of vertices
 */
inline int GeoShape::getnVertices() const { return nVertices; }

/**
 * Get some vertex of shape.
 * \param[in] num - the number of vertex
 */
inline Node3D* GeoShape::getVertex(int num) const {
  require(num >= 0 && num < nVertices, "Incorrect input parameter!", "GeoShape::getVertex");
  return &vertices[num];
}

/**
 * Get some length.
 * \param[in] num - what length do you want
 */
inline double GeoShape::getLen(int num) const {
  require(num >= 0 && num < 3, "Incorrect input parameter!", "GeoShape::getLen");
  return lengths[num];
}

/**
 * Get some component of rotation vector
 * \param[in] num - the number of component
 */
inline double GeoShape::getRVec(int num) const {
  require(num >= 0 && num < 3, "Incorrect input parameter!", "GeoShape::getRVec");
  return rvector[num];
}

/**
 * Get the angle of rotation
 */
inline double GeoShape::getAngle() const { return angle; }

/**
 * Get absolute characteristic value
 */
inline double GeoShape::getCL() const { return cl; }

/**
 * Get the number of control points
 */
inline int GeoShape::getnControlPoints() const { return nControlPoints; }

/**
 * Get control point.
 * \param[in] num - the number of control point
 */
inline Node3D* GeoShape::getControlPoint(int num) const {
  require(num >= 0 && num < nControlPoints, "Incorrect input parameter!", "GeoShape::getControlPoint");
  return &controlPoints[num];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \param[out] min - array of minimal limits
 * \param[out] max - array of maximal limits
 */
void GeoShape::getLimits(bool extendedShape, double min[], double max[]) {
  Node3D *array; // array of nodes (vertices or controlPoints)
  int dim; // dimension of array of nodes
  if (!extendedShape) { // limits of real shape
    array = &vertices[0];
    dim = nVertices;
  }
  else { // limits of extended shape
    array = &controlPoints[0];
    dim = nControlPoints;
  }
  min[0] = max[0] = array[0].getX();
  min[1] = max[1] = array[0].getY();
  min[2] = max[2] = array[0].getZ();
  double x, y, z;
  for (int i = 1; i < dim; i++) {
    x = array[i].getX(); // coordinates of the current vertex
    y = array[i].getY();
    z = array[i].getZ();
    if (x < min[0]) min[0] = x; // define max and min coordinate for every direction
    if (x > max[0]) max[0] = x;
    if (y < min[1]) min[1] = y;
    if (y > max[1]) max[1] = y;
    if (z < min[2]) min[2] = z;
    if (z > max[2]) max[2] = z;
  }
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return x-component of minimal limit
 */
double GeoShape::getMinX(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return min[0];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return x-component of maximal limit
 */
double GeoShape::getMaxX(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return max[0];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return y-component of minimal limit
 */
double GeoShape::getMinY(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return min[1];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return y-component of maximal limit
 */
double GeoShape::getMaxY(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return max[1];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return z-component of minimal limit
 */
double GeoShape::getMinZ(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return min[2];
}

/**
 * After rotation minimal and maximal coordinates of shape change,
 * therefore we need to define new limits of shape
 * \param[in] extendedShape - use additional layer to define new limits (true) or define real limits (false)
 * \return z-component of maximal limit
 */
double GeoShape::getMaxZ(bool extendedShape) {
  double min[3], max[3];
  getLimits(extendedShape, min, max);
  return max[2];
}

/**
 * Check the intersection with some outside bigger box (superelement)
 * \param[in] startPoint - start point of master brick
 * \param[in] sidesLen - lengths of sides of master brick
 */
bool GeoShape::checkIntersection(Node3D &startPoint, double sidesLen[]) {
  // superelement limits
  const double x0 = startPoint.getX();
  const double y0 = startPoint.getY();
  const double z0 = startPoint.getZ();
  const double x1 = x0 + sidesLen[0];
  const double y1 = y0 + sidesLen[1];
  const double z1 = z0 + sidesLen[2];

  // limits of extended shape
  double min[3], max[3];
  getLimits(true, min, max);
  
  // if all limits are inside superelement, then there is no intersection.
  // we can compare in such a manner because superelement's faces are parallel to coordinate axes
  if (max[0] > x1 - CHECK_INTERSECTION_TOLERANCE || \
      min[0] < x0 + CHECK_INTERSECTION_TOLERANCE || \
      max[1] > y1 - CHECK_INTERSECTION_TOLERANCE || \
      min[1] < y0 + CHECK_INTERSECTION_TOLERANCE || \
      max[2] > z1 - CHECK_INTERSECTION_TOLERANCE || \
      min[2] < z0 + CHECK_INTERSECTION_TOLERANCE)
    return true; // there is an intersection

  return false; // there is no intersection
}

/**
 * Check the intersection with another geometric shape
 * \param[in] sh - another geometric shape
 */
bool GeoShape::checkIntersection(GeoShape *sh) {
  // to check this intersection we need to
  // have an array of points (control points)
  // and check every point if it is in our shape

  GeoShape *one = 0, *two = 0;
  // we check every point from the second shape with the first one,
  // so the first shape should be bigger to avoid some situations
  if (this->volume() >= sh->volume()) {
    one = this;
    two = sh;
  } else {
    one = sh;
    two = this;
  }
  
  for (int cp = 0; cp < two->getnControlPoints(); cp++) {
    if (one->hasPoint(two->getControlPoint(cp), CHECK_INTERSECTION_TOLERANCE))
                   // if our shape contain any control point of 'sh' shape,
      return true; // then this is an intersection
  }

  return false; // there is no intersection
}

/**
 * \brief Initialization of the control points.
 * \param[in] masterElementFileName - the name of .msh file where control points are saved.
 * For every shape we should have such file.
 * .msh file is defined for standard shape only (i.e. for sphere: diameter = 1, center at (0, 0, 0))
 */
void GeoShape::cpointsInit(std::string masterElementFileName) {
  // NOTE! file with control points (masterElementFileName) MUST exist for correct calculation.
  // file with control points contains the mesh of unit shape,
  // i.e. lengths[] = { 1, 1, 1 }, and with a center at the origin of coordinates
  
  require(boost::filesystem::exists(masterElementFileName),
          "File with control points '" + masterElementFileName + "' doesn't exist!" + \
          "You can create it using static function 'createUnitElementCPoints'.", "GeoShape::cpointsInit");
  
  TriangularMesh sm; // mesh containing shape triangulation
  sm.readFromGmsh(masterElementFileName); // read mesh from file
  
  int nMeshNodes = sm.getnNodes(); // the number of mesh nodes

  nControlPoints = 2 * nMeshNodes; // the number of control points is in 2 times more than the number of mesh nodes
  controlPoints = new Node3D[nControlPoints];
  
  double coord[4], coord_new[4];
  coord[3] = 1.0; // constant component
  Node3D *node;
  for (int cp = 0; cp < nMeshNodes; cp++) {
    node = sm.getNode(cp); // current node
    // NOTE! control points correspond a little bit bigger shape than it exists in real
    coord[0] = node->getX() * lengths[0] * (1.0 + 0.01 * ADDITIONAL_LAYER_THICKNESS); // scaling
    coord[1] = node->getY() * lengths[1] * (1.0 + 0.01 * ADDITIONAL_LAYER_THICKNESS);
    coord[2] = node->getZ() * lengths[2] * (1.0 + 0.01 * ADDITIONAL_LAYER_THICKNESS);
    toNewCoord(coord, coord_new); // rotation and translation
    controlPoints[cp].init(coord_new, cp); // control point initialization
    // NOTE! control points correspond a little bit smaller shape than it exists in real
    coord[0] = node->getX() * lengths[0] * (1.0 - 0.01 * ADDITIONAL_LAYER_THICKNESS); // scaling
    coord[1] = node->getY() * lengths[1] * (1.0 - 0.01 * ADDITIONAL_LAYER_THICKNESS);
    coord[2] = node->getZ() * lengths[2] * (1.0 - 0.01 * ADDITIONAL_LAYER_THICKNESS);
    toNewCoord(coord, coord_new); // rotation and translation
    controlPoints[cp + nMeshNodes].init(coord_new, cp + nMeshNodes); // control point initialization
  }
  
#ifdef DEBUG
//  printControlPoints(masterElementFileName + "_check.msh", &sm);
#endif
}

/**
 * Print control points to Gmsh file.
 * \param[in] fileName - the name of file for writing
 * \param[in] mesh - the mesh from file with control points of unit shape
 */
void GeoShape::printControlPoints(std::string fileName, TriangularMesh *mesh) {
  std::ofstream out(fileName.c_str());
  frequire(out, fileName, "GeoShape::printControlPoints");
  
  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << nControlPoints << "\n";
  for (int i = 0; i < nControlPoints; i++)
    out << controlPoints[i].getNumber() + 1 << " " \
        << controlPoints[i].getX() << " " \
        << controlPoints[i].getY() << " " \
        << controlPoints[i].getZ() << "\n";
  out << "$EndNodes\n$Elements\n" << mesh->getnTris() << "\n";
  MeshTriangle3D *tri;
  for (int i = 0; i < mesh->getnTris(); i++) {
    tri = mesh->getTri(i); // current triangle
    out << i + 1 << " 2 2 " \
        << tri->getDomain() << " " \
        << tri->getDomain() + 1 << " " \
        << tri->getVertex(0) + 1 << " " \
        << tri->getVertex(1) + 1 << " " \
        << tri->getVertex(2) + 1 << "\n";
  }
  out << "$EndElements\n";
  
  out.close();
}

/**
 * Translate the shape to new place in such a manner that the center would equal to new point.
 * \param[in] cen - new center
 */
void GeoShape::translate(Node3D &cen) {
  // lengths, rotation vector, angle and characteristic length are the same.
    
  // new center
  center = cen;
  
  // new transformation matrix
  tMatrixInit();
  
  // new inverse transformation matrix
  tInvMatInit();
  
#ifdef DEBUG
//  checkTMatrices();
#endif

  // new initialization of the vertices
  verticesInit();
  
  // new initialization of the control points
  cpointsInit(cpointsFileName);
}

/**
 * Print information about geometric shape for Gmsh
 * \param[in] surfNumber - the number of surface
 * \param[in] volNumber - the number of volume
 */
std::string GeoShape::printGeo(int surfNumber, int volNumber) {
  std::string str = "";
  for (int i = 0; i < nVertices; i++)
    str += "x" + d2s(i) + " = " + d2s(vertices[i].getX()) + "; " + \
           "y" + d2s(i) + " = " + d2s(vertices[i].getY()) + "; " + \
           "z" + d2s(i) + " = " + d2s(vertices[i].getZ()) + ";\n";
  str += "cl = " + d2s(cl) + ";\n";
  str += "surfNumber = " + d2s(surfNumber) + ";\n";
  str += "volNumber = " + d2s(volNumber) + ";\n";
  str += "Call " + templateName + ";\n";
  return str;
}


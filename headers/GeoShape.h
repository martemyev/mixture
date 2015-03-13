//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#ifndef GEOSHAPE_H
#define GEOSHAPE_H

#include "Node3D.h"
#include "Mathematics.h"

/**
 * Constant value for checking intersection
 */
#define CHECK_INTERSECTION_TOLERANCE 1e-12

/**
 * Additional layer is needed to guarantee that inclusions won't intersect each other and Gmsh will be able to build mesh
 * This layer is defined in percents of length of shape side (or axis)
 */
#define ADDITIONAL_LAYER_THICKNESS 0.1

class TriangularMesh;

/**
 * Class for geometric shapes
 */
class GeoShape {
protected:
  Node3D center; // center of mass of geometric shape (or another point that can be chosen as center)
  int nVertices; // number of vertices (or some other points that can be chosen as vertices)
  Node3D *vertices; // array of vertices
  double lengths[3]; // length of shape for every direction
  double rvector[3]; // vector of rotation (normalized)
  double angle; // angle of rotation (in degrees)
  double cl; // characteristic length (for mesh creation), absolute value according to average length of shape
  double tMatrix[4][4]; // transformation matrix (rotation + translation)
  double tInvMat[4][4]; // inverse transformation matrix
  std::string cpointsFileName; // name of the file with mesh, that is used to initialize control points.
                               // this name MUST be initialized in constructor of derived class
  int nControlPoints; // the number of control points
  Node3D *controlPoints; // points for checking the intersection
  std::string templateName; // name of template in Gsmh file
  
  void tMatrixInit(); // initialization of transformation matrix
  void tInvMatInit(); // initialization of inverse transformation matrix
  void checkTMatrices(); // check product of tMatrix and tInvMat (must equal to unit matrix)
  
public:
  //GeoShape(const GeoShape&);
  //GeoShape& operator =(const GeoShape&);

  virtual ~GeoShape();

  void init(const Node3D &, double[], double[], double, double); // initialization procedure
  
  void toNewCoord(double[], double[]); // means: transformation matrix multiplies on vector of coordinates to get new coordinates
  void toOldCoord(double[], double[]); // means: inverse transformation matrix multiplies on vector of coordinates to get "old" coordinates
  void translate(Node3D&); // translate the shape to new place in such a manner that the center would equal to new point

  Node3D* getCenter(); // get the center of shape
  int getnVertices() const; // get the number of vertices
  Node3D* getVertex(int) const; // get some vertex of shape
  double getLen(int) const; // get some length
  double getRVec(int) const; // get some component of rotation vector
  double getAngle() const; // get the angle of rotation
  double getCL() const; // get absolute characteristic value
  int getnControlPoints() const; // get the number of control points
  Node3D* getControlPoint(int) const; // get control point
  
  void getLimits(bool, double[], double[]); // after rotation minimal and maximal coordinates of shape change,
                                            // therefore we need a procedure to define new limits of shape
  double getMinX(bool);
  double getMaxX(bool);
  double getMinY(bool);
  double getMaxY(bool);
  double getMinZ(bool);
  double getMaxZ(bool);
  
  virtual void verticesInit() = 0; // initialization of the vertices
  void cpointsInit(std::string); // initialization of the control points
  virtual double volume() const = 0; // calculate volume of the shape
  bool checkIntersection(Node3D&, double[]); // check the intersection with some outside bigger box
  bool checkIntersection(GeoShape*); // check the intersection with another geometric shape
  virtual bool hasPoint(Node3D*, double) = 0; // check if the shape contains some point with accuracy in tolerance
  std::string printGeo(int, int); // print information about geometric shape for Gmsh
  virtual void createUnitElementCPoints(double) = 0; // create control points of unit element for every geometric shape
  
  void printControlPoints(std::string, TriangularMesh*); // print control points to Gmsh file
};
  
#endif

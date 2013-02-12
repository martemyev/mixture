#ifndef CYLINDER_H
#define CYLINDER_H

#include "GeoShape.h"

class Node3D;

/**
 * Cylinder is a geometric shape
 */
class Cylinder : public GeoShape {
public:
  Cylinder();
  ~Cylinder();

  void verticesInit(); // initialization of the vertices
  double volume() const; // calculate volume of the cylinder
  bool hasPoint(Node3D*, double); // check the containing of the point
  void createUnitElementCPoints(double); // create control points of unit cylinder
};

#endif

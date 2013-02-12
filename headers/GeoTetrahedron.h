#ifndef GEOTETRAHEDRON_H
#define GEOTETRAHEDRON_H

#include "GeoShape.h"

class Node3D;

/**
 * Tetrahedron is a geometric shape
 */
class GeoTetrahedron : public GeoShape {
public:
  GeoTetrahedron();
  ~GeoTetrahedron();
  
  void verticesInit(); // initialization of the vertices
  double volume() const; // calculate volume of the tetrahedron
  bool hasPoint(Node3D*, double); // check if tetrahedron contains some point with accuracy in tolerance
  void createUnitElementCPoints(double); // create control points of unit tetrahedron
};

#endif

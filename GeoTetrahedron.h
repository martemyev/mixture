#ifndef GEOTETRAHEDRON_H
#define GEOTETRAHEDRON_H

#include "GeoShape.h"

class Node3D;

// geometric shape - tetrahedron
class GeoTetrahedron : public GeoShape {
protected:
public:
  GeoTetrahedron();
  ~GeoTetrahedron();
  
  void verticesInit(); // initialization of the vertices
  double volume(); // calculate volume of the tetrahedron
  bool hasPoint(Node3D*, double); // check if tetrahedron contains some point with accuracy in tolerance
};

#endif

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "GeoShape.h"

class Node3D;

// geometric shape - ellipsoid
class Ellipsoid : public GeoShape {
protected:
public:
  Ellipsoid();
  ~Ellipsoid();
  
  void verticesInit(); // initialization of the vertices
  double volume(); // calculate volume of the ellipsoid
  bool hasPoint(Node3D*, double); // check if ellipsoid contains some point with accuracy in tolerance
};

#endif

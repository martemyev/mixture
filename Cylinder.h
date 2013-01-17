#ifndef CYLINDER_H
#define CYLINDER_H

#include "GeoShape.h"

class Node3D;

// geometric shape - cylinder
class Cylinder : public GeoShape {
protected:
public:
  Cylinder();
  ~Cylinder();
  
  void verticesInit(); // initialization of the vertices
  double volume(); // calculate volume of the cylinder
  bool hasPoint(Node3D*, double); // check the containing of the point
};

#endif

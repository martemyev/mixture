#ifndef ORTHOBRICK_H
#define ORTHOBRICK_H

#include "GeoShape.h"

class Node3D;

// geometric shape - hexahedron with faces that are orthogonal to each other
class OrthoBrick : public GeoShape {
protected:
public:
  OrthoBrick();
  ~OrthoBrick();
  
  void verticesInit(); // initialization of the vertices
  double volume(); // calculate volume of the brick
  bool hasPoint(Node3D*, double); // check the containing of the point
};

#endif

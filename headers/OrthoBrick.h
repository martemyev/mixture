//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#ifndef ORTHOBRICK_H
#define ORTHOBRICK_H

#include "GeoShape.h"

class Node3D;

/**
 * Geometric shape - hexahedron with faces that are orthogonal to each other
 */
class OrthoBrick : public GeoShape {
public:
  OrthoBrick();
  ~OrthoBrick();
  
  void verticesInit(); // initialization of the vertices
  double volume() const; // calculate volume of the brick
  bool hasPoint(Node3D*, double); // check the containing of the point
  void createUnitElementCPoints(double); // create control points of unit orthobrick
};

#endif

//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "GeoShape.h"

class Node3D;

/**
 * Ellipsoid is a geometric shape.
 * This class is used for spheres also, because sphere is a particular case of ellipsoid (when all axis are equal).
 */
class Ellipsoid : public GeoShape {
public:
  Ellipsoid();
  ~Ellipsoid();
  
  void verticesInit(); // initialization of the vertices
  double volume() const; // calculate volume of the ellipsoid
  bool hasPoint(Node3D*, double); // check if ellipsoid contains some point with accuracy in tolerance
  void createUnitElementCPoints(double); // create control points of unit ellipsoid (i.e. sphere)
};

#endif

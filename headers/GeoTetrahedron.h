/*
 * mixture - Copyright (c) 2012-2013 Mikhail Artemiev
 *
 * http://code.google.com/p/mixture
 *
 * This library is provided under the terms of MIT license.
 * See the LICENSE file for license information.
 */

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

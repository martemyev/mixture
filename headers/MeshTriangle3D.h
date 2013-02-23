/*
 * mixture - Copyright (c) 2012-2013 Mikhail Artemiev
 *
 * http://code.google.com/p/mixture
 *
 * This library is provided under the terms of MIT license.
 * See the LICENSE file for license information.
 */

#ifndef MESHTRIANGLE3D_H
#define MESHTRIANGLE3D_H

#include "Node3D.h"
#include "Require.h"
#include "Convert.h"
#include "Mathematics.h"
#include <iostream>

#define ORDERED_TRI3D_VERTICES // triangle vertices are ordered

/**
 * Triangle for meshes
 */
class MeshTriangle3D {
private:

  static const int nVertices = 3; // the number of vertices
  int vertices[nVertices]; // the vertices of the triangle
  int domain; // physical domain
  
  double x[nVertices], y[nVertices], z[nVertices]; // the coordinates of the vertices
  double area; // the area of the triangle

public:
  MeshTriangle3D()
    : domain(0)
  { }

  MeshTriangle3D(const MeshTriangle3D &t) {
    for (int i = 0; i < nVertices; i++) {
      vertices[i] = t.vertices[i];
      x[i] = t.x[i];
      y[i] = t.y[i];
      z[i] = t.z[i];
    }
    area = t.area;
    domain = t.domain;
  }

  MeshTriangle3D& operator =(const MeshTriangle3D &t) {
    for (int i = 0; i < nVertices; i++) {
      vertices[i] = t.vertices[i];
      x[i] = t.x[i];
      y[i] = t.y[i];
      z[i] = t.z[i];
    }
    area = t.area;
    domain = t.domain;
    return *this;
  }

  void init(int ver[], int dom, Node3D *meshNodes) {
    for (int i = 0; i < nVertices; i++)
      vertices[i] = ver[i];
#ifdef ORDERED_TRI3D_VERTICES
    int tmp;
    for (int i = 0; i < nVertices - 1; i++) {
      for (int j = i + 1; j < nVertices; j++) {
        if (vertices[i] > vertices[j]) { // ascending ordering of vertices
          tmp = vertices[i];
          vertices[i] = vertices[j];
          vertices[j] = tmp;
        }
      }
    }
#endif
    for (int i = 0; i < nVertices; i++) {
      x[i] = meshNodes[vertices[i]].getX();
      y[i] = meshNodes[vertices[i]].getY();
      z[i] = meshNodes[vertices[i]].getZ();
    }
    double d1 = Det(x[0], x[1], x[2], y[0], y[1], y[2], 1.0, 1.0, 1.0);
    double d2 = Det(y[0], y[1], y[2], z[0], z[1], z[2], 1.0, 1.0, 1.0);
    double d3 = Det(z[0], z[1], z[2], x[0], x[1], x[2], 1.0, 1.0, 1.0);
    area = 0.5 * sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    domain = dom;
  }

  static int getnVertices() { return nVertices; }
  int getDomain() const { return domain; }

  int getVertex(int k) const {
    require(k >= 0 && k < nVertices, "Incorrect input parameter!", "MeshTriangle3D::getVertex");
    return vertices[k];
  }

  double getArea() const { return area; }
};

#endif

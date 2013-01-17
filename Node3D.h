#ifndef NODE3D_H
#define NODE3D_H

#include <iostream>

// node in 3D space
class Node3D {
protected:
  double x, y, z; // Cartesian coordinates
  int number; // the number of node

public:
  Node3D() {
    x = y = z = 0.0;
    number = -1;
  }

  Node3D(const Node3D &n) {
    x = n.x;
    y = n.y;
    z = n.z;
    number = n.number;
  }

  Node3D& operator =(const Node3D &n) {
    x = n.x;
    y = n.y;
    z = n.z;
    number = n.number;
    return *this;
  }

  double getX() const { return x; } // get X-coordinate
  double getY() const { return y; } // get Y-coordinate
  double getZ() const { return z; } // get Z-coordinate
  int getNumber() const { return number; } // get the number of node

  void init(double _x, double _y, double _z, int _number = -1) {
    x = _x;
    y = _y;
    z = _z;
    number = _number;
  }

  void init(double coord[], int num = -1) {
    init(coord[0], coord[1], coord[2], num);
  }
};

#endif

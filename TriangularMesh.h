#ifndef TRIANGULARMESH_H
#define TRIANGULARMESH_H

#include "MeshTriangle3D.h"
#include <vector>
#include <string>

class Node3D;

class TriangularMesh {
private:
  int nNodes; // the number of mesh nodes
  int nTris; // the number of boundary triangles

  Node3D *nodes; // the mesh nodes
  std::vector<MeshTriangle3D> tris; // the boundary triangles
  
  int getnElemNodes(int); // get the number of nodes that correspond the type of element

public:
  TriangularMesh();
  ~TriangularMesh();

  int getnNodes(); // only mesh nodes
  int getnTris();

  Node3D* getNode(int);
  MeshTriangle3D* getTri(int);

  void readFromGmsh(std::string);
  void readFromGmshBinary(std::string);
  void readFromGmshASCII(std::string);
};

#endif

//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#ifndef TRIANGULARMESH_H
#define TRIANGULARMESH_H

#include "MeshTriangle3D.h"
#include <vector>
#include <string>

class Node3D;

/**
 * Class for triangular surface mesh built by Gmsh
 */
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

  int getnNodes() const; // get the number of mesh nodes
  int getnTris() const; // get the number of mesh triangles

  Node3D* getNode(int); // get the node
  MeshTriangle3D* getTri(int); // get the triangle

  int readFromGmsh(std::string); // read from .msh file
  int readFromGmshBinary(std::string); // read from binary (version 2) format
  int readFromGmshASCII(std::string); // read from ascii (version 2) format
};

#endif

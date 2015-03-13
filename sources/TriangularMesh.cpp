//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#include "TriangularMesh.h"
#include "Node3D.h"
#include "Require.h"
#include "Convert.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

/**
 * Constructor of TriangularMesh class
 */
TriangularMesh::TriangularMesh()
  : nNodes(-1), nTris(-1)
  { }

/**
 * Destructor of TriangularMesh class
 */
TriangularMesh::~TriangularMesh() {
  delete[] nodes;
  tris.clear();
}

/**
 * Get the number of mesh nodes
 */
int TriangularMesh::getnNodes() const { return nNodes; }

/**
 * Get the number of mesh triangles
 */
int TriangularMesh::getnTris() const { return nTris; }

/**
 * Get the mesh node.
 * \param[in] num - the number of node
 */
Node3D* TriangularMesh::getNode(int num) {
  require(num >= 0 && num < nNodes, "Incorrect input parameter!", "TriangularMesh::getNode");
  return &nodes[num];
}

/**
 * Get the mesh triangle.
 * \param[in] num - the number of triangle
 */
MeshTriangle3D* TriangularMesh::getTri(int num) {
  require(num >= 0 && num < nTris, "Incorrect input parameter!", "TriangularMesh::getTri");
  return &tris[num];
}

/**
 * Read mesh from Gmsh's .msh file.
 * \param[in] fileName - the name of .msh file
 */
int TriangularMesh::readFromGmsh(std::string fileName) {

  clock_t t1, t2;
  t1 = clock();

  std::string procName = "TriangularMesh::readFromGmsh";

  FILE *in = fopen(fileName.c_str(), "r");
  frequire(in, fileName, procName);

  char string[256];
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $MeshFormat
  require(!strncmp(string, "$MeshFormat", 11), "Incorrect mesh format ($MeshFormat)!", procName);

  double version;
  int binary, dsize, scn;
  scn = fscanf(in, "%lf %d %d", &version, &binary, &dsize);
  fclose(in);

  int ret = -1; // unsuccessful reading by default
  if (binary)
    ret = readFromGmshBinary(fileName); // read mesh in binary format
  else
    ret = readFromGmshASCII(fileName); // read mesh in text format

  t2 = clock();

#ifdef DEBUG
  std::cout << procName << ": time (sec) = " << (t2 - t1) / (double)CLOCKS_PER_SEC << std::endl;
#endif

  return ret;
}

/**
 * Read the mesh in Gmsh's ASCII (version 2) format.
 * \param[in] fileName - the name of .msh file
 */
int TriangularMesh::readFromGmshASCII(std::string fileName) {

  std::string procName = "TriangularMesh::readFromGmshASCII";

  std::ifstream in(fileName.c_str());
  frequire(in, fileName, procName);

  std::string tmp;

  getline(in, tmp); // $MeshFormat

  double version;
  int binary, dsize;

  in >> version >> binary >> dsize;
  getline(in, tmp); // the rest of the line

  require(binary == 0, "This mesh file is not in ASCII format!", procName);

  getline(in, tmp); // $EndMeshFormat
  getline(in, tmp); // $Nodes

  in >> nNodes;
  getline(in, tmp); // the rest of the line

  nodes = new Node3D[nNodes];
  int number;
  double x, y, z;

  for (int i = 0; i < nNodes; i++) {
    in >> number >> x >> y >> z;
    number--;
    require(number == i, "The sequence of numbers of nodes is not dense!", procName);
    nodes[i].init(x, y, z, number);
  }
  getline(in, tmp); // the rest of the line
  getline(in, tmp); // $EndNodes
  getline(in, tmp); // $Elements

  int nElements; // all elements
  in >> nElements;
  getline(in, tmp); // the rest of the line

  int elType, nTags;
  int nElemNodes, *nod;
  int physDomain, elemDomain, partition;
  MeshTriangle3D tr;

  int k;
  for (int elem = 0; elem < nElements; elem++) {
    in >> k >> elType >> nTags;
    if (nTags > 0) in >> physDomain; // physical domain
    if (nTags > 1) in >> elemDomain; // elementary domain
    if (nTags > 2) in >> partition; // domain decomposition (Chaco, Metis, etc)
    nElemNodes = getnElemNodes(elType); // how much nodes we need
    nod = new int[nElemNodes];
    for (int i = 0; i < nElemNodes; i++) {
      in >> nod[i]; // read the numbers of nodes
      nod[i]--; // Gmsh numerates the nodes from 1
    }
    switch (elType) {
    case 2: // 3-node 3D triangle
      tr.init(nod, physDomain, nodes);
      tris.push_back(tr);
      break;
    default:
      require(false, "Unknown type of the element " + d2s(elType) + "!", procName);
    }
    delete[] nod;
  } // for all elements

  in.close();

  nTris = (int)tris.size();

  return 0;
}

/**
 * Read the mesh in Gmsh's binary (version 2) format.
 * \param[in] fileName - the name of .msh file
 */
int TriangularMesh::readFromGmshBinary(std::string fileName) {

  std::string procName = "TriangularMesh::readFromGmshBinary";

#ifdef DEBUG
  std::cout << procName << std::endl;
#endif

  FILE *in = fopen(fileName.c_str(), "rb");
  frequire(in, fileName, procName);

  char string[256];
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $MeshFormat
  require(!strncmp(string, "$MeshFormat", 11), "Incorrect file format ($MeshFormat)!", procName);

  double version;
  int binary, dsize, scn;
  scn = fscanf(in, "%lf %d %d", &version, &binary, &dsize);
  require(binary == 1, "This mesh file is not in binary format!", procName);

  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // the rest of the line (the symbol of new line)

  int one;
  scn = fread(&one, sizeof(int), 1, in);
  require(one == 1, "Incorrect file format (binary one is wrong)!", procName);

  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // the rest of the line
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $EndMeshFormat
  require(!strncmp(string, "$EndMeshFormat", 14), "Incorrect file format ($EndMeshFormat)!", procName);

  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $Nodes
  require(!strncmp(string, "$Nodes", 6), "Incorrect file format ($Nodes)!", procName);

  scn = fscanf(in, "%d", &nNodes);
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // the rest of the line

  nodes = new Node3D[nNodes];

  int number;
  double coordinates[3];
  for (int i = 0; i < nNodes; i++) {
    scn = fread(&number, sizeof(int), 1, in); // the number of node
    number--; // Gmsh numerates the nodes from 1
    require(number == i, "The sequence of numbers of nodes is not dense!", procName);
    scn = fread(coordinates, sizeof(double), 3, in); // Cartesian coordinates
    nodes[i].init(coordinates, number);
  }

  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // the rest of the line
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $EndNodes
  require(!strncmp(string, "$EndNodes", 9), "Incorrect file format ($EndNodes)!", procName);
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // $Elements
  require(!strncmp(string, "$Elements", 9), "Incorrect file format ($Elements)!", procName);

  int nElements;
  scn = fscanf(in, "%d", &nElements);
  require(fgets(string, 256, in) != NULL, "fgets returned NULL", procName); // the rest of the line

  int nElemPart = 0;
  int header[3], elType, amount, nTags;
  int nElemNodes, nData, *data;
  int physDomain, elemDomain, partition;
  int *nod;

  while (nElemPart < nElements) {
    scn = fread(header, sizeof(int), 3, in);
    elType = header[0]; // the type of the element
    amount = header[1]; // the number of the elements of this type
    nTags = header[2]; // the number of tags
    
    nElemPart += amount;
    nElemNodes = getnElemNodes(elType); // how much nodes we need
    nData = 1 + nTags + nElemNodes;
    data = new int[nData];
    nod = new int[nElemNodes];
    MeshTriangle3D tr;

    for (int oneElem = 0; oneElem < amount; oneElem++) {
      scn = fread(data, sizeof(int), nData, in);
      number = data[0];
      physDomain = (nTags > 0) ? data[1] : 0; // physical domain - the most important value
      elemDomain = (nTags > 1) ? data[2] : 0; // elementary domain
      partition = (nTags > 2) ? data[3] : 0; // partition (Metis, Chaco, etc)
      for (int i = 0; i < nElemNodes; i++)
        nod[i] = data[nTags + 1 + i] - 1; // Gmsh numerates the nodes from 1

      switch (elType) {
      case 2: // 3-node 3D triangle
        tr.init(nod, physDomain, nodes);
        tris.push_back(tr);
        break;
      default:
        require(false, "Unknown type of the element " + d2s(elType) + "!", procName);
      }
    }
    delete[] data;
    delete[] nod;
  }

  require(nElemPart == nElements, "We read more(less) than nElements elements!", procName);

  fclose(in);

  nTris = (int)tris.size();

  return 0;
}

// get the number of nodes that correspond the type of element
int TriangularMesh::getnElemNodes(int type) {
  int n = -1;
  switch (type) {
  case 1: // 2-node line
    n = 2;
    break;
  case 2: // 3-node triangle
    n = 3;
    break;
  case 4: // 4-node tetrahedron
    n = 4;
    break;
  default:
    require(false, "Unknown type of the element " + d2s(type) + "!", "TriangularMesh::getnElemNodes");
  }
  return n;
}

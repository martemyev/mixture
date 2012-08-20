#include "TriangularMesh.h"
#include "Node3D.h"
#include "Require.h"
#include "Convert.h"
#include "Param.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

TriangularMesh::TriangularMesh() {
  nNodes = nTris = -1;
}

TriangularMesh::~TriangularMesh() {
  delete[] nodes;
  tris.clear();
}

int TriangularMesh::getnNodes() { return nNodes; }

int TriangularMesh::getnTris() { return nTris; }

Node3D* TriangularMesh::getNode(int num) {
  require(num >= 0 && num < nNodes, "Incorrect input parameter!", "TriangularMesh::getNode");
  return &nodes[num];
}

MeshTriangle3D* TriangularMesh::getTri(int num) {
  require(num >= 0 && num < nTris, "Incorrect input parameter!", "TriangularMesh::getTri");
  return &tris[num];
}

void TriangularMesh::readFromGmsh(std::string fileName) {

  clock_t t1, t2;
  t1 = clock();

  std::string procName = "TriangularMesh::readFromGmsh";

  FILE *in = fopen(fileName.c_str(), "r");
  frequire(in, fileName, procName);

  char string[256];
  fgets(string, 256, in); // $MeshFormat
  require(!strncmp(string, "$MeshFormat", 11), "Incorrect file format ($MeshFormat)!", procName);

  double version;
  int binary, dsize;
  fscanf(in, "%lf %d %d", &version, &binary, &dsize);
  fclose(in);
  if (binary)
    readFromGmshBinary(fileName); // read mesh in binary format
  else
    readFromGmshASCII(fileName); // read mesh in text format

  t2 = clock();

#ifdef DEBUG
  std::cout << procName << ": time (sec) = " << (t2 - t1) / (double)CLOCKS_PER_SEC << std::endl;
#endif
}

void TriangularMesh::readFromGmshASCII(std::string fileName) {

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
}

void TriangularMesh::readFromGmshBinary(std::string fileName) {

  std::string procName = "TriangularMesh::readFromGmshBinary";

#ifdef DEBUG
  std::cout << procName << std::endl;
#endif

  FILE *in = fopen(fileName.c_str(), "rb");
  frequire(in, fileName, procName);

  char string[256];
  fgets(string, 256, in); // $MeshFormat
  require(!strncmp(string, "$MeshFormat", 11), "Incorrect file format ($MeshFormat)!", procName);

  double version;
  int binary, dsize;
  fscanf(in, "%lf %d %d", &version, &binary, &dsize);
  require(binary == 1, "This mesh file is not in binary format!", procName);

  fgets(string, 256, in); // the rest of the line (the symbol of new line)

  int one;
  fread(&one, sizeof(int), 1, in);
  require(one == 1, "Incorrect file format (binary one is wrong)!", procName);

  fgets(string, 256, in); // the rest of the line
  fgets(string, 256, in); // $EndMeshFormat
  require(!strncmp(string, "$EndMeshFormat", 14), "Incorrect file format ($EndMeshFormat)!", procName);

  fgets(string, 256, in); // $Nodes
  require(!strncmp(string, "$Nodes", 6), "Incorrect file format ($Nodes)!", procName);

  fscanf(in, "%d", &nNodes);
  fgets(string, 256, in); // the rest of the line

  nodes = new Node3D[nNodes];

  int number;
  double coordinates[3];
  for (int i = 0; i < nNodes; i++) {
    fread(&number, sizeof(int), 1, in); // the number of node
    number--; // Gmsh numerates the nodes from 1
    require(number == i, "The sequence of numbers of nodes is not dense!", procName);
    fread(coordinates, sizeof(double), 3, in); // Cartesian coordinates
    nodes[i].init(coordinates, number);
  }

  fgets(string, 256, in); // the rest of the line
  fgets(string, 256, in); // $EndNodes
  require(!strncmp(string, "$EndNodes", 9), "Incorrect file format ($EndNodes)!", procName);
  fgets(string, 256, in); // $Elements
  require(!strncmp(string, "$Elements", 9), "Incorrect file format ($Elements)!", procName);

  int nElements;
  fscanf(in, "%d", &nElements);
  fgets(string, 256, in); // the rest of the line

  int nElemPart = 0;
  int header[3], elType, amount, nTags;
  int nElemNodes, nData, *data;
  int physDomain, elemDomain, partition;
  int *nod;

  while (nElemPart < nElements) {
    fread(header, sizeof(int), 3, in);
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
      fread(data, sizeof(int), nData, in);
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
    std::string tmp = "Unknown type of the element " + d2s(type) + "!";
    require(false, tmp, "TriangularMesh::getnElemNodes");
  }
  return n;
}

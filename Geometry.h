#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "GeoShape.h"
#include <string>
#include <vector>
#include "pugixml.hpp"

#define MARGIN 0.5 // margin in percents relatively length of superelement side
#define N_TRIES_TO_CREATE_ELEMENT 10 // we try several times to create an element, because new element could intersect the previous one

class Geometry : public pugi::xml_tree_walker { // xml_tree_walker need to walk through parameter file
private:
  int nElements; // the number of all elements (inclusions)
  std::vector<GeoShape*> elements; // all inclusions
  
  int nOrthoBricks; // the number of orthobricks that you want to have
  int nRealOrthoBricks; // the number of orthobricks that you really have
  int orthobricksInit(pugi::xml_node&); // orthobricks initialization
  
  int nSpheres; // the number of spheres that you want to have
  int nRealSpheres; // the number of spheres that you really have
  int spheresInit(pugi::xml_node&); // spheres initialization
  
  int nEllipsoids; // the number of ellipsoids that you want to have
  int nRealEllipsoids; // the number of ellipsoids that you really have
  int ellipsoidsInit(pugi::xml_node&); // ellipsoids initialization

  int nCylinders; // the number of cylinders that you want to have
  int nRealCylinders; // the number of cylinders that you really have
  int cylindersInit(pugi::xml_node&); // cylinders initialization

  int nTetrahedra; // the number of tetrahedra that you want to have
  int nRealTetrahedra; // the number of tetrahedra that you really have
  int tetrahedraInit(pugi::xml_node&); // tetrahedra initialization
  
  Node3D masterBrickStartPoint; // start point of the master brick (superelement)
  double masterBrickLengths[3]; // lengths of sides of the master brick
  double clMasterBrick; // characteristic length for mesh inside master brick
  int masterbrickInit(pugi::xml_node&); // initialization of the master brick

  void writeResults(std::string); // write geometry to the corresponding file
  void writeEmptyMasterBrick(std::string); // write geometry of empty superelement to the corresponding file
  
  // to walk on xml document with geometry parameters
  virtual bool for_each(pugi::xml_node&);
  
  // get some information from xml file
  void getInfoFromParamFile(std::string, std::string, std::string, std::string, \
                            std::string, pugi::xml_node&, int, double, double*);
  void getNumber(pugi::xml_node&, int*, std::string);
  int getValueType(pugi::xml_node&, std::string);
  void getLength(pugi::xml_node&, std::string, double[], int, double, std::string);
  int getRVec(pugi::xml_node&, double[], std::string);
  void getAngle(pugi::xml_node&, double[], std::string);
  double getCL(pugi::xml_node&, std::string);
  void getCenter(pugi::xml_node&, int, int, Node3D*, Node3D, double[], double[], std::string);
  
public:
  Geometry();
  ~Geometry();
  
  int init(double, double, double, double, double, double, std::string, double, std::string); // initialization
  int init(std::string, std::string); // initialization
  static void printTemplates(std::string); // print Gmsh-templates to the file
  
  int getnOrthoBricks();
  int getnRealOrthoBricks();
  int getnCylinders();
  int getnRealCylinders();
  int getnEllipsoids();
  int getnRealEllipsoids();
  int getnSpheres();
  int getnRealSpheres();
  int getnTetrahedra();
  int getnRealTetrahedra();

};
#endif

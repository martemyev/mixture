/*
 * mixture - Copyright (c) 2012-2013 Mikhail Artemiev
 *
 * http://code.google.com/p/mixture
 *
 * This library is provided under the terms of MIT license.
 * See the LICENSE file for license information.
 */

#include "Geometry.h"
#include "OrthoBrick.h"
#include "Ellipsoid.h"
#include "Cylinder.h"
#include "GeoTetrahedron.h"
#include "Config.h"
#include "Require.h"
#include "Convert.h"
#include "Mathematics.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <ctime>
#include <boost/filesystem.hpp>

/**
 * Constructor of Geometry class.
 * All n(umber)* parameters are initialized by 0.
 * Lengths of master brick are initialized by -1.
 */
Geometry::Geometry() {
  nElements = 0; // the number of all inclusions
  nOrthoBricks = nRealOrthoBricks = 0;
  nCylinders = nRealCylinders = 0;
  nEllipsoids = nRealEllipsoids = 0;
  nSpheres = nRealSpheres = 0;
  nTetrahedra = nRealTetrahedra = 0;
  masterBrickLengths[0] = masterBrickLengths[1] = masterBrickLengths[2] = -1;
}

/**
 * Destructor of Geometry class.
 * The list of elements (inclusions) is cleaned.
 */
Geometry::~Geometry() {
  for (int i = 0; i < (int)elements.size(); i++)
    delete elements[i]; // free the memory
  elements.clear(); // clear the list of inclusions
}

/**
 * Initialization of geometry when parameters of master brick is defined directly.
 * \param[in] x0, y0, z0 - start point of brick
 * \param[in] lenX, lenY, lenZ - lengths of sides
 * \param[in] clMain - part of characteristic length of mesh in master brick
 *            clMasterBrick = clMain * (lenX + lenY + lenZ) / 3.0
 * \param[in] geoParamFileName - the name of xml file containing the parameters of geometry
 * \param[in] geoFileName - the name of .geo file that will have the resulting geometry in Gmsh's compatible format
 */
int Geometry::init(double x0, double y0, double z0, \
                   double lenX, double lenY, double lenZ, double clMain, \
                   const std::string &geoParamFileName, const std::string &geofileName) {

  masterBrickStartPoint.init(x0, y0, z0); // start point of the master brick
  masterBrickLengths[0] = lenX; // lengths of sides of the master brick
  masterBrickLengths[1] = lenY;
  masterBrickLengths[2] = lenZ;

  clMasterBrick = clMain * (lenX + lenY + lenZ) / 3.0; // characteristic length

  return init(geoParamFileName, geofileName);
}

/**
 * Initialization of geometry based on xml file only.
 * \param[in] geoParamFileName - the name of xml file containing the parameters of geometry
 * \param[in] geoFileName - the name of .geo file that will have the resulting geometry in Gmsh's compatible format
 */
int Geometry::init(const std::string &geoParamFileName, const std::string &geofileName) {
  std::string procName = "Geometry::init";

  // check that directory with control points files exists.
  // if not we should create it
  if (!boost::filesystem::exists(CPOINTS_DIR))
    boost::filesystem::create_directory(CPOINTS_DIR);

  // check that the templates file exists.
  // if not we should create it
  if (!boost::filesystem::exists(TEMPLATES_DIR + TEMPLATES_FILENAME))
    printTemplates(TEMPLATES_DIR);

#ifdef DEBUG
  std::cout << procName << " for parameters file " << geoParamFileName << " and resulting file " << geofileName << std::endl;
  std::cout << "Superbrick = [" << masterBrickLengths[0] << ", " << masterBrickLengths[1] << ", " << masterBrickLengths[2] << "]" << std::endl;
#endif

//  srand() procedure Main class should call (Main - means that from which this procedure was called)
  srand(time(0) + (int)masterBrickLengths[0] + (int)masterBrickLengths[1] + (int)masterBrickLengths[2]);

  require(boost::filesystem::exists(geoParamFileName), "File " + geoParamFileName + " doesn't exist!", procName);

  // parse an xml file with parameters
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(geoParamFileName.c_str());

  require(result, "File " + geoParamFileName + " is wrong! Please, check it.", procName);

#ifdef DEBUG
  std::cout << "Geometry::init : parse result = " << result << std::endl;
#endif

  doc.traverse(*this); // walk through the file

  nElements = elements.size(); // the number of created inclusions

//  if (nElements == 0)
//    return 1; // error
  
  if (nElements == 0) // if there is no inclusions in superelement, we create empty element (homogenized)
    writeEmptyMasterBrick(geofileName);
  else
    // write resulting geometry
    writeResults(geofileName);

  return 0; // all is good
}

/**
 * Initialization of the master brick.
 * /param[in] node - node of xml tree containing parameters of master brick
 */
int Geometry::masterbrickInit(pugi::xml_node &node) {
  std::string procName = "Geometry::masterbrickInit";

  masterBrickStartPoint.init(0, 0, 0);

  pugi::xml_attribute attr_x = node.attribute("lenX"); // get attributes from xml file
  pugi::xml_attribute attr_y = node.attribute("lenY");
  pugi::xml_attribute attr_z = node.attribute("lenZ");
  require(attr_x && attr_y && attr_z, "Unknown lengths of super brick!", procName);

  masterBrickLengths[0] = attr_x.as_double();
  masterBrickLengths[1] = attr_y.as_double();
  masterBrickLengths[2] = attr_z.as_double();
  require(masterBrickLengths[0] > 0 && masterBrickLengths[1] > 0 && masterBrickLengths[2] > 0, \
          "Lengths of master brick are not positive!", procName);

  clMasterBrick = getCL(node, procName) * (masterBrickLengths[0] + masterBrickLengths[1] + masterBrickLengths[2]) / 3.0;
  return 0;
}

/**
 * Tetrahedra initialization.
 * \param[in] node - node of xml tree containing parameters of tetrahedra
 */
int Geometry::tetrahedraInit(pugi::xml_node &node) {

  std::string procName = "Geometry::tetrahedraInit";

  getNumber(node, &nTetrahedra, procName); // get the number of tetrahedra
  int valueType = getValueType(node, procName); // get the type of values of all
                                                // parameters: absolute (0) or relative (1) (exept cl, that is always absolute)

  // tetrahedron parameters
  double edgeLen[2];
  getLength(node, "edgelen", edgeLen, valueType, masterBrickLengths[0], procName);
  
  double rvector[3]; // rotation vector
  int rvec = getRVec(node, rvector, procName); // -2 - vector is not defined,
                                               // -1 - vector is defined exactly,
                                               // 0 - random, 1/2/3 - along x/y/z-direction
  double angle[] = { 0, 0 }; // angle of rotation
  if (rvec != -2) // only if rotation vector is defined
    getAngle(node, angle, procName); // get angle (min and max values, or exact one)

  double cl = getCL(node, procName); // characteristic length for mesh building

#ifdef DEBUG
  std::cout << "nTetrahedra = " << nTetrahedra << std::endl;
  std::cout << "valueType = " << valueType << std::endl;
  std::cout << "edgeLen = " << edgeLen[0] << ", " << edgeLen[1] << std::endl;
  std::cout << "rvec = " << rvec << std::endl;
  std::cout << "angle = " << angle[0] << ", " << angle[1] << std::endl;
  std::cout << "cl = " << cl << std::endl;
#endif

  // tetrahedra initialization
  int nTries = 0; // the number of unsuccessful tries
  for (int el = 0; el < nTetrahedra && nTries < N_TRIES_TO_CREATE_ELEMENT;) {
    // center of tetrahedron is
    // at center of superelement
    double xcen = masterBrickStartPoint.getX() + 0.5 * masterBrickLengths[0];
    double ycen = masterBrickStartPoint.getY() + 0.5 * masterBrickLengths[1];
    double zcen = masterBrickStartPoint.getZ() + 0.5 * masterBrickLengths[2];
    // length of the edge.
    // if first element is negative, then the second element is exact value
    double elen = (edgeLen[0] > 0 ? randomBetween(edgeLen[0], edgeLen[1]) : edgeLen[1]);
    // rotation vector
    double rv[] = { 0.0, 0.0, 0.0 };
    switch (rvec) {
    case -2: // vector is not defined
      rv[0] = 1.0; // some unit vector (it doesn't change anything because angle is zero in this case)
      break;
    case -1: // vector is defined exactly (by components)
      rv[0] = rvector[0];
      rv[1] = rvector[1];
      rv[2] = rvector[2];
      break;
    case 0: // random vector
      rv[0] = randomBetween(0.0, 1.0);
      rv[1] = randomBetween(0.0, 1.0);
      rv[2] = randomBetween(0.0, 1.0);
      break;
    case 1: // along x-direction
      rv[0] = xcen;
      break;
    case 2: // along y-direction
      rv[1] = ycen;
      break;
    case 3: // along z-direction
      rv[2] = zcen;
      break;
    default:
      require(false, "Unknown type of rotation vector!", procName);
    }
    // angle
    double an = 0.0; // there is no rotation by default
    if (rvec != -2) // if rotation vector is defined
      an = (angle[0] >= 0 ? randomBetween(angle[0], angle[1]) : angle[1]);

#ifdef DEBUG
    std::cout << "rv = " << rv[0] << ", " << rv[1] << ", " << rv[2] << std::endl;
    std::cout << "an = " << an << std::endl;
#endif

    // initialization
    Node3D cen;
    cen.init(xcen, ycen, zcen); // center of tetrahedron
    double len[] = { elen, elen, elen }; // lengths
    GeoTetrahedron *tetrahedron = new GeoTetrahedron;
    tetrahedron->init(cen, len, rv, an, cl); // tetrahedron initialization

    // translation of the tetrahedron:
    // we multiply MARGIN on 0.01 because MARGIN is written in percents regarding to length of superelement side.
    Node3D newCenter;
    double minLim[3], maxLim[3];
    tetrahedron->getLimits(false, minLim, maxLim); // get real limits
    double limits[] = { maxLim[0] - minLim[0],
                        maxLim[1] - minLim[1],
                        maxLim[2] - minLim[2] }; // lengths of real limits

    getCenter(node, nTetrahedra, valueType, &newCenter, \
              masterBrickStartPoint, masterBrickLengths, limits, procName);
    tetrahedron->translate(newCenter); // translate tetrahedron to new center
    
    // check if this tetrahedron intersects superelement faces or other geometric shapes
    bool intersection = false; // there is no intersection by default
    // intersection with superelement faces
    intersection = tetrahedron->checkIntersection(masterBrickStartPoint, masterBrickLengths);
    // intersection with other geometric shapes
    for (int i = 0; i < (int)elements.size() && !intersection; i++)
      intersection = tetrahedron->checkIntersection(elements.at(i));
    if (!intersection) { // if there is no intersection
      elements.push_back(tetrahedron); // add to list
      nTries = 0; // reset counter of bad tries
      el++; // the number of tetrahedra increases
      nRealTetrahedra++; // real number of tetrahedra
    }
    else
      nTries++; // regular unsuccessful attempt
  }

  return 0; // all is good
}

/**
 * Cylinders initialization.
 * \param[in] node - node of xml tree containing parameters of cylinders
 */
int Geometry::cylindersInit(pugi::xml_node &node) {

  std::string procName = "Geometry::cylindersInit";

  getNumber(node, &nCylinders, procName); // get the number of cylinders
  int valueType = getValueType(node, procName); // get the type of values of all
                                                // parameters: absolute (0) or relative (1) (exept cl, that is always absolute)

  // cylinder parameters
  double dMajor[2], dMinor[2], height[2];
  getLength(node, "dMajor", dMajor, valueType, masterBrickLengths[0], procName);
  getLength(node, "dMinor", dMinor, valueType, masterBrickLengths[1], procName);
  getLength(node, "height", height, valueType, masterBrickLengths[2], procName);
  
  double rvector[3]; // rotation vector
  int rvec = getRVec(node, rvector, procName); // -2 - vector is not defined,
                                               // -1 - vector is defined exactly,
                                               // 0 - random, 1/2/3 - along x/y/z-direction
  double angle[] = { 0, 0 }; // angle of rotation
  if (rvec != -2) // only if rotation vector is defined
    getAngle(node, angle, procName); // get angle (min and max values, or exact one)

  double cl = getCL(node, procName); // characteristic length for mesh building

#ifdef DEBUG
  std::cout << "nCylinders = " << nCylinders << std::endl;
  std::cout << "valueType = " << valueType << std::endl;
  std::cout << "dMajor = " << dMajor[0] << ", " << dMajor[1] << std::endl;
  std::cout << "dMinor = " << dMinor[0] << ", " << dMinor[1] << std::endl;
  std::cout << "height = " << height[0] << ", " << height[1] << std::endl;
  std::cout << "rvec = " << rvec << std::endl;
  std::cout << "angle = " << angle[0] << ", " << angle[1] << std::endl;
  std::cout << "cl = " << cl << std::endl;
#endif

  // cylinders initialization
  int nTries = 0; // the number of unsuccessful tries
  for (int cy = 0; cy < nCylinders && nTries < N_TRIES_TO_CREATE_ELEMENT;) {
    // center of cylinder is
    // at center of superelement
    double xcen = masterBrickStartPoint.getX() + 0.5 * masterBrickLengths[0];
    double ycen = masterBrickStartPoint.getY() + 0.5 * masterBrickLengths[1];
    double zcen = masterBrickStartPoint.getZ() + 0.5 * masterBrickLengths[2];
    // lengths.
    // if first element is negative, then the second element is exact value
    double dMaj = (dMajor[0] > 0 ? randomBetween(dMajor[0], dMajor[1]) : dMajor[1]);
    double dMin = (dMinor[0] > 0 ? randomBetween(dMinor[0], dMinor[1]) : dMinor[1]);
    double h = (height[0] > 0 ? randomBetween(height[0], height[1]) : height[1]);
    // rotation vector
    double rv[] = { 0.0, 0.0, 0.0 };
    switch (rvec) {
    case -2: // vector is not defined
      rv[0] = 1.0; // some unit vector (it doesn't change anything because angle is zero in this case)
      break;
    case -1: // vector is defined exactly (by components)
      rv[0] = rvector[0];
      rv[1] = rvector[1];
      rv[2] = rvector[2];
      break;
    case 0: // random vector
      rv[0] = randomBetween(0.0, 1.0);
      rv[1] = randomBetween(0.0, 1.0);
      rv[2] = randomBetween(0.0, 1.0);
      break;
    case 1: // along x-direction
      rv[0] = xcen;
      break;
    case 2: // along y-direction
      rv[1] = ycen;
      break;
    case 3: // along z-direction
      rv[2] = zcen;
      break;
    default:
      require(false, "Unknown type of rotation vector!", procName);
    }
    // angle
    double an = 0.0; // there is no rotation by default
    if (rvec != -2) // if rotation vector is defined
      an = (angle[0] >= 0 ? randomBetween(angle[0], angle[1]) : angle[1]);

#ifdef DEBUG
    std::cout << "rv = " << rv[0] << ", " << rv[1] << ", " << rv[2] << std::endl;
    std::cout << "an = " << an << std::endl;
#endif

    // initialization
    Node3D cen;
    cen.init(xcen, ycen, zcen); // center of cylinder
    double len[] = { dMaj, dMin, h }; // lengths
    Cylinder *cylinder = new Cylinder;
    cylinder->init(cen, len, rv, an, cl); // cylinder initialization

    // translation of the cylinder:
    // we multiply MARGIN on 0.01 because MARGIN is written in percents regarding to length of superelement side.
    Node3D newCenter;
    double minLim[3], maxLim[3];
    cylinder->getLimits(false, minLim, maxLim); // get real limits
    double limits[] = { maxLim[0] - minLim[0],
                        maxLim[1] - minLim[1],
                        maxLim[2] - minLim[2] }; // lengths of real limits

    getCenter(node, nCylinders, valueType, &newCenter, \
              masterBrickStartPoint, masterBrickLengths, limits, procName);
    cylinder->translate(newCenter); // translate cylinder to new center
    
    // check if this cylinder intersects superelement faces or other geometric shapes
    bool intersection = false; // there is no intersection by default
    // intersection with superelement faces
    intersection = cylinder->checkIntersection(masterBrickStartPoint, masterBrickLengths);
    // intersection with other geometric shapes
    for (int i = 0; i < (int)elements.size() && !intersection; i++)
      intersection = cylinder->checkIntersection(elements.at(i));
    if (!intersection) { // if there is no intersection
      elements.push_back(cylinder); // add to list
      nTries = 0; // reset counter of bad tries
      cy++; // the number of cylinders increases
      nRealCylinders++; // real number of cylinders
    }
    else
      nTries++; // regular unsuccessful attempt
  }

  return 0; // all is good
}

/**
 * Ellipsoids initialization.
 * \param[in] node - node of xml tree containing parameters of ellipsoids
 */
int Geometry::ellipsoidsInit(pugi::xml_node &node) {

  std::string procName = "Geometry::ellipsoidsInit";

  getNumber(node, &nEllipsoids, procName); // get the number of ellipsoids
  int valueType = getValueType(node, procName); // get the type of values of all
                                                // parameters: absolute (0) or relative (1) (except cl, that is always absolute)

  // ellipsoid parameters
  double lenX[2], lenY[2], lenZ[2];
  getLength(node, "lenX", lenX, valueType, masterBrickLengths[0], procName);
  getLength(node, "lenY", lenY, valueType, masterBrickLengths[1], procName);
  getLength(node, "lenZ", lenZ, valueType, masterBrickLengths[2], procName);
  
  double rvector[3]; // rotation vector
  int rvec = getRVec(node, rvector, procName); // -2 - vector is not defined,
                                               // -1 - vector is defined exactly,
                                               // 0 - random, 1/2/3 - along x/y/z-direction
  double angle[] = { 0, 0 }; // angle of rotation
  if (rvec != -2) // only if rotation vector is defined
    getAngle(node, angle, procName); // get angle (min and max values, or exact one)

  double cl = getCL(node, procName); // characteristic length for mesh building

#ifdef DEBUG
  std::cout << "nEllipsoids = " << nEllipsoids << std::endl;
  std::cout << "valueType = " << valueType << std::endl;
  std::cout << "lenX = " << lenX[0] << ", " << lenX[1] << std::endl;
  std::cout << "lenY = " << lenY[0] << ", " << lenY[1] << std::endl;
  std::cout << "lenZ = " << lenZ[0] << ", " << lenZ[1] << std::endl;
  std::cout << "rvec = " << rvec << std::endl;
  std::cout << "angle = " << angle[0] << ", " << angle[1] << std::endl;
  std::cout << "cl = " << cl << std::endl;
#endif

  // ellipsoids initialization
  int nTries = 0; // the number of unsuccessful tries
  for (int el = 0; el < nEllipsoids && nTries < N_TRIES_TO_CREATE_ELEMENT;) {
    // center of ellipsoid is
    // at center of superelement
    double xcen = masterBrickStartPoint.getX() + 0.5 * masterBrickLengths[0];
    double ycen = masterBrickStartPoint.getY() + 0.5 * masterBrickLengths[1];
    double zcen = masterBrickStartPoint.getZ() + 0.5 * masterBrickLengths[2];
    // lengths.
    // if first element is negative, then the second element is exact value
    double lx = (lenX[0] > 0 ? randomBetween(lenX[0], lenX[1]) : lenX[1]);
    double ly = (lenY[0] > 0 ? randomBetween(lenY[0], lenY[1]) : lenY[1]);
    double lz = (lenZ[0] > 0 ? randomBetween(lenZ[0], lenZ[1]) : lenZ[1]);
    // rotation vector
    double rv[] = { 0.0, 0.0, 0.0 };
    switch (rvec) {
    case -2: // vector is not defined
      rv[0] = 1.0; // some unit vector (it doesn't change anything because angle is zero in this case)
      break;
    case -1: // vector is defined exactly (by components)
      rv[0] = rvector[0];
      rv[1] = rvector[1];
      rv[2] = rvector[2];
      break;
    case 0: // random vector
      rv[0] = randomBetween(0.0, 1.0);
      rv[1] = randomBetween(0.0, 1.0);
      rv[2] = randomBetween(0.0, 1.0);
      break;
    case 1: // along x-direction
      rv[0] = xcen;
      break;
    case 2: // along y-direction
      rv[1] = ycen;
      break;
    case 3: // along z-direction
      rv[2] = zcen;
      break;
    default:
      require(false, "Unknown type of rotation vector!", procName);
    }
    // angle
    double an = 0.0; // there is no rotation by default
    if (rvec != -2) // if rotation vector is defined
      an = (angle[0] >= 0 ? randomBetween(angle[0], angle[1]) : angle[1]);

    // initialization
    Node3D cen;
    cen.init(xcen, ycen, zcen); // center of ellipsoid
    double len[] = { lx, ly, lz }; // lengths
    Ellipsoid *ellipsoid = new Ellipsoid;
    ellipsoid->init(cen, len, rv, an, cl); // ellipsoid initialization

    // translation of the ellipsoid:
    // we multiply MARGIN on 0.01 because MARGIN is written in percents regarding to length of superelement side.
    Node3D newCenter;
    double minLim[3], maxLim[3];
    ellipsoid->getLimits(false, minLim, maxLim); // get real limits
    double limits[] = { maxLim[0] - minLim[0],
                        maxLim[1] - minLim[1],
                        maxLim[2] - minLim[2] }; // lengths of real limits

    getCenter(node, nEllipsoids, valueType, &newCenter, \
              masterBrickStartPoint, masterBrickLengths, limits, procName);
    ellipsoid->translate(newCenter); // translate ellipsoid to new center
    
    // check if this ellipsoid intersects superelement faces or other geometric shapes
    bool intersection = false; // there is no intersection by default
    // intersection with superelement faces
    intersection = ellipsoid->checkIntersection(masterBrickStartPoint, masterBrickLengths);
    // intersection with other geometric shapes
    for (int i = 0; i < (int)elements.size() && !intersection; i++)
      intersection = ellipsoid->checkIntersection(elements.at(i));
    if (!intersection) { // if there is no intersection
      elements.push_back(ellipsoid); // add to list
      nTries = 0; // reset counter of bad tries
      el++; // the number of ellipsoids increases
      nRealEllipsoids++; // real number of ellipsoids
    }
    else
      nTries++; // regular unsuccessful attempt
  }

  return 0; // all is good
}

/**
 * Orthobricks initialization.
 * \param[in] node - node of xml tree containing parameters of orthobricks
 */
int Geometry::orthobricksInit(pugi::xml_node &node) {

  std::string procName = "Geometry::orthobricksInit";

  getNumber(node, &nOrthoBricks, procName); // get the number of orthobricks
  int valueType = getValueType(node, procName); // get the type of values of all
                                                // parameters: absolute (0) or relative (1) (exept cl, that is always absolute)

  // orthobrick parameters
  double lenX[2], lenY[2], lenZ[2];
  getLength(node, "lenX", lenX, valueType, masterBrickLengths[0], procName);
  getLength(node, "lenY", lenY, valueType, masterBrickLengths[1], procName);
  getLength(node, "lenZ", lenZ, valueType, masterBrickLengths[2], procName);
  
  double rvector[3]; // rotation vector
  int rvec = getRVec(node, rvector, procName); // -2 - vector is not defined,
                                               // -1 - vector is defined exactly,
                                               // 0 - random, 1/2/3 - along x/y/z-direction
  double angle[] = { 0, 0 }; // angle of rotation
  if (rvec != -2) // only if rotation vector is defined
    getAngle(node, angle, procName); // get angle (min and max values, or exact one)

  double cl = getCL(node, procName); // characteristic length for mesh building

#ifdef DEBUG
  std::cout << "nOrthoBricks = " << nOrthoBricks << std::endl;
  std::cout << "valueType = " << valueType << std::endl;
  std::cout << "lenX = " << lenX[0] << ", " << lenX[1] << std::endl;
  std::cout << "lenY = " << lenY[0] << ", " << lenY[1] << std::endl;
  std::cout << "lenZ = " << lenZ[0] << ", " << lenZ[1] << std::endl;
  std::cout << "rvec = " << rvec << std::endl;
  std::cout << "angle = " << angle[0] << ", " << angle[1] << std::endl;
  std::cout << "cl = " << cl << std::endl;
#endif

  // orthobricks initialization
  int nTries = 0; // the number of unsuccessful tries
  for (int el = 0; el < nOrthoBricks && nTries < N_TRIES_TO_CREATE_ELEMENT;) {
    // center of orthobrick is
    // at center of superelement
    double xcen = masterBrickStartPoint.getX() + 0.5 * masterBrickLengths[0];
    double ycen = masterBrickStartPoint.getY() + 0.5 * masterBrickLengths[1];
    double zcen = masterBrickStartPoint.getZ() + 0.5 * masterBrickLengths[2];
    // lengths.
    // if first element is negative, then the second element is exact value
    double lx = (lenX[0] > 0 ? randomBetween(lenX[0], lenX[1]) : lenX[1]);
    double ly = (lenY[0] > 0 ? randomBetween(lenY[0], lenY[1]) : lenY[1]);
    double lz = (lenZ[0] > 0 ? randomBetween(lenZ[0], lenZ[1]) : lenZ[1]);
    // rotation vector
    double rv[] = { 0.0, 0.0, 0.0 };
    switch (rvec) {
    case -2: // vector is not defined
      rv[0] = 1.0; // some unit vector (it doesn't change anything because angle is zero in this case)
      break;
    case -1: // vector is defined exactly (by components)
      rv[0] = rvector[0];
      rv[1] = rvector[1];
      rv[2] = rvector[2];
      break;
    case 0: // random vector
      rv[0] = randomBetween(0.0, 1.0);
      rv[1] = randomBetween(0.0, 1.0);
      rv[2] = randomBetween(0.0, 1.0);
      break;
    case 1: // along x-direction
      rv[0] = xcen;
      break;
    case 2: // along y-direction
      rv[1] = ycen;
      break;
    case 3: // along z-direction
      rv[2] = zcen;
      break;
    default:
      require(false, "Unknown type of rotation vector!", procName);
    }
    // angle
    double an = 0.0; // there is no rotation by default
    if (rvec != -2) // if rotation vector is defined
      an = (angle[0] >= 0 ? randomBetween(angle[0], angle[1]) : angle[1]);

    // initialization
    Node3D cen;
    cen.init(xcen, ycen, zcen); // center of orthobrick
    double len[] = { lx, ly, lz }; // lengths
    OrthoBrick *orthobrick = new OrthoBrick;
    orthobrick->init(cen, len, rv, an, cl); // orthobrick initialization

    // translation of the orthobrick:
    // we multiply MARGIN on 0.01 because MARGIN is written in percents regarding to length of superelement side.
    Node3D newCenter;
    double minLim[3], maxLim[3];
    orthobrick->getLimits(false, minLim, maxLim); // get real limits
    double limits[] = { maxLim[0] - minLim[0],
                        maxLim[1] - minLim[1],
                        maxLim[2] - minLim[2] }; // lengths of real limits

    getCenter(node, nOrthoBricks, valueType, &newCenter, \
              masterBrickStartPoint, masterBrickLengths, limits, procName);
    orthobrick->translate(newCenter); // translate orthobrick to new center
    
    // check if this orthobrick intersects superelement faces or other geometric shapes
    bool intersection = false; // there is no intersection by default
    // intersection with superelement faces
    intersection = orthobrick->checkIntersection(masterBrickStartPoint, masterBrickLengths);
    // intersection with other geometric shapes
    for (int i = 0; i < (int)elements.size() && !intersection; i++)
      intersection = orthobrick->checkIntersection(elements.at(i));
    if (!intersection) { // if there is no intersection
      elements.push_back(orthobrick); // add to list
      nTries = 0; // reset counter of bad tries
      el++; // the number of orthobricks increases
      nRealOrthoBricks++; // real number of orthobricks
    }
    else
      nTries++; // regular unsuccessful attempt

#ifdef DEBUG
    std::cout << "center = " << orthobrick->getCenter()->getX() << ", " << orthobrick->getCenter()->getY() << ", " << orthobrick->getCenter()->getZ() << std::endl;
#endif
  }

  return 0; // all is good
}

/**
 * Spheres initialization.
 * \param[in] node - node of xml tree containing parameters of spheres
 */
int Geometry::spheresInit(pugi::xml_node &node) {

  std::string procName = "Geometry::spheresInit";

  getNumber(node, &nSpheres, procName); // get the number of spheres
  int valueType = getValueType(node, procName); // get the type of values of all
                                                // parameters: absolute (0) or relative (1) (exept cl, that is always absolute)

  // sphere parameters
  double diam[2];
  getLength(node, "diam", diam, valueType, \
            (masterBrickLengths[0] + masterBrickLengths[1] + masterBrickLengths[2]) / 3.0, procName);

  double cl = getCL(node, procName); // characteristic length for mesh building

#ifdef DEBUG
  std::cout << "nSpheres = " << nSpheres << std::endl;
  std::cout << "valueType = " << valueType << std::endl;
  std::cout << "diam = " << diam[0] << ", " << diam[1] << std::endl;
  std::cout << "cl = " << cl << std::endl;
#endif

  // spheres initialization
  int nTries = 0; // the number of unsuccessful tries
  for (int el = 0; el < nSpheres && nTries < N_TRIES_TO_CREATE_ELEMENT;) {
    // diameter: 
    // if first element is negative, then the second element is exact value
    double d = (diam[0] > 0 ? randomBetween(diam[0], diam[1]) : diam[1]);
    // rotation vector
    double rv[] = { 1.0, 0.0, 0.0 };
    // angle
    double an = 0.0;

    // initialization
    Node3D center; // center of sphere
    double len[] = { d, d, d }; // lengths
    getCenter(node, nSpheres, valueType, &center, \
              masterBrickStartPoint, masterBrickLengths, len, procName);
    Ellipsoid *sphere = new Ellipsoid;
    sphere->init(center, len, rv, an, cl); // sphere initialization
    
    // check if this sphere intersects superelement faces or other geometric shapes
    bool intersection = false; // there is no intersection by default
    // intersection with superelement faces
    intersection = sphere->checkIntersection(masterBrickStartPoint, masterBrickLengths);
    // intersection with other geometric shapes
    for (int i = 0; i < (int)elements.size() && !intersection; i++)
      intersection = sphere->checkIntersection(elements.at(i));
    if (!intersection) { // if there is no intersection
      elements.push_back(sphere); // add to list
      nTries = 0; // reset counter of bad tries
      el++; // the number of spheres increases
      nRealSpheres++; // real number of spheres
    }
    else
      nTries++; // regular unsuccessful attempt
  }

  return 0; // all is good
}

/**
 * Write geometry to the corresponding file.
 * \param[in] fileName - the name of .geo file where the geometry will be saved
 */
void Geometry::writeResults(std::string fileName) {

  require(masterBrickLengths[0] > 0, "Super brick is not initialized! Please check the parameters file!", "Geometry::writeResults");

  std::ofstream out(fileName.c_str());
  frequire(out, fileName, "Geometry::writeResults");
  
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n\n";
  
  out << "x0 = " << masterBrickStartPoint.getX() << ";\n";
  out << "y0 = " << masterBrickStartPoint.getY() << ";\n";
  out << "z0 = " << masterBrickStartPoint.getZ() << ";\n";
  out << "lenX = " << masterBrickLengths[0] << ";\n";
  out << "lenY = " << masterBrickLengths[1] << ";\n";
  out << "lenZ = " << masterBrickLengths[2] << ";\n";
  out << "clBrick = " << clMasterBrick << ";\n";
  out << "surfNumber = 0;\n";
  out << "volNumber = 0;\n";
  out << "Call MasterBrickTemplate;\n\n";
  
  for (int i = 0; i < nElements; i++) {
    GeoShape *el = elements.at(i); // current element
    out << el->printGeo(i + 1, i) << "\n\n"; // write some information to geo file
  }
  
  out << "vol = newv;\n";
  out << "Volume(vol) = { surfaceLoops[] };\n";
  out << "Physical Surface(" << LEFT_BOUNDARY << ") = { masterBrickLeftFace };\n";
  out << "Physical Surface(" << RIGHT_BOUNDARY << ") = { masterBrickRightFace };\n";
  out << "Physical Surface(" << FRONT_BOUNDARY << ") = { masterBrickFrontFace };\n";
  out << "Physical Surface(" << BACK_BOUNDARY << ") = { masterBrickBackFace };\n";
  out << "Physical Surface(" << BOTTOM_BOUNDARY << ") = { masterBrickBottomFace };\n";
  out << "Physical Surface(" << TOP_BOUNDARY << ") = { masterBrickTopFace };\n";
  out << "Physical Surface(" << INCL_SURFACE << ") = { orthoBrickSurfaces[], ellSurfaces[], cylinderSurfaces[], tetSurfaces[] };\n";
  out << "Physical Volume(" << MAIN_DOMAIN << ") = { vol };\n";
  out << "Physical Volume(" << INCLUSION_DOMAIN << ") = { inclVolumes[] };\n";

  out.close();
}

/**
 * Write geometry of empty superelement to the corresponding file.
 * \param[in] fileName - the name of .geo file where the geometry will be saved
 */
void Geometry::writeEmptyMasterBrick(std::string fileName) {

  require(masterBrickLengths[0] > 0, "Super brick is not initialized! Please check the parameters file!", "Geometry::writeResults");

  std::ofstream out(fileName.c_str());
  frequire(out, fileName, "Geometry::writeEmptyMasterBrick");
  
  out << "Include \"" << TEMPLATES_DIR + TEMPLATES_FILENAME << "\";\n\n";
  
  out << "x0 = " << masterBrickStartPoint.getX() << ";\n";
  out << "y0 = " << masterBrickStartPoint.getY() << ";\n";
  out << "z0 = " << masterBrickStartPoint.getZ() << ";\n";
  out << "lenX = " << masterBrickLengths[0] << ";\n";
  out << "lenY = " << masterBrickLengths[1] << ";\n";
  out << "lenZ = " << masterBrickLengths[2] << ";\n";
  out << "clBrick = " << clMasterBrick << ";\n";
  out << "surfNumber = 0;\n";
  out << "volNumber = 0;\n";
  out << "Call MasterBrickTemplate;\n\n";

  out << "Physical Surface(" << LEFT_BOUNDARY << ") = { masterBrickLeftFace };\n";
  out << "Physical Surface(" << RIGHT_BOUNDARY << ") = { masterBrickRightFace };\n";
  out << "Physical Surface(" << FRONT_BOUNDARY << ") = { masterBrickFrontFace };\n";
  out << "Physical Surface(" << BACK_BOUNDARY << ") = { masterBrickBackFace };\n";
  out << "Physical Surface(" << BOTTOM_BOUNDARY << ") = { masterBrickBottomFace };\n";
  out << "Physical Surface(" << TOP_BOUNDARY << ") = { masterBrickTopFace };\n";
  out << "Physical Volume(" << MAIN_DOMAIN << ") = { masterBrickVolumes[volNumber] };\n";

  out.close();
}

/**
 * Print Gmsh-templates to the file "templates.geo" in 'dir' directory
 * \param[in] dir - the path to the "templates.geo" file
 */
void Geometry::printTemplates(const std::string &dir) {
  std::string fileName = dir + "/" + TEMPLATES_FILENAME;
  std::cout << "  printTemplates" << std::endl;
  std::ofstream out(fileName.c_str());
  frequire(out, fileName, "Geometry::printTemplates");
  
  out << "// template for orthobrick\n";
  out << "// input parameters:\n";
  out << "//           6 --------- 7\n";
  out << "//          /|          /|\n";
  out << "//         / |         / |\n";
  out << "//        2 --------- 3  |\n";
  out << "//        |  4 -------|- 5              Z\n";
  out << "//        | /         | /               |  Y\n";
  out << "//        |/          |/                | /\n";
  out << "//        0 --------- 1                 |/____X\n";
  out << "// xi, yi, zi, i=0..7 - coordinates of the vertices\n";
  out << "// clOrthoBrick - characteristic length\n";
  out << "// surfNumber - the number in surfaceLoops array\n";
  out << "// volNumber - the number in orthoBrickVolumes array\n";

  out << "Function OrthoBrickTemplate\n";

  out << "  p1 = newp; Point(p1) = { x0, y0, z0, cl };\n";
  out << "  p2 = newp; Point(p2) = { x1, y1, z1, cl };\n";
  out << "  p3 = newp; Point(p3) = { x2, y2, z2, cl };\n";
  out << "  p4 = newp; Point(p4) = { x3, y3, z3, cl };\n";
  out << "  p5 = newp; Point(p5) = { x4, y4, z4, cl };\n";
  out << "  p6 = newp; Point(p6) = { x5, y5, z5, cl };\n";
  out << "  p7 = newp; Point(p7) = { x6, y6, z6, cl };\n";
  out << "  p8 = newp; Point(p8) = { x7, y7, z7, cl };\n";

  out << "  l1 = newl; Line(l1) = { p1, p2 };\n";
  out << "  l2 = newl; Line(l2) = { p1, p3 };\n";
  out << "  l3 = newl; Line(l3) = { p2, p4 };\n";
  out << "  l4 = newl; Line(l4) = { p3, p4 };\n";
  out << "  l5 = newl; Line(l5) = { p5, p6 };\n";
  out << "  l6 = newl; Line(l6) = { p5, p7 };\n";
  out << "  l7 = newl; Line(l7) = { p6, p8 };\n";
  out << "  l8 = newl; Line(l8) = { p7, p8 };\n";
  out << "  l9 = newl; Line(l9) = { p1, p5 };\n";
  out << "  l10= newl; Line(l10)= { p2, p6 };\n";
  out << "  l11= newl; Line(l11)= { p3, p7 };\n";
  out << "  l12= newl; Line(l12)= { p4, p8 };\n";

  out << "  ll1 = newll; Line Loop(ll1) = { l2, l11, -l6, -l9 };\n";
  out << "  ll2 = newll; Line Loop(ll2) = { l3, l12, -l7, -l10 };\n";
  out << "  ll3 = newll; Line Loop(ll3) = { l1, l3, -l4, -l2 };\n";
  out << "  ll4 = newll; Line Loop(ll4) = { l5, l7, -l8, -l6 };\n";
  out << "  ll5 = newll; Line Loop(ll5) = { l1, l10, -l5, -l9 };\n";
  out << "  ll6 = newll; Line Loop(ll6) = { l4, l12, -l8, -l11 };\n";

  out << "  s1 = news; Plane Surface(s1) = { ll1 };\n";
  out << "  s2 = news; Plane Surface(s2) = { ll2 };\n";
  out << "  s3 = news; Plane Surface(s3) = { ll3 };\n";
  out << "  s4 = news; Plane Surface(s4) = { ll4 };\n";
  out << "  s5 = news; Plane Surface(s5) = { ll5 };\n";
  out << "  s6 = news; Plane Surface(s6) = { ll6 };\n";

  out << "  sl1 = newsl;\n";
  out << "  Surface Loop(sl1) = { s1, s2, s3, s4, s5, s6 };\n";
  out << "  surfaceLoops[surfNumber] = sl1;\n";

  out << "  v1 = newv;\n";
  out << "  Volume(v1) = { sl1 };\n";
  out << "  inclVolumes[volNumber] = v1;\n";
  
  out << "  orthoBrickSurfaces[6 * volNumber] = s1;\n";
  out << "  orthoBrickSurfaces[6 * volNumber + 1] = s2;\n";
  out << "  orthoBrickSurfaces[6 * volNumber + 2] = s3;\n";
  out << "  orthoBrickSurfaces[6 * volNumber + 3] = s4;\n";
  out << "  orthoBrickSurfaces[6 * volNumber + 4] = s5;\n";
  out << "  orthoBrickSurfaces[6 * volNumber + 5] = s6;\n";

  out << "Return\n";

  out << "// ----------------------------------------------\n";

  out << "// template for master brick (brick with edges parallel to Cartesian coordinates axes)\n";
  out << "// input parameters:\n";
  out << "// x0, y0, z0 - brick start point\n";
  out << "// lenX, lenY, LenZ - lengths of brick sides\n";
  out << "// clBrick - characteristic length\n";
  out << "// surfNumber - the number in surfaceLoops array\n";
  out << "// volNumber - the number in brickVolumes array\n";

  out << "Function MasterBrickTemplate\n";

  out << "  p1 = newp; Point(p1) = { x0, y0, z0, clBrick };\n";
  out << "  p2 = newp; Point(p2) = { x0 + lenX, y0, z0, clBrick };\n";
  out << "  p3 = newp; Point(p3) = { x0, y0, z0 + lenZ, clBrick };\n";
  out << "  p4 = newp; Point(p4) = { x0 + lenX, y0, z0 + lenZ, clBrick };\n";
  out << "  p5 = newp; Point(p5) = { x0, y0 + lenY, z0, clBrick };\n";
  out << "  p6 = newp; Point(p6) = { x0 + lenX, y0 + lenY, z0, clBrick };\n";
  out << "  p7 = newp; Point(p7) = { x0, y0 + lenY, z0 + lenZ, clBrick };\n";
  out << "  p8 = newp; Point(p8) = { x0 + lenX, y0 + lenY, z0 + lenZ, clBrick };\n";

  out << "  l1 = newl; Line(l1) = { p1, p2 };\n";
  out << "  l2 = newl; Line(l2) = { p1, p3 };\n";
  out << "  l3 = newl; Line(l3) = { p2, p4 };\n";
  out << "  l4 = newl; Line(l4) = { p3, p4 };\n";
  out << "  l5 = newl; Line(l5) = { p5, p6 };\n";
  out << "  l6 = newl; Line(l6) = { p5, p7 };\n";
  out << "  l7 = newl; Line(l7) = { p6, p8 };\n";
  out << "  l8 = newl; Line(l8) = { p7, p8 };\n";
  out << "  l9 = newl; Line(l9) = { p1, p5 };\n";
  out << "  l10= newl; Line(l10)= { p2, p6 };\n";
  out << "  l11= newl; Line(l11)= { p3, p7 };\n";
  out << "  l12= newl; Line(l12)= { p4, p8 };\n";

  out << "  ll1 = newll; Line Loop(ll1) = { l2, l11, -l6, -l9 };\n";
  out << "  ll2 = newll; Line Loop(ll2) = { l3, l12, -l7, -l10 };\n";
  out << "  ll3 = newll; Line Loop(ll3) = { l1, l3, -l4, -l2 };\n";
  out << "  ll4 = newll; Line Loop(ll4) = { l5, l7, -l8, -l6 };\n";
  out << "  ll5 = newll; Line Loop(ll5) = { l1, l10, -l5, -l9 };\n";
  out << "  ll6 = newll; Line Loop(ll6) = { l4, l12, -l8, -l11 };\n";

  out << "  s1 = news; Plane Surface(s1) = { ll1 };\n";
  out << "  s2 = news; Plane Surface(s2) = { ll2 };\n";
  out << "  s3 = news; Plane Surface(s3) = { ll3 };\n";
  out << "  s4 = news; Plane Surface(s4) = { ll4 };\n";
  out << "  s5 = news; Plane Surface(s5) = { ll5 };\n";
  out << "  s6 = news; Plane Surface(s6) = { ll6 };\n";

  out << "  sl1 = newsl;\n";
  out << "  Surface Loop(sl1) = { s1, s2, s3, s4, s5, s6 };\n";
  out << "  surfaceLoops[surfNumber] = sl1;\n";

  out << "  v1 = newv;\n";
  out << "  Volume(v1) = { sl1 };\n";
  out << "  masterBrickVolumes[volNumber] = v1;\n";
  
  out << "  masterBrickLeftFace = s1;\n";
  out << "  masterBrickRightFace = s2;\n";
  out << "  masterBrickFrontFace = s3;\n";
  out << "  masterBrickBackFace = s4;\n";
  out << "  masterBrickBottomFace = s5;\n";
  out << "  masterBrickTopFace = s6;\n";

  out << "Return\n";

  out << "// ----------------------------------------------\n";

  out << "// template for tetrahedron\n";
  out << "// input parameters:\n";
  out << "// x0, y0, z0 - vertex 1\n";
  out << "// x1, y1, z1 - vertex 2\n";
  out << "// x2, y2, z2 - vertex 3\n";
  out << "// x3, y3, z3 - vertex 4\n";
  out << "// cl - characteristic length\n";
  out << "// surfNumber - the number in surfaceLoops array\n";
  out << "// volNumber - the number in tetVolumes array\n";

  out << "Function TetrahedronTemplate\n";
  out << "  p1 = newp; Point(p1) = { x0, y0, z0, cl };\n";
  out << "  p2 = newp; Point(p2) = { x1, y1, z1, cl };\n";
  out << "  p3 = newp; Point(p3) = { x2, y2, z2, cl };\n";
  out << "  p4 = newp; Point(p4) = { x3, y3, z3, cl };\n";

  out << "  l1 = newl; Line(l1) = { p1, p2 };\n";
  out << "  l2 = newl; Line(l2) = { p1, p3 };\n";
  out << "  l3 = newl; Line(l3) = { p1, p4 };\n";
  out << "  l4 = newl; Line(l4) = { p2, p3 };\n";
  out << "  l5 = newl; Line(l5) = { p2, p4 };\n";
  out << "  l6 = newl; Line(l6) = { p3, p4 };\n";

  out << "  ll1 = newll; Line Loop(ll1) = { l1, l4, -l2 };\n";
  out << "  ll2 = newll; Line Loop(ll2) = { l1, l5, -l3 };\n";
  out << "  ll3 = newll; Line Loop(ll3) = { l2, l6, -l3 };\n";
  out << "  ll4 = newll; Line Loop(ll4) = { l4, l6, -l5 };\n";

  out << "  s1 = news; Plane Surface(s1) = { ll1 };\n";
  out << "  s2 = news; Plane Surface(s2) = { ll2 };\n";
  out << "  s3 = news; Plane Surface(s3) = { ll3 };\n";
  out << "  s4 = news; Plane Surface(s4) = { ll4 };\n";

  out << "  sl1 = newsl;\n";
  out << "  Surface Loop(sl1) = { s1, s2, s3, s4 };\n";
  out << "  surfaceLoops[surfNumber] = sl1;\n";

  out << "  v1 = newv;\n";
  out << "  Volume(v1) = { sl1 };\n";
  out << "  inclVolumes[volNumber] = v1;\n";
  
  out << "  tetSurfaces[4 * volNumber] = s1;\n";
  out << "  tetSurfaces[4 * volNumber + 1] = s2;\n";
  out << "  tetSurfaces[4 * volNumber + 2] = s3;\n";
  out << "  tetSurfaces[4 * volNumber + 3] = s4;\n";

  out << "Return\n";
  
  out << "// ----------------------------------------------\n";
  out << "// template for ellipsoid\n";
  out << "// input parameters:\n";
  out << "// x0, y0, z0 - left point\n";
  out << "// x1, y1, z1 - right point\n";
  out << "// x2, y2, z2 - front point\n";
  out << "// x3, y3, z3 - back point\n";
  out << "// x4, y4, z4 - bottom point\n";
  out << "// x5, y5, z5 - top point\n";
  out << "// clEll - characteristic length\n";
  out << "// surfNumber - the number in surfaceLoops array\n";
  out << "// volNumber - the number in sphVolumes array\n";

  out << "Function EllipsoidTemplate\n";

  out << "  xcen = (x0 + x1 + x2 + x3 + x4 + x5) / 6.0;\n";
  out << "  ycen = (y0 + y1 + y2 + y3 + y4 + y5) / 6.0;\n";
  out << "  zcen = (z0 + z1 + z2 + z3 + z4 + z5) / 6.0;\n";

  out << "  p1 = newp; Point(p1) = { xcen, ycen, zcen, cl };\n";
  out << "  p2 = newp; Point(p2) = { x5, y5, z5, cl };\n";
  out << "  p3 = newp; Point(p3) = { x0, y0, z0, cl };\n";
  out << "  p4 = newp; Point(p4) = { x4, y4, z4, cl };\n";
  out << "  p5 = newp; Point(p5) = { x1, y1, z1, cl };\n";
  out << "  p6 = newp; Point(p6) = { x2, y2, z2, cl };\n";
  out << "  p7 = newp; Point(p7) = { x3, y3, z3, cl };\n";

  out << "  l1 = newl; Ellipse(l1) = { p2, p1, p5, p3 };\n";
  out << "  l2 = newl; Ellipse(l2) = { p3, p1, p5, p4 };\n";
  out << "  l3 = newl; Ellipse(l3) = { p4, p1, p5, p5 };\n";
  out << "  l4 = newl; Ellipse(l4) = { p5, p1, p5, p2 };\n";
  out << "  l5 = newl; Ellipse(l5) = { p3, p1, p5, p6 };\n";
  out << "  l6 = newl; Ellipse(l6) = { p6, p1, p5, p5 };\n";
  out << "  l7 = newl; Ellipse(l7) = { p5, p1, p5, p7 };\n";
  out << "  l8 = newl; Ellipse(l8) = { p7, p1, p5, p3 };\n";
  out << "  l9 = newl; Ellipse(l9) = { p2, p1, p6, p6 };\n";
  out << "  l10= newl; Ellipse(l10)= { p6, p1, p6, p4 };\n";
  out << "  l11= newl; Ellipse(l11)= { p4, p1, p6, p7 };\n";
  out << "  l12= newl; Ellipse(l12)= { p7, p1, p6, p2 };\n";

  out << "  ll1 = newll; Line Loop(ll1) = { l1, l5, -l9 };\n";
  out << "  ll2 = newll; Line Loop(ll2) = { l4, l9, l6 };\n";
  out << "  ll3 = newll; Line Loop(ll3) = { l2, -l10, -l5 };\n";
  out << "  ll4 = newll; Line Loop(ll4) = { l3, -l6, l10 };\n";
  out << "  ll5 = newll; Line Loop(ll5) = { l1, -l8, l12 };\n";
  out << "  ll6 = newll; Line Loop(ll6) = { l4, -l12, -l7 };\n";
  out << "  ll7 = newll; Line Loop(ll7) = { l2, l11, l8 };\n";
  out << "  ll8 = newll; Line Loop(ll8) = { l3, l7, -l11 };\n";

  out << "  s1 = news; Ruled Surface(s1) = { ll1 };\n";
  out << "  s2 = news; Ruled Surface(s2) = { ll2 };\n";
  out << "  s3 = news; Ruled Surface(s3) = { ll3 };\n";
  out << "  s4 = news; Ruled Surface(s4) = { ll4 };\n";
  out << "  s5 = news; Ruled Surface(s5) = { ll5 };\n";
  out << "  s6 = news; Ruled Surface(s6) = { ll6 };\n";
  out << "  s7 = news; Ruled Surface(s7) = { ll7 };\n";
  out << "  s8 = news; Ruled Surface(s8) = { ll8 };\n";

  out << "  sl1 = newsl;\n";
  out << "  Surface Loop(sl1) = { s1, s2, s3, s4, s5, s6, s7, s8 };\n";
  out << "  surfaceLoops[surfNumber] = sl1;\n";

  out << "  v1 = newv;\n";
  out << "  Volume(v1) = { sl1 };\n";
  out << "  inclVolumes[volNumber] = v1;\n";

  out << "  ellSurfaces[8 * volNumber] = s1;\n";
  out << "  ellSurfaces[8 * volNumber + 1] = s2;\n";
  out << "  ellSurfaces[8 * volNumber + 2] = s3;\n";
  out << "  ellSurfaces[8 * volNumber + 3] = s4;\n";
  out << "  ellSurfaces[8 * volNumber + 4] = s5;\n";
  out << "  ellSurfaces[8 * volNumber + 5] = s6;\n";
  out << "  ellSurfaces[8 * volNumber + 6] = s7;\n";
  out << "  ellSurfaces[8 * volNumber + 7] = s8;\n";

  out << "Return\n";
  
  out << "// ----------------------------------------------\n";
  out << "// template for cylinder\n";
  out << "// xi, yi, zi, i = 1,..,8\n";
  out << "// clCylinder\n";
  out << "// surfNumber\n";
  out << "// volNumber\n";
  out << "Function CylinderTemplate\n";
  out << "  xcen1 = (x0 + x1 + x2 + x3) / 4.0;\n";
  out << "  ycen1 = (y0 + y1 + y2 + y3) / 4.0;\n";
  out << "  zcen1 = (z0 + z1 + z2 + z3) / 4.0;\n";
  out << "  xcen2 = (x4 + x5 + x6 + x7) / 4.0;\n";
  out << "  ycen2 = (y4 + y5 + y6 + y7) / 4.0;\n";
  out << "  zcen2 = (z4 + z5 + z6 + z7) / 4.0;\n";
  out << "  p0 = newp; Point(p0) = { x0, y0, z0, cl };\n";
  out << "  p1 = newp; Point(p1) = { x1, y1, z1, cl };\n";
  out << "  p2 = newp; Point(p2) = { x2, y2, z2, cl };\n";
  out << "  p3 = newp; Point(p3) = { x3, y3, z3, cl };\n";
  out << "  p4 = newp; Point(p4) = { x4, y4, z4, cl };\n";
  out << "  p5 = newp; Point(p5) = { x5, y5, z5, cl };\n";
  out << "  p6 = newp; Point(p6) = { x6, y6, z6, cl };\n";
  out << "  p7 = newp; Point(p7) = { x7, y7, z7, cl };\n";
  out << "  p11 =newp; Point(p11)= { xcen1, ycen1, zcen1, cl };\n";
  out << "  p21 =newp; Point(p21)= { xcen2, ycen2, zcen2, cl };\n";
  out << "  l1 = newl; Ellipse(l1) = { p0, p11, p0, p2 };\n";
  out << "  l2 = newl; Ellipse(l2) = { p2, p11, p0, p1 };\n";
  out << "  l3 = newl; Ellipse(l3) = { p1, p11, p0, p3 };\n";
  out << "  l4 = newl; Ellipse(l4) = { p3, p11, p0, p0 };\n";
  out << "  l5 = newl; Ellipse(l5) = { p4, p21, p4, p6 };\n";
  out << "  l6 = newl; Ellipse(l6) = { p6, p21, p4, p5 };\n";
  out << "  l7 = newl; Ellipse(l7) = { p5, p21, p4, p7 };\n";
  out << "  l8 = newl; Ellipse(l8) = { p7, p21, p4, p4 };\n";
  out << "  l9 = newl; Line(l9) = { p0, p4 };\n";
  out << "  l10= newl; Line(l10)= { p2, p6 };\n";
  out << "  l11= newl; Line(l11)= { p1, p5 };\n";
  out << "  l12= newl; Line(l12)= { p3, p7 };\n";
  out << "  ll1 = newll; Line Loop(ll1) = { l1, l10, -l5, -l9 };\n";
  out << "  ll2 = newll; Line Loop(ll2) = { l2, l11, -l6, -l10 };\n";
  out << "  ll3 = newll; Line Loop(ll3) = { l3, l12, -l7, -l11 };\n";
  out << "  ll4 = newll; Line Loop(ll4) = { l4, l9, -l8, -l12 };\n";
  out << "  ll5 = newll; Line Loop(ll5) = { l1, l2, l3, l4 };\n";
  out << "  ll6 = newll; Line Loop(ll6) = { l5, l6, l7, l8 };\n";
  out << "  s1 = news; Ruled Surface(s1) = { ll1 };\n";
  out << "  s2 = news; Ruled Surface(s2) = { ll2 };\n";
  out << "  s3 = news; Ruled Surface(s3) = { ll3 };\n";
  out << "  s4 = news; Ruled Surface(s4) = { ll4 };\n";
  out << "  s5 = news; Plane Surface(s5) = { ll5 };\n";
  out << "  s6 = news; Plane Surface(s6) = { ll6 };\n";
  out << "  sl1 = newsl;\n";
  out << "  Surface Loop(sl1) = { s1, s2, s3, s4, s5, s6 };\n";
  out << "  surfaceLoops[surfNumber] = sl1;\n";
  out << "  v1 = newv;\n";
  out << "  Volume(v1) = { sl1 };\n";
  out << "  inclVolumes[volNumber] = v1;\n";
  out << "  cylinderSurfaces[6 * volNumber] = s1;\n";
  out << "  cylinderSurfaces[6 * volNumber + 1] = s2;\n";
  out << "  cylinderSurfaces[6 * volNumber + 2] = s3;\n";
  out << "  cylinderSurfaces[6 * volNumber + 3] = s4;\n";
  out << "  cylinderSurfaces[6 * volNumber + 4] = s5;\n";
  out << "  cylinderSurfaces[6 * volNumber + 5] = s6;\n";
  
  out << "Return\n";
  
  out.close();
}

/**
 * @brief Geometry::getnOrthoBricks
 * @return The number of orthobricks that we would like to have (what was written in xml file).
 */
int Geometry::getnOrthoBricks() const { return nOrthoBricks; }

/**
 * @brief Geometry::getnRealOrthoBricks
 * @return The number of orthobricks that we really have. This number may differ from the one written in xml file,
 * if the desire number is too big, and we can't find a place where to insert new inclusion.
 */
int Geometry::getnRealOrthoBricks() const { return nRealOrthoBricks; }

/**
 * @brief Geometry::getnCylinders
 * @return The number of cylinders that we would like to have (what was written in xml file).
 */
int Geometry::getnCylinders() const { return nCylinders; }

/**
 * @brief Geometry::getnRealCylinders
 * @return The number of cylinders that we really have. This number may differ from the one written in xml file,
 * if the desire number is too big, and we can't find a place where to insert new inclusion.
 */
int Geometry::getnRealCylinders() const { return nRealCylinders; }

/**
 * @brief Geometry::getnEllipsoids
 * @return The number of ellipsoids that we would like to have (what was written in xml file).
 */
int Geometry::getnEllipsoids() const { return nEllipsoids; }

/**
 * @brief Geometry::getnRealEllipsoids
 * @return The number of ellipsoids that we really have. This number may differ from the one written in xml file,
 * if the desire number is too big, and we can't find a place where to insert new inclusion.
 */
int Geometry::getnRealEllipsoids() const { return nRealEllipsoids; }

/**
 * @brief Geometry::getnSpheres
 * @return The number of spheres that we would like to have (what was written in xml file).
 */
int Geometry::getnSpheres() const { return nSpheres; }

/**
 * @brief Geometry::getnRealSpheres
 * @return The number of spheres that we really have. This number may differ from the one written in xml file,
 * if the desire number is too big, and we can't find a place where to insert new inclusion.
 */
int Geometry::getnRealSpheres() const { return nRealSpheres; }

/**
 * @brief Geometry::getnTetrahedra
 * @return The number of tetrahedra that we would like to have (what was written in xml file).
 */
int Geometry::getnTetrahedra() const { return nTetrahedra; }

/**
 * @brief Geometry::getnRealTetrahedra
 * @return The number of tetrahedra that we really have. This number may differ from the one written in xml file,
 * if the desire number is too big, and we can't find a place where to insert new inclusion.
 */
int Geometry::getnRealTetrahedra() const { return nRealTetrahedra; }

/**
 * To walk on xml document checking what type of inclusions we should create.
 */
bool Geometry::for_each(pugi::xml_node &node) {

  std::string geoType = (std::string)node.name(); // type of geometric shape

#ifdef DEBUG
  std::cout << "type = " << node.type() << " name = " << node.name() << std::endl;
#endif

  if (geoType == "orthobrick")
    orthobricksInit(node);
  else if (geoType == "cylinder")
    cylindersInit(node);
  else if (geoType == "ellipsoid")
    ellipsoidsInit(node);
  else if (geoType == "sphere")
    spheresInit(node);
  else if (geoType == "tetrahedron")
    tetrahedraInit(node);
  else if (geoType == "superbrick")
    masterbrickInit(node);
//  else
//    require(false, "Unknown type of geometric shape ('" + geoType + "')!", "Geometry::for_each");
  return true;
}

/**
 * Get the number of elements (of partiacular inclusion type).
 * It's necessary attribute of inclusion in xml file.
 * \param[in] node - reference to node of xml tree
 * \param[out] number - the number itself (that follows from the name of function)
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
void Geometry::getNumber(pugi::xml_node &node, int *number, std::string procName) {
  pugi::xml_attribute attr = node.attribute("number"); // necessary attribute
  require(attr != 0, "There is no attribute 'number'!", procName);
  *number = attr.as_int(); // desirable number of the elements of this type
  require(*number > 0, "The number is less than 1!", procName);
}

/**
 * Get the type of values (absolute or relative).
 * It's optional attribute: absolute values are set by default
 * \param[in] node - reference to node of xml tree
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
int Geometry::getValueType(pugi::xml_node &node, std::string procName) {
  pugi::xml_attribute attr = node.attribute("vtype"); // optional attribute
  int vType = 0; // absolute values by default
  if (attr) { // if such attribute exists
    if((std::string)attr.value() == "absolute") // and values are absolute
      vType = 0; // so we have absolute values
    else if((std::string)attr.value() == "relative") // or they are relative
      vType = 1; // so we have relative values
    else
      require(false, "Unknown type of values", procName);
  }
  return vType;
}

/**
 * Get the length under different names.
 * It's necessary attribute
 * \param[in] node - reference to node of xml tree
 * \param[in] lengthName - lengths might be called variously, so we need to know the name of attribute
 * \param[out] lengths - the array of lengths
 * \param[in] valueType - global type of values (absolute or relative)
 * \param[in] superLength - the length of master brick
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
void Geometry::getLength(pugi::xml_node &node, std::string lengthName, double lengths[], \
                         int valueType, double superLength, std::string procName) {
  pugi::xml_node child = node.child(lengthName.c_str());
  require(child, "There is no child '" + lengthName + "'!", procName);

  // we should have min and max values to create length between them, or exact value
  pugi::xml_attribute attr_min = child.attribute("min"); // minimal value for this length. necessary
  pugi::xml_attribute attr_max = child.attribute("max"); // maximal value for this length. necessary
  pugi::xml_attribute attr_exact = child.attribute("exact"); // exact value. necessary
  // but if all (min, max and exact) values are defined, this is an error
  if (attr_exact)
    require(!attr_min && !attr_max, "Set of attributes is incorrect (min and/or max attribute exists with exact value) for child '" + lengthName + "'!", procName);
  else
    require(attr_min && attr_max, "There are too few values (min, max or exact) for child '" + lengthName + "'!", procName);

  pugi::xml_attribute attr_vt = child.attribute("vtype"); // value type specifically for this length. optional - valueType by default
  int vt = valueType; // general type of values by default
  if (attr_vt) {
    if ((std::string)attr_vt.value() == "absolute")
      vt = 0;
    else if ((std::string)attr_vt.value() == "relative")
      vt = 1;
    else
      require(false, "Unknown value type for child '" + lengthName + "'!", procName);
  }

  if (attr_exact) {
    double exact = attr_exact.as_double(); // exact value
    require(exact > 0, "Exact value is negative for child '" + lengthName + "'!", procName);
    lengths[0] = -1; // it means that the second element in array is exact value
    lengths[1] = (vt ? exact * 0.01 * superLength : exact); // exact length
  }
  else {
    double vmin = attr_min.as_double(); // minimal value
    double vmax = attr_max.as_double(); // maximal value
    require(vmin <= vmax && vmin > 0, "'min' is more than 'max' for child '" + lengthName + "' or is negative!", procName);
    lengths[0] = (vt ? vmin * 0.01 * superLength : vmin); // minimal length
    lengths[1] = (vt ? vmax * 0.01 * superLength : vmax); // maximal length
  }
}

/**
 * Get the direction of rotation vector:
 * -2 - vector is not defined
 * -1 - vector is defined exactly (by components)
 * 0 - random direction
 * ( 1 | 2 | 3 ) - along ( x | y | z ) - direction
 * \param[in] node - reference to node of xml tree
 * \param[out] vector - rotation vector (array of components)
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
int Geometry::getRVec(pugi::xml_node &node, double vector[], std::string procName) {
  pugi::xml_node child = node.child("rvector");
  if (!child)
    return -2; // rotation vector is not defined

  pugi::xml_attribute attr_x = child.attribute("x"); // x-component of the vector
  pugi::xml_attribute attr_y = child.attribute("y"); // y-component of the vector
  pugi::xml_attribute attr_z = child.attribute("z"); // z-component of the vector
  pugi::xml_attribute attr_dir = child.attribute("direction"); // some general direction
  if (attr_x && attr_y && attr_z) {
    require(!attr_dir, "Set of attributes is incorrect (x, y, z -components exist with direction value) for child 'rvector'!", procName);
    vector[0] = attr_x.as_double();
    vector[1] = attr_y.as_double();
    vector[2] = attr_z.as_double();
    return -1; // rotation vector is defined exactly
  }

  require(attr_dir, "There is no necessary set of parameters for child 'rvector'!", procName);
  require(!attr_x && !attr_y && !attr_z, "Set of attributes is incorrect (x, y, z -components exist with direction value) for child 'rvector'!", procName);
  int rv = attr_dir.as_int(); // direction of vector (0 - random, 1/2/3 - along x/y/z axis)
  require(rv >= 0 && rv <= 3, "Unknown value of 'direction' attribute for child 'rvector'!", procName);
  return rv;
}

/**
 * Get the values of angle of rotation.
 * We call this procedure only if rotation vector is defined,
 * so one can say that it's necessary attribute
 * \param[in] node - reference to node of xml tree
 * \param[out] angle - angle of rotation (it's array because angle might be random value, and we have to know the limits of the distribution)
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
void Geometry::getAngle(pugi::xml_node &node, double angle[], std::string procName) {
  pugi::xml_node child = node.child("angle");
  require(child, "There is no child 'angle'!", procName);

  // we should have min and max values to create angle between them, or exact value
  pugi::xml_attribute attr_min = child.attribute("min"); // minimal value for this angle. necessary
  pugi::xml_attribute attr_max = child.attribute("max"); // maximal value for this angle. necessary
  pugi::xml_attribute attr_exact = child.attribute("exact"); // exact value. necessary
  // but if all (min, max and exact) values are defined, this is an error
  if (attr_exact)
    require(!attr_min && !attr_max, "Set of attributes is incorrect (min and/or max attribute exists with exact value) for child 'angle'!", procName);
  else
    require(attr_min && attr_max, "There are too few values (min, max or exact) for child 'angle'!", procName);

  if (attr_exact) {
    double exact = attr_exact.as_double(); // exact value
    angle[0] = -1; // it means that the second element in array is exact value
    angle[1] = exact; // exact angle
  }
  else {
    double vmin = attr_min.as_double(); // minimal value
    double vmax = attr_max.as_double(); // maximal value
    require(vmin <= vmax && vmin >= 0, \
            "'min' is more than 'max' for child 'angle' or is negative: \
            (when you use min and max values to define angle, you should use only positive values)!", procName);
    angle[0] = vmin; // minimal angle
    angle[1] = vmax; // maximal angle
  }
}

/**
 * Get characteristic length for mesh building.
 * It's necessary value and always absolute (but really it becames relative in GeoShape.init procedure, see code)
 * \param[in] node - reference to node of xml tree
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
double Geometry::getCL(pugi::xml_node &node, std::string procName) {
  pugi::xml_attribute attr = node.attribute("cl");
  require(attr, "There is no 'cl' attribute!", procName);
  double cl = attr.as_double();
  require(cl > 1e-6 && cl <= 1, "'cl' is wrong (must be '> 1e-6' and '<= 1', but here is " + d2s(cl) + ")!", procName);
  return cl;
}

/**
 * Get the center of element.
 * \param[in] node - reference to node of xml tree
 * \param[in] nElements - the number of inclusions of this type (in one xml node)
 * \param[in] valueType - global type of values (absolute or relative)
 * \param[out] center - the center of element itself
 * \param[in] startPoint - the start point of master brick
 * \param[in] superLengths - the lengths of sides of master brick (array)
 * \param[in] limits - the lengths of area taken by inclusion
 * \param[in] procName - the name of procedure from where the function was called (for errors handling)
 */
void Geometry::getCenter(pugi::xml_node &node, int nElements, int valueType, Node3D *center, \
                         Node3D startPoint, double superLengths[], double limits[], std::string procName) {
  pugi::xml_node child = node.child("center");
  double xc, yc, zc; // center coordinates
  // all coordinates of new center we create randomly by default
  xc = randomBetween(startPoint.getX(), startPoint.getX() + superLengths[0], \
                     0.01 * MARGIN * superLengths[0] + 0.5 * limits[0]);
  yc = randomBetween(startPoint.getY(), startPoint.getY() + superLengths[1], \
                     0.01 * MARGIN * superLengths[1] + 0.5 * limits[1]);
  zc = randomBetween(startPoint.getZ(), startPoint.getZ() + superLengths[2], \
                     0.01 * MARGIN * superLengths[2] + 0.5 * limits[2]);
  if (!child) { // if center is not defined
    center->init(xc, yc, zc);
    return;
  }

  pugi::xml_attribute attr_x = child.attribute("x");
  pugi::xml_attribute attr_y = child.attribute("y");
  pugi::xml_attribute attr_z = child.attribute("z");
  // if all coordinates are defined and number of elements is more than 1, there is an error
  require(!(attr_x && attr_y && attr_z && nElements > 1), \
          "Exact center could be defined only for one element (but here are " + d2s(nElements) + " ones)!", procName);

  pugi::xml_attribute attr_vt = child.attribute("vtype"); // value type specifically for the center point. optional - valueType by default
  int vt = valueType; // general type of values by default
  if (attr_vt) {
    if ((std::string)attr_vt.value() == "absolute")
      vt = 0;
    else if ((std::string)attr_vt.value() == "relative")
      vt = 1;
    else
      require(false, "Unknown value type for child 'center'!", procName);
  }

  if (attr_x) // if x-coordinate is defined
    xc = (vt ? startPoint.getX() + attr_x.as_double() * 0.01 * superLengths[0] : attr_x.as_double());
  if (attr_y) // if y-coordinate is defined
    yc = (vt ? startPoint.getY() + attr_y.as_double() * 0.01 * superLengths[1] : attr_y.as_double());
  if (attr_z) // if z-coordinate is defined
    zc = (vt ? startPoint.getZ() + attr_z.as_double() * 0.01 * superLengths[2] : attr_z.as_double());

  center->init(xc, yc, zc); // center initialization
}

std::string getLibraryVersion() {
  //return d2s(mixture_VERSION_MAJOR) + "." + d2s(mixture_VERSION_MINOR) + "." + d2s(mixture_VERSION_PATCH);
  return (std::string)MIXTURE_VERSION_MAJOR + "." + (std::string)MIXTURE_VERSION_MINOR + "." + (std::string)MIXTURE_VERSION_PATCH;
}

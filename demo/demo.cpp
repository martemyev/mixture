// example of using 'mixture' library
// demo.cpp

#include <iostream>
#include <string>

// this header contains the description of all functions that we'll need
// to create the .geo file
#include "Geometry.h"

int main(int argc, char** argv)
{
  // parameters of the parallelepiped
  const double x0   = 0; // start point
  const double y0   = 0;
  const double z0   = 0;
  const double lenX = 10; // size of the shape
  const double lenY = 20;
  const double lenZ = 30;
  // characteristic length of the mesh in the parallelepiped
  const double cl   = 0.1;

  // the name of file containing information about inclusions
  std::string incl_description;
  if (argc >= 2)
    incl_description = std::string(argv[1]);
  else
    incl_description = "xml/inclusions.xml";

  // the resulting file that can be used as input file for Gmsh
  std::string res_geometry;
  if (argc >= 3)
    res_geometry = std::string(argv[2]);
  else
    res_geometry = "geometry.geo";

  Geometry geo;
  geo.init(x0, y0, z0, lenX, lenY, lenZ, cl, incl_description, res_geometry);
  // after creating inclusions it's good practice to check - how many inclusions were created
  std::cout << "we have "           << geo.getnRealSpheres() << \
               " spheres from "     << geo.getnSpheres() << \
               " that we've wanted" << std::endl;
  std::cout << "we have "           << geo.getnRealCylinders() << \
               " cylinders from "   << geo.getnCylinders() << \
               " that we've wanted" << std::endl;
  std::cout << "we have "           << geo.getnRealOrthoBricks() << \
               " bricks from "      << geo.getnOrthoBricks() << \
               " that we've wanted" << std::endl;
  return 0;
}

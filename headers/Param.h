#ifndef PARAM_H
#define PARAM_H

#include <iostream>
#include <string>

// parameters for task
class Param {
private:
  void testParam(); // testing input parameters

public:
  // domains in mesh
  // points are numbered from 1001
  // lines are numbered from 2001
  // surfaces are numbered from 3001
  // volumes are numbered from 4001
  static const int LEFT_FRONT_LINE = 2001;
  static const int LEFT_BACK_LINE = 2002;
  static const int LEFT_BOTTOM_LINE = 2003;
  static const int LEFT_TOP_LINE = 2004;
  
  static const int RIGHT_FRONT_LINE = 2011;
  static const int RIGHT_BACK_LINE = 2012;
  static const int RIGHT_BOTTOM_LINE = 2013;
  static const int RIGHT_TOP_LINE = 2014;
  
  static const int FRONT_BOTTOM_LINE = 2021;
  static const int FRONT_TOP_LINE = 2022;
  static const int BACK_BOTTOM_LINE = 2031;
  static const int BACK_TOP_LINE = 2032;
  
  static const int LEFT_FRONT_LINE_INCL = 2101;
  static const int LEFT_BACK_LINE_INCL = 2102;
  static const int LEFT_BOTTOM_LINE_INCL = 2103;
  static const int LEFT_TOP_LINE_INCL = 2104;
  
  static const int RIGHT_FRONT_LINE_INCL = 2111;
  static const int RIGHT_BACK_LINE_INCL = 2112;
  static const int RIGHT_BOTTOM_LINE_INCL = 2113;
  static const int RIGHT_TOP_LINE_INCL = 2114;
  
  static const int FRONT_BOTTOM_LINE_INCL = 2121;
  static const int FRONT_TOP_LINE_INCL = 2122;
  static const int BACK_BOTTOM_LINE_INCL = 2131;
  static const int BACK_TOP_LINE_INCL = 2132;
  
  static const int LEFT_BOUNDARY = 3001;
  static const int RIGHT_BOUNDARY = 3002;
  static const int FRONT_BOUNDARY = 3003;
  static const int BACK_BOUNDARY = 3004;
  static const int BOTTOM_BOUNDARY = 3005;
  static const int TOP_BOUNDARY = 3006;
  
  static const int LEFT_INCL = 3011;
  static const int RIGHT_INCL = 3012;
  static const int FRONT_INCL = 3013;
  static const int BACK_INCL = 3014;
  static const int BOTTOM_INCL = 3015;
  static const int TOP_INCL = 3016;
  
  static const int INCL_SURFACE = 3021;
  
  static const int MAIN_DOMAIN = 4001;
  static const int INCLUSION_DOMAIN = 4002;

  // constant parameters
  static const int N_DOMAINS = 2; // the number of domains is always 2 (main domain and inclusions)
  int DOMAINS[N_DOMAINS]; // domains (numbers)
  int SLAE_MAXITER; // for SLAE solver
  double SLAE_TOLERANCE; // for SLAE solver
  std::string REPORT_FILE_NAME;
  std::string INCLUSIONS_FILE_NAME;
  double TET_SEARCH_TOLERANCE;
  double BRICK_SEARCH_TOLERANCE;
  std::string GEO_PARAM_NAME; // the name of the geometry parameters file without extension and path to a geometry directory
  std::string GEO_DIR; // directory for geometry directories, templates, etc
  std::string GEO_FULL_DIR; // directory for geometry files
  std::string MSH_DIR; // directory for mesh directories
  std::string MSH_FULL_DIR; // directory for mesh files
  std::string OUTPUT_DIR; // directory for output files

  // changeable parameters
  std::string GEO_PARAM_FILE; // the name of file where geometry is defined (with a path to geometry directory)
  int CREATE_MESH; // 0 - use old mesh, 1 - create new mesh, 2 - use old mesh, but if some files are absent, create them
  double CL_MAIN; // characteristic length of the mesh in the main domain
  double RHO[N_DOMAINS]; // electrical resistivity for every domain ( [1.0, 1.0] by default )
  int ORDER_3D; // the order of the 3D basis funtions defined of the fine mesh (tetrahedra) ( 1 by default )
  int NUM_X, NUM_Y, NUM_Z; // the number of coarse bricks for every direction ( [1, 1, 1] by default)
  int INCL_PERCENT; // the percent of inclusions that exist in the specimen ( 100 by default)
  double X_BEG, X_END, Y_BEG, Y_END, Z_BEG, Z_END; // the computational region ( [0, 1] x [0, 1] x [0, 1] by default )
  double POTENTIAL_1; // the value of potential on the first plate ( 1.0 by default )
  double POTENTIAL_2; // the value of potential on the second plate ( -1.0 by default )
  int POTENTIAL_TYPE; // 1: left = Potential1, right = Potential2; 2: front = Potential1, back = Potential2; 3: top = Potential1, bottom = Potential2 ( 3 by default )
  int LMAT_COARSE; // integration scheme for local matrices calculation on the coarse grid ( 8 (gauss1) by default )
  int LMAT_FINE; // integration scheme for local matrices calculation on the fine grid (for ORDER_3D > 1 only!) (simplex integration only!) ( 1 (simplex1_1) by default )
  int RHO_INT; // integration scheme for calculation of the effective resistivity ( 8 (gauss1) by default )
  int GAUSS_LMAT; // parameter for Gauss integration for local matrices calculation in every coarse brick ( 1 by default )
  int GAUSS_RHO; // parameter for Gauss integration in electrical resistivity calculation in every coarse brick ( 1 by default )
  double X0_2D, X1_2D, Y0_2D, Y1_2D; // restrict domain for 2D view ( [0, 0] x [0, 0] by default )

  // for MPI
  int *MPI_WORK;

  Param(int argc = 0, char **argv = 0); // constructor

  std::string print(); // return all parameters as a string
};

#endif

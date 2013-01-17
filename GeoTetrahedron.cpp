#include "GeoTetrahedron.h"

GeoTetrahedron::GeoTetrahedron() {
  nVertices = 4;
  vertices = new Node3D[nVertices];
  cpointsFileName = "cpoints/cpointsTetSurf.msh";
  templateName = "TetrahedronTemplate";
}

GeoTetrahedron::~GeoTetrahedron() {
  delete[] vertices;
}

// initialization of the vertices
void GeoTetrahedron::verticesInit() {
  // at first, we create standard tetrahedron with center at origin of coordinates
  double a = lengths[0]; // the length of tetrahedron edge
  double coord[][4] = { { 0.0, 0.0, sqrt(6.0) * a / 4.0, 1.0 }, \
                        { 0.0, sqrt(3.0) * a / 3.0, -sqrt(6.0) * a / 12.0, 1.0 }, \
                        { a / 2.0, -sqrt(3.0) * a / 6.0, -sqrt(6.0) * a / 12.0, 1.0 }, \
                        { -a / 2.0, -sqrt(3.0) * a / 6.0, -sqrt(6.0) * a / 12.0, 1.0 } };
  double coord_new[4];
  for (int i = 0; i < nVertices; i++) {
    toNewCoord(coord[i], coord_new); // from standard tetrahedron at origin to tetrahedron with rotation and translation
    vertices[i].init(coord_new, i);
  }
}

// calculate volume of the tetrahedron
double GeoTetrahedron::volume() {
  double x[4], y[4], z[4]; // Cartesian coordinates of the tetrahedron vertices
  for (int i = 0; i < nVertices; i++) {
    x[i] = vertices[i].getX();
    y[i] = vertices[i].getY();
    z[i] = vertices[i].getZ();
  }
  double detD = Det(x[0], x[1], x[2], x[3],
                    y[0], y[1], y[2], y[3],
                    z[0], z[1], z[2], z[3],
                    1.0, 1.0, 1.0, 1.0); // some key determinant
  double vol = fabs(detD) / 6.0; // tetrahedron volume
  return vol;
}

// check if tetrahedron contains some point with accuracy in tolerance
bool GeoTetrahedron::hasPoint(Node3D *point, double tolerance) {
  double x[4], y[4], z[4]; // Cartesian coordinates of the tetrahedron vertices
  for (int i = 0; i < nVertices; i++) {
    x[i] = vertices[i].getX();
    y[i] = vertices[i].getY();
    z[i] = vertices[i].getZ();
  }
  double px = point->getX(); // Cartesian coordinates of the point under investigation
  double py = point->getY();
  double pz = point->getZ();

  double V = volume(); // the volume of the main tetrahedron
  double v[4]; // the volumes of the 'small' tetrahedra

  double S = Det(px, py, pz, 1.0,
                 x[1], y[1], z[1], 1.0,
                 x[2], y[2], z[2], 1.0,
                 x[3], y[3], z[3], 1.0);
  v[0] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          px, py, pz, 1.0,
          x[2], y[2], z[2], 1.0,
          x[3], y[3], z[3], 1.0);
  v[1] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          x[1], y[1], z[1], 1.0,
          px, py, pz, 1.0,
          x[3], y[3], z[3], 1.0);
  v[2] = fabs(S) / 6.0;
  S = Det(x[0], y[0], z[0], 1.0,
          x[1], y[1], z[1], 1.0,
          x[2], y[2], z[2], 1.0,
          px, py, pz, 1.0);
  v[3] = fabs(S) / 6.0;
  
  if (fabs(V - v[0] - v[1] - v[2] - v[3]) / V < tolerance)
    return true; // tetrahedron contains the point, because the sum of small volumes coicides with the main volume

  return false;
}

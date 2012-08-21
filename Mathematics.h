#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <cstdlib>
#include <vector>

#define PI 3.141592654

#define DISCR_STEP 1000

// the determinant of the third order matrix
// |x1 x2 x3|
// |y1 y2 y3|
// |z1 z2 z3|
inline double Det(double x1, double x2, double x3,
                  double y1, double y2, double y3,
                  double z1, double z2, double z3) {
  return x1 * y2 * z3 + x3 * y1 * z2 + x2 * y3 * z1 -\
         x3 * y2 * z1 - x1 * y3 * z2 - x2 * y1 * z3;
}

// the determinant of the fourth order matrix
// |a1 a2 a3 a4|
// |b1 b2 b3 b4|
// |c1 c2 c3 c4|
// |d1 d2 d3 d4|
inline double Det(double a1, double a2, double a3, double a4,
                  double b1, double b2, double b3, double b4,
                  double c1, double c2, double c3, double c4,
                  double d1, double d2, double d3, double d4) {
  // |a1 a2 a3 a4|        |b2 b3 b4|        |b1 b3 b4|        |b1 b2 b4|        |b1 b2 b3|
  // |b1 b2 b3 b4| = a1 * |c2 c3 c4| - a2 * |c1 c3 c4| + a3 * |c1 c2 c4| - a4 * |c1 c2 c3|
  // |c1 c2 c3 c4|        |d2 d3 d4|        |d1 d3 d4|        |d1 d2 d4|        |d1 d2 d3|
  // |d1 d2 d3 d4|
  double x1 = Det(b2, b3, b4, c2, c3, c4, d2, d3, d4);
  double x2 = Det(b1, b3, b4, c1, c3, c4, d1, d3, d4);
  double x3 = Det(b1, b2, b4, c1, c2, c4, d1, d2, d4);
  double x4 = Det(b1, b2, b3, c1, c2, c3, d1, d2, d3);
  return a1 * x1 - a2 * x2 + a3 * x3 - a4 * x4;
}

// dot product of 2 real vectors
// (a, b) = a1 * b1 + a2 * b2 + ... + an * bn
inline double DotProduct(double a[], double b[], int n) {
  double res = 0.0; // the result
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}

// the norm of vector in the real space
inline double Norm(double a[], int n) {
  return sqrt(DotProduct(a, a, n));
}

// A*b = c,
// where A is a dense square matrix of n-dimension,
// b, c are vectors of n-dimension
inline void MatVec(double *A, double b[], double c[], int n) {
  for (int i = 0; i < n; i++) {
    c[i] = 0.0;
    for (int j = 0; j < n; j++)
      c[i] += A[i * n + j] * b[j];
  }
}

// generate random number between 'begin' and 'end' numbers according to 'margin'
inline double randomBetween(double begin, double end, double margin = 0.0) {
  double value = begin + margin + 1.0 / DISCR_STEP * (rand() % (DISCR_STEP + 1)) * (end - begin - 2.0 * margin);
  return value;
}

// find the position of number in array
inline int findPos(int *array, int dim, int number) {
  for (int i = 0; i < dim; i++)
    if (array[i] == number)
      return i;
  return -1; // if array doesn't contain the number
}

// find the position of number in list
inline int findPos(std::vector<int> list, int number) {
  for (int i = 0; i < (int)list.size(); i++)
    if (list.at(i) == number)
      return i;
  return -1; // if list doesn't contain the number
}

// areatangent
inline double arth(double x) {
  return 0.5 * log((1.0 + x) / (1.0 - x));
}

#endif

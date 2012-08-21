#ifndef CONVERT_H
#define CONVERT_H

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const& s)
    : std::runtime_error(s)
  { }
};

// convert data to string
inline std::string d2s(double x, bool scientific = false, int precision = 6) {
  std::ostringstream o;
  if (scientific) {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }
  if (!(o << x))
    throw BadConversion("d2s(double)");
  return o.str();
}

inline std::string d2s(int x) {
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion("d2s(int)");
  return o.str();
}

#endif

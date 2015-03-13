//==========================================================//
// mixture - Copyright (c) 2012-2013 Mikhail Artemiev       //
//                                                          //
// This library is provided under the terms of MIT license. //
// See the LICENSE file for license information.            //
//==========================================================//

#ifndef CONVERT_H
#define CONVERT_H

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

/**
 * Conversion of integer and float point numbers to string
 */

class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const& s)
    : std::runtime_error(s)
  { }
};

/**
 * Convert data of double-type to string
 */
inline std::string d2s(double x, bool scientific = false, int precision = 6) {
  std::ostringstream o;
  if (scientific) {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }
  if (!(o << x))
    throw BadConversion("Bad conversion from (double) to (string)!");
  return o.str();
}

/**
 * Convert data of integer-type to string
 */
inline std::string d2s(int x) {
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion("Bad conversion from (int) to (string)!");
  return o.str();
}

#endif

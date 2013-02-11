#ifndef REQUIRE_H
#define REQUIRE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

inline void frequire(FILE *f, std::string file = "Unknown file", std::string procedure = "Unknown procedure") {
  if (f == 0) {
    std::cerr << "\n" << procedure << " : File " << file << " cannot be opened!\n" << std::endl;
    exit(1);
  }
}

inline void frequire(std::ifstream &in, std::string file = "Unknown file", std::string procedure = "Unknown procedure") {
  if (!in) {
    std::cerr << "\n" << procedure << " : File " << file << " cannot be opened for reading!\n" << std::endl;
    exit(1);
  }
}

inline void frequire(std::ofstream &out, std::string file = "Unknown file", std::string procedure = "Unknown procedure") {
  if (!out) {
    std::cerr << "\n" << procedure << " : File " << file << " cannot be opened for writing!\n" << std::endl;
    exit(1);
  }
}

inline void require(bool requirement, std::string error = "Requirement failed!", std::string procedure = "Unknown procedure") {
  if (!requirement) {
    std::cerr << "\n" << procedure << " : " << error << "\n" << std::endl;
    exit(1);
  }
}

#endif

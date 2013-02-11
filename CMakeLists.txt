cmake_minimum_required(VERSION 2.8)
project(mixture)

set(mixture_VERSION_MAJOR 2)
set(mixture_VERSION_MINOR 4)
set(mixture_VERSION_PATCH 2)

set(gmsh_BIN gmsh)
set(gmsh_build_cpoints_options "-2 -bin -algo front3d")
set(cpoints_DIR ${PROJECT_SOURCE_DIR}/cpoints)
set(templates_DIR ${PROJECT_SOURCE_DIR})

add_definitions(-g)

configure_file(
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h")

include_directories(${PROJECT_BINARY_DIR})

#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_LINK ON)
#set(Boost_USE_MULTITHREADED ON)
#find_package(Boost COMPONENTS filesystem system)

#if(Boost_FOUND)
#  message ("Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
#  message ("Boost_LIBRARIES = ${Boost_LIBRARIES}")
#  include_directories (${Boost_INCLUDE_DIRS})
#endif(Boost_FOUND)

set(HEADERS
      Convert.h
      Cylinder.h
      Ellipsoid.h
      Geometry.h
      GeoShape.h
      GeoTetrahedron.h
      Mathematics.h
      MeshTriangle3D.h
      Node3D.h
      OrthoBrick.h
      pugiconfig.hpp
      pugixml.hpp
      Require.h
      TriangularMesh.h)

set(SOURCES
      Cylinder.cpp
      Ellipsoid.cpp
      Geometry.cpp
      GeoShape.cpp
      GeoTetrahedron.cpp
      OrthoBrick.cpp
      pugixml.cpp
      TriangularMesh.cpp)

add_library(${PROJECT_NAME} ${HEADERS} ${SOURCES})

#if(Boost_FOUND)
#  target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
#endif(Boost_FOUND)

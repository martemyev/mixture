cmake_minimum_required (VERSION 2.8)

set (PROJECT mixture)

project (${PROJECT})

set (HEADERS
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
       Param.h
       pugiconfig.hpp
       pugixml.hpp
       Require.h
       TriangularMesh.h)

set (SOURCES
       Cylinder.cpp
       Ellipsoid.cpp
       Geometry.cpp
       GeoShape.cpp
       GeoTetrahedron.cpp
       OrthoBrick.cpp
       pugixml.cpp
       TriangularMesh.cpp)

add_library (${PROJECT} ${HEADERS} ${SOURCES})


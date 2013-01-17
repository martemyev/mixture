set (HEADERS_DIR "${CMAKE_SOURCE_DIR}/headers")

message ("HEADERS_DIR = ${HEADERS_DIR}")
include_directories (${HEADERS_DIR})

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

set (PROJECT_NAME mixture)
project (${PROJECT_NAME})

add_library (${PROJECT_NAME} ${SOURCES})



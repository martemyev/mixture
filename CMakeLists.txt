cmake_minimum_required(VERSION 2.8)

# --- compiler options ---
# taken from Gmsh's CMakeLists.txt:
# if CMAKE_BUILD_TYPE is specified use it; otherwise set the default
# build type to "RelWithDebInfo" ("-O2 -g" with gcc) prior to calling
# project()
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif(DEFINED CMAKE_BUILD_TYPE)
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# --- project ---
project(mixture CXX)

# --- set version ---
set(mixture_VERSION_MAJOR 2)
set(mixture_VERSION_MINOR 4)
set(mixture_VERSION_PATCH 3)

# --- necessary paths and options ---
set(gmsh_BIN "gmsh" CACHE PATH "Path to Gmsh binary file")
set(gmsh_build_cpoints_options "-2 -bin -algo front3d" CACHE STRING "Gmsh meshing options")
set(cpoints_DIR ${PROJECT_SOURCE_DIR}/cpoints CACHE PATH "Path to directory with cpoints (control points) files")
set(templates_DIR ${PROJECT_SOURCE_DIR} CACHE PATH "Path to directory where the templates.geo file is")

#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_LINK ON)
#set(Boost_USE_MULTITHREADED ON)
#find_package(Boost COMPONENTS filesystem system)

#if(Boost_FOUND)
#  message ("Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
#  message ("Boost_LIBRARIES = ${Boost_LIBRARIES}")
#  include_directories (${Boost_INCLUDE_DIRS})
#endif(Boost_FOUND)

# --- sources ---
set(SOURCES_DIR ${PROJECT_SOURCE_DIR}/sources)
aux_source_directory(${SOURCES_DIR} SOURCES_LIST)
#message("SOURCES_LIST = ${SOURCES_LIST}")

# --- headers ---
set(HEADERS_DIR ${PROJECT_SOURCE_DIR}/headers)
include_directories(${HEADERS_DIR})
file(GLOB HEADERS_LIST "${HEADERS_DIR}/*.h" "${HEADERS_DIR}/*.hpp")
#message("HEADERS_LIST = ${HEADERS_LIST}")

# --- write the variables initialized here to Config.h file ---
configure_file(
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h")
include_directories(${PROJECT_BINARY_DIR}) # because Config.h file is in binary directory we need to include it

add_library(${PROJECT_NAME} ${HEADERS_LIST} ${SOURCES_LIST})

#if(Boost_FOUND)
#  target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
#endif(Boost_FOUND)

# --- packaging of the library ---
# -- make package --
set(CPACK_PACKAGE_VERSION_MAJOR "${mixture_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${mixture_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${mixture_VERSION_PATCH}")
set(CPACK_GENERATOR "TGZ")
# -- make package_source --
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "./*.user;build;doc;cpoints;.git;.gitignore")
include(CPack)

# --- autodocumentation using doxygen ---
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(
    "${CMAKE_SOURCE_DIR}/Doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
  add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
else(DOXYGEN_FOUND)
  message("WARNING: Doxygen is not found - Autodocumentation will not be created")
endif(DOXYGEN_FOUND)


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
set(mixture_VERSION_PATCH 5)

# --- necessary paths and options ---
if(WIN32)
  set(gmsh_BIN "C:/Program Files/gmsh_directory/gmsh.exe" CACHE PATH "Path to Gmsh binary file")
  set(WIN32_OS 1)
else(WIN32)
  set(gmsh_BIN "gmsh" CACHE PATH "Path to Gmsh binary file")
  set(WIN32_OS 0)
endif(WIN32)
set(gmsh_build_cpoints_options "-2 -bin -algo front3d" CACHE STRING "Gmsh meshing options")
set(cpoints_DIR ${PROJECT_SOURCE_DIR}/cpoints CACHE PATH "Path to directory with cpoints (control points) files")
set(templates_DIR ${PROJECT_SOURCE_DIR} CACHE PATH "Path to directory where the templates.geo file is")

set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_LINK ON)
#set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS filesystem system)
if(Boost_FOUND)
  message("Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
  message("Boost_LIBRARIES = ${Boost_LIBRARIES}")
  include_directories(${Boost_INCLUDE_DIRS})
  set(BOOST_FS 1)
else(Boost_FOUND)
  if(WIN32)
    set(boost_INCL_DIR "C:/Users/artemiev/Download/boost_1_53_0/boost_1_53_0" CACHE PATH "The path to the boost headers")
    set(boost_LIB_DIR "C:/Users/artemiev/Download/boost_1_53_0/boost_1_53_0/stage/lib" CACHE PATH "The path to the boost libraries")
  else(WIN32)
    set(boost_INCL_DIR "/home/artemiev/install/boost_1_53_0" CACHE PATH "The path to the boost headers")
    set(boost_LIB_DIR "/home/artemiev/install/boost_1_53_0/stage/lib" CACHE PATH "The path to the boost libraries")
  endif(WIN32)
  include_directories(${boost_INCL_DIR})
  find_library(BOOST_FILESYSTEM boost_filesystem libboost_filesystem ${boost_LIB_DIR})
  find_library(BOOST_SYSTEM boost_system libboost_system ${boost_LIB_DIR})
  message("BOOST_FILESYSTEM = ${BOOST_FILESYSTEM}")
  message("BOOST_SYSTEM = ${BOOST_SYSTEM}")
  if(BOOST_FILESYSTEM)
    set(BOOST_FS 1)
  else(BOOST_FILESYSTEM)
    set(BOOST_FS 0)
  endif(BOOST_FILESYSTEM)
endif(Boost_FOUND)

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

# --- check that Gmsh works ---
message("check that Gmsh works: its version is")
execute_process(COMMAND ${gmsh_BIN} --version RESULT_VARIABLE gmsh_CHECKING)
message("gmsh_CHECKING = ${gmsh_CHECKING}")
if(${gmsh_CHECKING} EQUAL 0)
  message("Gmsh binary file is working properly!")
else(${gmsh_CHECKING} EQUAL 0)
  #message(FATAL_ERROR "Gmsh binary file is not working!")
  message("")
  message("PLEASE PAY ATTENTION!")
  message("Gmsh binary file is not working!")
  message("You can miss this warning only if you're strongly sure that Gmsh binary file is working!")
endif(${gmsh_CHECKING} EQUAL 0)


# --- we want to build static library ---
set(BUILD_SHARED_LIBS OFF)

add_library(${PROJECT_NAME} ${HEADERS_LIST} ${SOURCES_LIST})
if(Boost_FOUND)
  target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
else(Boost_FOUND)
  if(BOOST_FS)
    target_link_libraries(${PROJECT_NAME} ${BOOST_FILESYSTEM} ${BOOST_SYSTEM})
  endif(BOOST_FS)
endif(Boost_FOUND)

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

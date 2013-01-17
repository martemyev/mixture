cmake_minimum_required (VERSION 2.8)

set (SOURCES_DIR "${CMAKE_SOURCE_DIR}/sources")
add_subdirectory (${SOURCES_DIR})

set (HEADERS_DIR "${CMAKE_SOURCE_DIR}/headers")
add_subdirectory (${HEADERS_DIR})


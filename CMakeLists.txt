# Specify cmake version
cmake_minimum_required(VERSION 3.2)
project(sympler)


# Default location of executables in top build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# Specify libraries to link executables
set(link_lib
-lgsl
-lgslcblas
-ldl
-lxml2
-lm
)


# Find required and optional libs

# XML2 is required
find_package(LibXml2 REQUIRED)

# GSL is required
find_package(GSL REQUIRED)		 				

# CPPUNIT is required
find_path(CPPUNIT_INCLUDE_DIR cppunit/Test.h)
find_library(CPPUNIT_LIBRARY NAMES cppunit)
if(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  message(STATUS "Found CppUnit: ${CPPUNIT_LIBRARY}")
else(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  message(FATAL_ERROR "Could not find CppUnit")
endif(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)

# SDL is optional
find_package(SDL)
if(SDL_FOUND)
  add_definitions(-DHAVE_SDL)
endif(SDL_FOUND)

# FREESTEAM is optional
execute_process(COMMAND freesteam-config --version OUTPUT_VARIABLE freesteam_version OUTPUT_STRIP_TRAILING_WHITESPACE)
if(${freesteam_version} GREATER 1.0)
  message("-- Found Freesteam (found version ${freesteam_version})")
  add_definitions(-DHAVE_FREESTEAM)
  set(link_lib ${link_lib} -lfreesteam)
else(${freesteam_version} GREATER 1.0)
  message("-- Could NOT find FREESTEAM")
endif(${freesteam_version} GREATER 1.0)


# Test Endian System
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DIS_BIG_ENDIAN)
endif()


# Add subdirectories to the build
add_subdirectory(source/src)
add_subdirectory(source/src_test)


#add_definitions(-DHAVE_JAMA_JAMA_LU_H)

#add_definitions(-DHAVE_TNT_TNT_H)

#add_definitions(-DWITH_ARRAY_TYPES)		???
#add_definitions(-D_MPI)			???
#......


# target_link_libraries(sympler GSL::gsl ${CMAKE_DL_LIBS} ${LIBXML2_LIBRARIES})		# without freesteam
# link SDL, JAMA, TNT, ... if set!!!




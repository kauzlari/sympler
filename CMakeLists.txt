# Specify cmake version
cmake_minimum_required (VERSION 3.2)
project(sympler)


# Default location of executables in top build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# Add subdirectories to the build.
add_subdirectory(source/src)
add_subdirectory(source/src_test)


# Find required and optional libs
find_package(LibXml2 REQUIRED)					# sympler definitely needs XML2

find_package(GSL REQUIRED)		 				# sympler definitely needs GSL

FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/Test.h)	# sympler definitely needs CPPUNIT
FIND_LIBRARY(CPPUNIT_LIBRARY NAMES cppunit)
IF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  MESSAGE(STATUS "Found CppUnit: ${CPPUNIT_LIBRARY}")
ELSE (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  MESSAGE(FATAL_ERROR "Could not find CppUnit")
ENDIF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)

find_package(SDL)								# looking for optional SDL
if(SDL_FOUND)
  add_definitions(-DHAVE_SDL)
endif(SDL_FOUND)

execute_process(COMMAND freesteam-config --version OUTPUT_VARIABLE freesteam_version)	# looking for optional FREESTEAM
if(${freesteam_version} GREATER 1.0)
  add_definitions(-DHAVE_FREESTEAM)
else(${freesteam_version} GREATER 1.0)
  message("-- Could NOT find FREESTEAM")
endif(${freesteam_version} GREATER 1.0)

include (TestBigEndian)							# Test Endian System
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DIS_BIG_ENDIAN)
endif()


#add_definitions(-DHAVE_JAMA_JAMA_LU_H)

#add_definitions(-DHAVE_TNT_TNT_H)

#add_definitions(-DWITH_ARRAY_TYPES)		???
#add_definitions(-D_MPI)			???
#......


# target_link_libraries(sympler GSL::gsl ${CMAKE_DL_LIBS} ${LIBXML2_LIBRARIES})		# without freesteam
# link SDL, JAMA, TNT, ... if set!!!




# Specify cmake version
cmake_minimum_required(VERSION 3.2)
project(sympler)


# Default location of executables in top build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# Include header files
include_directories(${PROJECT_SOURCE_DIR}/source/include/basic)
include_directories(${PROJECT_SOURCE_DIR}/source/include/boundary)
include_directories(${PROJECT_SOURCE_DIR}/source/include/calculator)
include_directories(${PROJECT_SOURCE_DIR}/source/include/force)
include_directories(${PROJECT_SOURCE_DIR}/source/include/function_parser)
include_directories(${PROJECT_SOURCE_DIR}/source/include/geometry)
include_directories(${PROJECT_SOURCE_DIR}/source/include/geometry/stl)
include_directories(${PROJECT_SOURCE_DIR}/source/include/geometry/sat)
include_directories(${PROJECT_SOURCE_DIR}/source/include/integrator)
include_directories(${PROJECT_SOURCE_DIR}/source/include/meter)
include_directories(${PROJECT_SOURCE_DIR}/source/include/meter/grid_meter)
include_directories(${PROJECT_SOURCE_DIR}/source/include/particle_creator)
include_directories(${PROJECT_SOURCE_DIR}/source/include/postprocessor)
include_directories(${PROJECT_SOURCE_DIR}/source/include/reflector)
include_directories(${PROJECT_SOURCE_DIR}/source/include/callable)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol/freesteam)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol/particle_cache)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol/val_calculator_pair)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol/val_calculator_part)
include_directories(${PROJECT_SOURCE_DIR}/source/include/symbol/other)
include_directories(${PROJECT_SOURCE_DIR}/source/include/tools)
include_directories(${PROJECT_SOURCE_DIR}/source/include/weighting_function)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/basic)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/boundary)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/calculator)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/force)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/function_parser)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/geometry)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/geometry/stl)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/geometry/sat)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/integrator)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/meter)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/meter/grid_meter)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/particle_creator)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/postprocessor)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/reflector)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/callable)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol/freesteam)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol/particle_cache)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol/val_calculator_pair)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol/val_calculator_part)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/symbol/other)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/tools)
include_directories(${PROJECT_SOURCE_DIR}/source/include_test/weighting_function)
include_directories(${SDKROOT}/usr/include/libxml2)
include_directories(/usr/local/lib)
include_directories(${LIBXML2_INCLUDE_DIR})


# Specify required libraries to link executables
set(link_lib
-lgsl
-lgslcblas
-ldl
-lxml2
-lm
-pthread
-L/usr/local/lib
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

# SDL is linked, if installed
find_package(SDL)
if(SDL_FOUND)
  if(${static})
    message("-- Dont build SDL, because of statically linking")
  else(${static})
    add_definitions(-DHAVE_SDL)
    set(link_lib ${link_lib} -lSDL)
  endif(${static})
endif(SDL_FOUND)

# FREESTEAM is linked, if installed
execute_process(COMMAND freesteam-config --version OUTPUT_VARIABLE freesteam_version OUTPUT_STRIP_TRAILING_WHITESPACE)
if(${freesteam_version})
  message("-- Found Freesteam (found version ${freesteam_version})")
  if(${static})
    message("-- Dont build Freesteam, because of statically linking")
  else(${static})
    add_definitions(-DHAVE_FREESTEAM)
    set(link_lib ${link_lib} -lfreesteam)
  endif(${static})
else(${freesteam_version})
  message("-- Could NOT find FREESTEAM")
endif(${freesteam_version})

# TNT & JAMA are optional, if -Dtnt=1 is set
if(${tnt})
  message("-- Use TNT & JAMA Libraries")
  add_definitions(-DHAVE_TNT_TNT_H)
  add_definitions(-DHAVE_JAMA_JAMA_LU_H)
endif(${tnt})

# Compile IntegratorStaticLSE using SuperLU (+ TNT & JAMA) Library, if -DIntegratorStaticLSE=1 is set
if(${IntegratorStaticLSE})
  message("-- Use SuperLU Libraries to compile IntegratorStaticLSE")
  add_definitions(-DHAVE_TNT_TNT_H)
  add_definitions(-DHAVE_JAMA_JAMA_LU_H)
  add_definitions(-DWITH_INTEGRATORSTATICLSE)
  set(link_lib ${link_lib} -lsuperlu)
endif(${IntegratorStaticLSE})


# Create static (standalone) executable, if -Dstatic=1 is set
if(${static})
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")
  message("-- Build static (standalone) binaries!")
  set(link_lib ${link_lib} -lz -llzma -licuuc -licudata)
endif(${static})


# Test Endian System
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DIS_BIG_ENDIAN)
endif()


# Add subdirectories to the build
add_subdirectory(source/src)
add_subdirectory(source/src_test)



# Specify cmake version
cmake_minimum_required(VERSION 2.8)
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
if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.1)
  find_package(GSL REQUIRED)
else(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.1)
  execute_process(COMMAND gsl-config --version OUTPUT_VARIABLE gsl_version OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(${gsl_version})
    message(STATUS "Found GSL (found version ${gsl_version})")
  else(${gsl_version})
    message(FATAL_ERROR "Could NOT find GSL")
  endif(${gsl_version})
endif(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.1)

# CPPUNIT is optional, Unittest compiling only if CPPUNIT Library existing
find_path(CPPUNIT_INCLUDE_DIR cppunit/Test.h)
find_library(CPPUNIT_LIBRARY NAMES cppunit)
if(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  message(STATUS "Found CppUnit: ${CPPUNIT_LIBRARY}")
  set(CPPUNIT_FOUND TRUE)
else(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
  message(STATUS "Could NOT find CppUnit: No Unittest compiling")
  set(CPPUNIT_FOUND FALSE)
endif(CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)

# SDL is linked, if installed (SDL not linked, if on MacOS or static executable is required)
if(NOT APPLE)
  find_package(SDL)
  if(SDL_FOUND)
    if(${static})
      message(STATUS "Dont build SDL, because of statically linking")
    else(${static})
      add_definitions(-DHAVE_SDL)
      set(link_lib ${link_lib} -lSDL)
    endif(${static})
  endif(SDL_FOUND)
else(NOT APPLE)
  message(STATUS "Dont check SDL library for MacOS")
endif(NOT APPLE)

# FREESTEAM is linked, if installed
execute_process(COMMAND freesteam-config --version OUTPUT_VARIABLE freesteam_version OUTPUT_STRIP_TRAILING_WHITESPACE)
if(${freesteam_version})
  message(STATUS "Found Freesteam (found version ${freesteam_version})")
  if(${static})
    message(STATUS "Dont build Freesteam, because of statically linking")
  else(${static})
    add_definitions(-DHAVE_FREESTEAM)
    set(link_lib ${link_lib} -lfreesteam)
  endif(${static})
else(${freesteam_version})
  message(STATUS "Could NOT find Freesteam")
endif(${freesteam_version})

# TNT & JAMA are optional, if -Dtnt=1 is set
if(${tnt})
  message(STATUS "Use TNT & JAMA Libraries")
  add_definitions(-DHAVE_TNT_TNT_H)
  add_definitions(-DHAVE_JAMA_JAMA_LU_H)
endif(${tnt})

# Compile IntegratorStaticLSE using SuperLU (+ TNT & JAMA) Library, if -DIntegratorStaticLSE=1 is set
if(${IntegratorStaticLSE})
  message(STATUS "Use SuperLU (+ TNT & JAMA) Libraries to compile IntegratorStaticLSE")
  add_definitions(-DHAVE_TNT_TNT_H)
  add_definitions(-DHAVE_JAMA_JAMA_LU_H)
  add_definitions(-DWITH_INTEGRATORSTATICLSE)
  set(link_lib ${link_lib} -lsuperlu)
endif(${IntegratorStaticLSE})


# Create static (standalone) executable, if -Dstatic=1 is set
if(${static})
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")
  message(STATUS "Build static (standalone) binaries!")
  set(link_lib ${link_lib} -lz -llzma -licuuc -licudata)
endif(${static})


# Test Endian System
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DIS_BIG_ENDIAN)
endif(IS_BIG_ENDIAN)


# Version control tracking
message(STATUS "Execute version control tracking")
# Check if git repository
execute_process(COMMAND git rev-parse --is-inside-work-tree WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} ERROR_QUIET OUTPUT_VARIABLE git_repository OUTPUT_STRIP_TRAILING_WHITESPACE)
if(${git_repository})
  # Get version number from last commit date by using git log history
  execute_process(COMMAND git log -1 --format=%ci WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE last_commit_date OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Create version number string
  string(SUBSTRING ${last_commit_date} 0 16 version_date)
  string(REPLACE ":" "-" version_date ${version_date})
  string(REPLACE " " "_" version_date ${version_date})
  # Get branch name using git branch
  execute_process(COMMAND git rev-parse --abbrev-ref HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE branch_name OUTPUT_STRIP_TRAILING_WHITESPACE)
else(${git_repository})
  if(APPLE)
    # Get version number from last changed date of file for macOS systems
    set(get_date_info_command stat -f %Sm -t %Y-%m-%d CMakeLists.txt)
    execute_process(COMMAND ${get_date_info_command} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE last_commit_date OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Create version number string
    string(SUBSTRING ${last_commit_date} 0 10 version_date)
  else(APPLE)
    # Get version number from last changed date of file in Linux systems
    set(get_date_info_command stat -c %y CMakeLists.txt)
    execute_process(COMMAND ${get_date_info_command} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE last_commit_date OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Create version number string
    string(SUBSTRING ${last_commit_date} 0 16 version_date)
    string(REPLACE ":" "-" version_date ${version_date})
    string(REPLACE " " "_" version_date ${version_date})
    # Get branch name using source directory name
    get_filename_component(directory_name ${CMAKE_CURRENT_LIST_DIR} NAME)
    string(FIND ${directory_name} "-" pos REVERSE)
    MATH(EXPR pos "${pos}+1")
    if(pos)
      string(SUBSTRING ${directory_name} ${pos} -1 branch_name)
    endif(pos)
  endif(APPLE)
endif(${git_repository})
# Concatenate branch name if it is not master
if(branch_name AND NOT(${branch_name} MATCHES "master"))
  string(CONCAT version_date ${version_date}_${branch_name})
endif(branch_name AND NOT(${branch_name} MATCHES "master"))
# Set macro with version number
add_definitions(-D_VERSIONNUMBER="v${version_date}")


# Add subdirectories to the build
add_subdirectory(source/src)
if(${CPPUNIT_FOUND})
  add_subdirectory(source/src_test)
endif(${CPPUNIT_FOUND})



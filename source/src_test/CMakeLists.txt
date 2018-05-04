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
include_directories(${LIBXML2_INCLUDE_DIR})


# Define source files to be compiled
set(basic_test 
basic/smart_enumTest.cpp
)

set(boundary_test
)

set(calculator_test
)

set(callable_test
)

set(force_test
)

set(function_parser_test
)

set(geometry_test 
geometry/vertex_listTest.cpp
)

set(integrator_test
integrator/integrator_velocity_verletTest.cpp
)

set(meter_test
)

set(particle_creator_test
)

set(postprocessor_test
)

set(reflector_test
)

set(symbol_test
symbol/freesteam/pca_iapws-if97_1varTest.cpp
symbol/freesteam/pca_iapws-if97_2varTest.cpp
symbol/freesteam/pca_iapws-if97_cpTest.cpp
symbol/freesteam/pca_iapws-if97_etaTest.cpp
symbol/freesteam/pca_iapws-if97_kappaTest.cpp
symbol/freesteam/pca_iapws-if97_pTest.cpp
symbol/freesteam/pca_iapws-if97_rhoTest.cpp
symbol/freesteam/pca_iapws-if97Test.cpp
)

set(tools_test
)

set(weighting_function_test
)


# Create executable symplerTest
add_executable(symplerTest 
mainTest.cpp
${basic_test}
${geometry_test}
${integrator_test}
${symbol_test}
)

# Link libraries to executable symplerTest
target_link_libraries(symplerTest
basic
boundary
calculator
force
function_parser
geometry
integrator
meter
particle_creator
postprocessor
reflector
callable
symbol
weighting_function
cppunit
-lgsl
-lgslcblas
-lfreesteam
${CMAKE_DL_LIBS} 
${LIBXML2_LIBRARIES}
-lxml2
-lm
)


# Execute tests after building symplerTest
enable_testing()
ADD_CUSTOM_COMMAND(TARGET symplerTest POST_BUILD COMMAND ${CMAKE_BINARY_DIR}/symplerTest)



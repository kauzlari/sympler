# Define source files to be compiled
set(basic_test 
basic/smart_enumTest.cpp
)

set(boundary_test
)

set(calculator_test
)

set(callable_test
callable/cbl_pair_particle_arbitraryTest.cpp
callable/cbl_pair_particle_tensorTest.cpp
callable/cbl_pair_particle_vectorTest.cpp
)

set(force_test
)

set(function_parser_test
)

set(geometry_test 
geometry/vertex_listTest.cpp
)

set(integrator_test
integrator/integrator_positionTest.cpp
integrator/integrator_pos_vel_step2Test.cpp
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
symbol/particle_cache/fake_function_particle.cpp
symbol/particle_cache/fake_particle_cache.cpp
symbol/particle_cache/fake_pca_2nd_sph_deriv_corr.cpp
symbol/particle_cache/fake_symbol_f_particle_scalar.cpp
symbol/particle_cache/particle_velsTest.cpp
symbol/particle_cache/particle_cacheTest.cpp
symbol/particle_cache/pca_2nd_sph_deriv_corrTest.cpp
symbol/particle_cache/symbol_f_particle_arbitraryTest.cpp
symbol/particle_cache/symbol_f_particle_scalarTest.cpp
symbol/particle_cache/symbol_f_particle_velsTest.cpp
)

set(tools_test
)

set(weighting_function_test
)


# Create executable symplerTest
add_executable(symplerTest 
mainTest.cpp
${basic_test}
${callable_test}
${geometry_test}
${integrator_test}
${symbol_test}
)

# The following list is complete (2018-05-25). Add directories to upper SUBDIRS as soon as you have at least one compilable *Test.cpp file
#add_executable(symplerTest 
#mainTest.cpp
#${basic_test}
#${boundary_test}
#${calculator_test}
#${callable_test}
#${force_test}
#${function_parser_test}
#${geometry_test}
#${integrator_test}
#${meter_test}
#${particle_creator_test}
#${postprocessor_test}
#${reflector_test}
#${symbol_test}
#${tools_test}
#${weighting_function_test}
#)


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
${link_lib}
)


# Execute tests after build
enable_testing()
#ADD_CUSTOM_COMMAND(TARGET symplerTest POST_BUILD COMMAND ${CMAKE_BINARY_DIR}/symplerTest)
ADD_CUSTOM_COMMAND(
	TARGET symplerTest 
	POST_BUILD 
		COMMAND 
			cp ${CMAKE_BINARY_DIR}/symplerTest 
					${CMAKE_BINARY_DIR}/symplerTestBackup 
		COMMAND ${CMAKE_BINARY_DIR}/symplerTest)



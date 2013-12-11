
noinst_LTLIBRARIES = libintegrator.la

libintegrator_la_SOURCES = \
	integrator_energy.cpp \
	integrator_lse.cpp \
	integrator_omelyan.cpp \
	integrator_position.cpp \
	integrator_scalar.cpp \
	integrator_scalar_lambda.cpp \
	integrator_static.cpp \
	integrator_static_lse.cpp \
	integrator_tensor.cpp \
	integrator_velocity_verlet_disp.cpp \
	integrator_velocity_verlet_disp_x.cpp \
	integrator_tensor_lambda.cpp \
	integrator_vector.cpp \
	integrator_vector_lambda.cpp \
	integrator_velocity_verlet.cpp \
	integrator_velocity_verlet_pressure.cpp \
	integrator_velocity_verlet_symbols.cpp \
	integrator_velocity_verlet_x.cpp \
	integrator_velocity_verlet_x_full.cpp

libintegrator_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head


noinst_LTLIBRARIES = libboundary.la

libboundary_la_SOURCES = \
	boundary_arbitrary.cpp \
	boundary_bubble_jet.cpp \
	boundary_cuboid.cpp \
	boundary_couette.cpp \
	boundary_diffusor.cpp \
	boundary_step.cpp \
	boundary_stl.cpp \
	boundary_stl_periodic.cpp\
	boundary_obstacle.cpp \
	boundary_with_inlet.cpp
	

libboundary_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/boundary/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/geometry/stl/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/particle_creator/head



noinst_LTLIBRARIES = libgeometry.la

libgeometry_la_SOURCES = \
	geometric_algorithms.cpp \
	geometric_primitives.cpp \
	vertex_list.cpp \
	wall.cpp \
	wall_triangle.cpp \
	wall_container.cpp

libgeometry_la_CXXFLAGS = @GSL_CFLAGS@ @XML_CFLAGS@

libgeometry_la_LDFLAGS = @LT_LDFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/calculator/head

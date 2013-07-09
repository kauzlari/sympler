
noinst_LTLIBRARIES = libstl.la

libstl_la_SOURCES = \
	stl_ascii.cpp \
	stl_binary.cpp

libstl_la_CXXFLAGS = @GSL_CFLAGS@ @XML_CFLAGS@

libstl_la_LDFLAGS = @LT_LDFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/calculator/head -I$(top_srcdir)/src/geometry/stl/head -I$(top_srcdir)/src/integrator/head

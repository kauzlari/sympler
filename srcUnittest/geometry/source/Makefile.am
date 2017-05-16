
noinst_LTLIBRARIES = libgeometryTest.la

libgeometryTest_la_SOURCES = \
	vertex_listTest.cpp

libgeometryTest_la_CXXFLAGS = @GSL_CFLAGS@ @XML_CFLAGS@

libgeometryTest_la_LDFLAGS = @LT_LDFLAGS@

INCLUDES = -I$(top_srcdir)/srcUnittest/geometry/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/calculator/head

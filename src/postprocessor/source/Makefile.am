
noinst_LTLIBRARIES = libpostprocessor.la

libpostprocessor_la_SOURCES = \
	output.cpp \
	output_file.cpp \
	output_mathematica.cpp \
	output_dcd.cpp \
	output_pdb.cpp \
	output_vtk.cpp

libpostprocessor_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/postprocessor/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/geometry/head

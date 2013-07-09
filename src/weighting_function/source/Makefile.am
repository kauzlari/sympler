
noinst_LTLIBRARIES = libweighting_function.la

libweighting_function_la_SOURCES = \
	wf_input.cpp \
	wf_linear.cpp \
	wf_lucy.cpp \
	wf_lucy_0th_order.cpp \
	wf_lucy_1st_order.cpp \
	wf_square.cpp

libweighting_function_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/weighting_function/head

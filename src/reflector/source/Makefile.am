
noinst_LTLIBRARIES = libreflector.la

libreflector_la_SOURCES = \
	reflector_bounce_back.cpp \
	reflector_mirror.cpp \
	reflector_stochastic.cpp \
	reflector_thermalize.cpp \
	reflector_thermalize_ie.cpp \
	reflector_with_rng.cpp

libreflector_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/reflector/head

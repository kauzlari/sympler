bin_PROGRAMS = sympler

sympler_SOURCES = main.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/boundary/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/force/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/meter/grid_meter/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/postprocessor/head -I$(top_srcdir)/src/reflector/head -I$(top_srcdir)/src/particle_creator/head -I$(top_srcdir)/src/thermostat/head -I$(top_srcdir)/src/weighting_function/head

sympler_LDADD = \
	-Lbasic/source \
	-Lboundary/source \
	-Lcalculator/source \
	-Lforce/source \
	-Lfunction_parser/source \
	-Lgeometry/source \
	-Lgeometry/stl/source \
	-Lintegrator/source \
	-Lmeter/source \
	-Lmeter/grid_meter/source \
	-Lparticle_creator/source \
	-Lpostprocessor/source \
	-Lreflector/source \
	-Lthermostat/source \
	-Lweighting_function/source \
	\
	-lboundary \
	-lforce \
	-lgeometry \
	-lstl \
	-lintegrator \
	-lmeter \
	-lgrid_meter \
	-lparticle_creator \
	-lpostprocessor \
	-lreflector \
	-lthermostat \
	-lweighting_function \
	-lbasic \
	-lcalculator \
	-lfunction_parser \
	\
	-ldl \
	\
        @GSL_LIBS@ @SDL_LIBS@ @XML_LIBS@ \
	-llzma -licui18n -licuuc -licudata
#	-lz -lpthread -llzma -licui18n -licuuc -licudata
# -lz seems to be automatically added by XML_LIBS
# -lpthread was replaced by -pthread in sympler_LDFLAGS below
# -llzma...-licudata are required for but not added by XML_LIBS when linking with --enable-all-static. In fact this is the output of 'icu-config --ldflags', but currently no idea how to include this in a more automated way


sympler_LDFLAGS = \
	@LT_LDFLAGS@ -pthread

# from good old subversion times
#sympler_CXXFLAGS = \
#	@GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@ -D'SVNREV="$(SVNREV)"'

sympler_CXXFLAGS = \
	@GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

sympler_DEPENDENCIES = \
	basic/source/libbasic.la \
	boundary/source/libboundary.la \
	calculator/source/libcalculator.la \
	force/source/libforce.la \
	function_parser/source/libfunction_parser.la \
	geometry/source/libgeometry.la \
	geometry/stl/source/libstl.la \
	integrator/source/libintegrator.la \
	meter/source/libmeter.la \
	meter/grid_meter/source/libgrid_meter.la \
	particle_creator/source/libparticle_creator.la \
	postprocessor/source/libpostprocessor.la \
	reflector/source/libreflector.la \
	thermostat/source/libthermostat.la \
	weighting_function/source/libweighting_function.la

SUBDIRS = \
	basic \
	boundary \
	calculator \
	force \
	function_parser \
	geometry \
	integrator \
	meter \
	particle_creator \
	postprocessor \
	reflector \
	thermostat \
	tools \
	weighting_function


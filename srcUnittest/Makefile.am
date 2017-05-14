bin_PROGRAMS = symplerTest

all:
	./symplerTest

symplerTest_SOURCES = mainTest.cpp

INCLUDES = -I$(top_srcdir)/srcUnittest/basic/head -I$(top_srcdir)/srcUnittest/boundary/head -I$(top_srcdir)/srcUnittest/calculator/head -I$(top_srcdir)/srcUnittest/force/head -I$(top_srcdir)/srcUnittest/function_parser/head -I$(top_srcdir)/srcUnittest/geometry/head -I$(top_srcdir)/srcUnittest/meter/head -I$(top_srcdir)/srcUnittest/meter/grid_meter/head -I$(top_srcdir)/srcUnittest/integrator/head -I$(top_srcdir)/srcUnittest/postprocessor/head -I$(top_srcdir)/srcUnittest/reflector/head -I$(top_srcdir)/srcUnittest/particle_creator/head -I$(top_srcdir)/srcUnittest/thermostat/head -I$(top_srcdir)/srcUnittest/weighting_function/head -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/boundary/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/force/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/meter/grid_meter/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/postprocessor/head -I$(top_srcdir)/src/reflector/head -I$(top_srcdir)/src/particle_creator/head -I$(top_srcdir)/src/thermostat/head -I$(top_srcdir)/src/weighting_function/head

symplerTest_LDADD = \
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
	-ldl \
	\
	-lstdc++ \
	\
	-lcppunit \
	\
	-lm \
	\
        @GSL_LIBS@ @SDL_LIBS@ @XML_LIBS@ \
	-llzma -licui18n -licuuc -licudata
#	-lz -lpthread -llzma -licui18n -licuuc -licudata
# -lz seems to be automatically added by XML_LIBS
# -lpthread was replaced by -pthread in sympler_LDFLAGS below
# -llzma...-licudata are required for but not added by XML_LIBS when linking with --enable-all-static. In fact this is the output of 'icu-config --ldflags', but currently no idea how to include this in a more automated way


symplerTest_LDFLAGS = \
	@LT_LDFLAGS@ -pthread \
	\
	../src/basic/source/libbasic.la \
	../src/boundary/source/libboundary.la \
	../src/calculator/source/libcalculator.la \
	../src/force/source/libforce.la \
	../src/function_parser/source/libfunction_parser.la \
	../src/geometry/source/libgeometry.la \
	../src/geometry/stl/source/libstl.la \
	../src/integrator/source/libintegrator.la \
	../src/meter/source/libmeter.la \
	../src/meter/grid_meter/source/libgrid_meter.la \
	../src/particle_creator/source/libparticle_creator.la \
	../src/postprocessor/source/libpostprocessor.la \
	../src/reflector/source/libreflector.la \
	../src/thermostat/source/libthermostat.la \
	../src/weighting_function/source/libweighting_function.la \
	\
	geometry/source/libgeometryTest_la-vertex_listTest.o

# from good old subversion times
#sympler_CXXFLAGS = \
#	@GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@ -D'SVNREV="$(SVNREV)"'

symplerTest_CXXFLAGS = \
	@GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

symplerTest_DEPENDENCIES = \
	mainTest.cpp \
	\
	basic/source/libbasicTest.la \
	boundary/source/libboundaryTest.la \
	calculator/source/libcalculatorTest.la \
	force/source/libforceTest.la \
	function_parser/source/libfunction_parserTest.la \
	geometry/source/libgeometryTest.la \
	geometry/stl/source/libstlTest.la \
	integrator/source/libintegratorTest.la \
	meter/source/libmeterTest.la \
	meter/grid_meter/source/libgrid_meterTest.la \
	particle_creator/source/libparticle_creatorTest.la \
	postprocessor/source/libpostprocessorTest.la \
	reflector/source/libreflectorTest.la \
	thermostat/source/libthermostatTest.la \
	weighting_function/source/libweighting_functionTest.la \
	\
	../src/basic/source/libbasic.la \
	../src/boundary/source/libboundary.la \
	../src/calculator/source/libcalculator.la \
	../src/force/source/libforce.la \
	../src/function_parser/source/libfunction_parser.la \
	../src/geometry/source/libgeometry.la \
	../src/geometry/stl/source/libstl.la \
	../src/integrator/source/libintegrator.la \
	../src/meter/source/libmeter.la \
	../src/meter/grid_meter/source/libgrid_meter.la \
	../src/particle_creator/source/libparticle_creator.la \
	../src/postprocessor/source/libpostprocessor.la \
	../src/reflector/source/libreflector.la \
	../src/thermostat/source/libthermostat.la \
	../src/weighting_function/source/libweighting_function.la

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


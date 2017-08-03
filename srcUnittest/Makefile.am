bin_PROGRAMS = symplerTest

all:
	./symplerTest

symplerTest_SOURCES = mainTest.cpp

INCLUDES = \
	-I$(top_srcdir)/srcUnittest/*/head \
	-I$(top_srcdir)/srcUnittest/*/*/head \
	-I$(top_srcdir)/src/*/head \
	-I$(top_srcdir)/src/*/*/head

symplerTest_LDADD = \
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
	../src/*/source/lib*.la \
	../src/*/*/source/lib*.la \
	\
	../src/*/source/*.o \
	../src/*/*/source/*.o \
	\
	*/*/*/*.o \
	*/*/*.o

symplerTest_CXXFLAGS = \
	@GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

symplerTest_DEPENDENCIES = \
	mainTest.cpp \
	\
	*/source/lib*.la \
	*/*/source/lib*.la \
	\
	../src/*/source/lib*.la \
	../src/*/*/source/lib*.la

EXTRA_DIST = README

SUBDIRS = \
	basic \
	geometry \
	integrator \
	meter

# The following list as complete (2017-08-03). Add directories to upper SUBDIRS as soon as you have at least one compilable test cpp file
#SUBDIRS = \
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


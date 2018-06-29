
noinst_LTLIBRARIES = libintegratorTest.la

libintegratorTest_la_SOURCES = 	\
	integrator_positionTest.cpp \
	integrator_pos_vel_step2Test.cpp \
	integrator_velocity_verletTest.cpp


libintegratorTest_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

libintegratorTest_la_LDFLAGS = @LT_LDFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/srcUnittest/integrator/head -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head


noinst_LTLIBRARIES = libmeter.la

libmeter_la_SOURCES = \
	grid_averager.cpp \
	grid_averager_circular.cpp \
	grid_averager_structured.cpp \
	grid_meter.cpp \
	meter_autocorrelation_vector.cpp \
	meter_autocorrelation_vector_t.cpp \
	meter_average.cpp \
	meter_bonded_crosscorrelation_vector.cpp \
	meter_distribution.cpp \
	meter_grouped.cpp \
	meter_live.cpp \
	meter_live_coloured.cpp \
	meter_live_geomview.cpp \
	meter_pair_distribution.cpp \
	meter_pair_distribution_with_walls.cpp \
	meter_pos_vel.cpp \
	meter_relative_velocity.cpp 

libmeter_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/meter/grid_meter/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head


noinst_LTLIBRARIES = libgrid_meter.la

libgrid_meter_la_SOURCES = \
	gm_array_attr.cpp \
	gm_density.cpp \
	gm_ekin.cpp \
	gm_scalar.cpp \
	gm_scalar_fourier_cos.cpp \
	gm_scalar_fourier_sin.cpp \
	gm_scalar_func.cpp \
	gm_scalar_sum.cpp \
	gm_scalar_times_scalar_sum.cpp \
	gm_local_density.cpp \
	gm_pressure.cpp \
	gm_self_diffusion.cpp \
	gm_surface_tension.cpp \
	gm_temperature.cpp \
	gm_temperature_ie.cpp \
	gm_vector.cpp \
	gm_velocity.cpp \
	gm_velocity_moment.cpp \
	gm_volume_fraction.cpp 

libgrid_meter_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/meter/grid_meter/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head

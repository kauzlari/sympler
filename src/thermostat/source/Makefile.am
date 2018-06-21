
noinst_LTLIBRARIES = libthermostat.la

libthermostat_la_SOURCES = \
	apply_scalar_field.cpp \
	apply_vector_field.cpp \
	apply_vector_field_file.cpp \
	apply_vel_field.cpp \
	cbl_pair_particle_arbitrary.cpp \
	cbl_pair_particle_tensor.cpp \
	cbl_pair_particle_vector.cpp \
	shift_particle.cpp \
	thermostat.cpp \
	thermostat_energy_rescaling.cpp \
	thermostat_la.cpp \
	thermostat_la_ec.cpp \
	thermostat_la_with_weight_ec.cpp \
	thermostat_peters.cpp \
	thermostat_peters_ec.cpp \
	thermostat_peters_iso.cpp \
	thermostat_vels.cpp \
	thermostat_wall_vel.cpp \
	thermostat_wall_vel_ec.cpp \
	thermostat_wall_vel_iso.cpp \
	thermostat_with_rng.cpp

libthermostat_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/thermostat/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/force/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head


noinst_LTLIBRARIES = libforce.la

libforce_la_SOURCES = \
	connector_basic.cpp \
	f_angular.cpp \
	f_curvature.cpp \
	f_centrifugal.cpp \
	f_coriolis.cpp \
	f_dpd.cpp \
	f_dpde.cpp \
	f_generic.cpp \
	f_kinetic.cpp \
	f_mdpd.cpp \
	f_pair.cpp \
	f_pair_arbitrary.cpp \
	f_pair_arbitrary_BC.cpp \
	f_pair_arbitrary_wf.cpp \
	f_pair_scalar.cpp \
	f_pair_tensor.cpp \
	f_pair_vector.cpp \
	f_pair_vels.cpp \
	f_pair_vels_BC.cpp \
	f_pair_vels_wf.cpp \
	f_pair_wf.cpp \
	f_particle.cpp \
	f_particle_scalar.cpp \
	f_particle_vector.cpp \
	f_particle_vector_matrix.cpp \
	f_particle_vector_rand_matrix.cpp \
	f_particle_vels.cpp \
	f_particle_vels_matrix.cpp \
        f_particle_tensor.cpp \
	f_rand.cpp \
	f_rand_var_noise.cpp \
	f_specific.cpp\
	f_wall_dpd.cpp \
	f_wall_repulsion.cpp \
	f_with_rng.cpp \
	gen_connector.cpp \
	gen_f_solo.cpp \
	gen_triplet.cpp \
	gen_quintet.cpp \
	ie_heat_conduction.cpp \
	ie_wall_heat_conduction.cpp \
	lennard_jones.cpp \
	lennard_jones_vc.cpp

libforce_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/force/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head

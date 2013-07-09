
noinst_LTLIBRARIES = libparticle_creator.la

libparticle_creator_la_SOURCES = \
    pc_connector.cpp \
	pc_file.cpp \
	pc_free.cpp \
	pc_free_f.cpp \
	pc_free_pcalc.cpp \
	pc_inlet.cpp \
	pc_lattice.cpp \
	pc_lattice_bcc.cpp \
	pc_lattice_frozen.cpp \
	pc_lse.cpp \
	pc_random.cpp \
	pc_tube.cpp \
	pc_wall.cpp \
	pc_wall_textured.cpp \
	pc_static.cpp \
	pc_with_rng.cpp

libparticle_creator_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/boundary/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/particle_creator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/force/head

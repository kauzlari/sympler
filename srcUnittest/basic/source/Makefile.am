noinst_LTLIBRARIES = libbasicTest.la
libbasicTest_la_SOURCES = \
	pca_iapws-if97Test.cpp \
	pca_iapws-if97_1varTest.cpp \
	pca_iapws-if97_2varTest.cpp \
	pca_iapws-if97_cpTest.cpp \
	pca_iapws-if97_etaTest.cpp \
	pca_iapws-if97_kappaTest.cpp \
	pca_iapws-if97_pTest.cpp \
	pca_iapws-if97_rhoTest.cpp \
	fake_function_particle.cpp \
	particle_velsTest.cpp \
	smart_enumTest.cpp \
	symbol_f_particle_velsTest.cpp

libbasicTest_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

libbasicTest_la_LDFLAGS = @LT_LDFLAGS@

AM_CPPFLAGS = -I$(top_srcdir)/srcUnittest/basic/head -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/particle_creator/head -I$(top_srcdir)/src/reflector/head

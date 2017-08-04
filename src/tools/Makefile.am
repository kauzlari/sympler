bin_PROGRAMS = vtk2profile

vtk2profile_SOURCES = \
	vtk2profile.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/geometry/head

vtk2profile_LDADD = \
	-L../basic/source \
	-lbasic \
	-lstdc++

#vtk2profile_DEPENDENCIES = \
#	../basic/source/libbasic.la


dnl 

dnl AC_PREREQ(2.57)
AC_INIT(sympler, 0.1, david.kauzlaric@imtek.uni-freiburg.de)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([src/main.cpp])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_LANG(C++)


AC_DISABLE_SHARED
AM_DISABLE_SHARED

AM_PROG_LIBTOOL
dnl AC_PROG_RANLIB

XML_CFLAGS="DO_NOT_EXIST"

# Checks for libraries.
AM_PATH_GSL(1.3, :, AC_MSG_ERROR("-----------> GSL not found but required for sympler"))
AM_PATH_XML2(2.0.0, :, AC_MSG_ERROR("-----------> libxml2 not found but required for sympler"))

if test "${XML_CFLAGS}" = "DO_NOT_EXIST" ; then
   XML_CFLAGS="${XML_CPPFLAGS}"
   AC_SUBST(XML_CFLAGS)
fi
AC_CHECK_LIB(z, gzread)

AC_C_BIGENDIAN(ENDIAN="big", ENDIAN="little")

#AC_SUBST(AM_CPPFLAGS)

SDL_CFLAGS=""
AC_SUBST(SDL_CFLAGS)
SDL_LIBS=""
AC_SUBST(SDL_LIBS)

CHECK_SDL="yes"
HAVE_SDL="no"

# Static compile??
AC_ARG_ENABLE(all-static,
  [  --enable-all-static     build completely static (standalone) binaries], 
  [ 
    if test "$enableval" = "yes" ; then 
      LT_LDFLAGS="-all-static $LT_LDFLAGS" 
      CHECK_SDL="no"
    fi
  ]) 

USER_SDL="no"

# SDL-activation by user (Since I don't understand autotools, this is the only solution I am able to implement)
AC_ARG_ENABLE(meter-lives,
  [  --enable-meter-lives     enables 2 MeterLive modules and requires SDL-library >=1.2.0 and <2.0], 
  [ 
    if test "$enableval" = "yes" ; then  
      USER_SDL="yes"
    fi
  ]) 

AC_SUBST(LT_LDFLAGS)

if test "${CHECK_SDL}" = "no"; then
   HAVE_SDL="no"
else
   if test "${USER_SDL}" = "yes"; then
      HAVE_SDL="yes"
# Since I don't understand autotools, the next does not work if SDL not installed: therefore commented out
#      AM_PATH_SDL(1.2.0, HAVE_SDL="yes", HAVE_SDL="no")
      # Next is hopefully true on most systems. Don't know better since I don't understand...
      SDL_LIBS="-lSDL"
      AC_SUBST(SDL_LIBS)
   fi
fi

if test "${HAVE_SDL}" = "yes"; then
  AC_DEFINE(HAVE_SDL, 1, Is the SDL library present?)
  AC_MSG_WARN(Compiling with SDL library for modules such as MeterLive.)
  AC_MSG_WARN(Make sure you have libsdl1.2 or newer but NOT libsdl2.0 or newer!)
else
  AC_MSG_WARN(Compiling without MeterLive)
fi


# Use pthreads?
#AC_ARG_ENABLE(pthreads,
#  [  --enable-pthreads       build multithreaded code using pthreads (currently unstable)],
#  [
#    if test "$enableval" = "yes"; then
#      ENABLE_PTHREADS="yes"
#    fi
#  ])
#
#if test "${ENABLE_PTHREADS}" = "yes"; then
#  AC_DEFINE(ENABLE_PTHREADS, 1, Compile multithreaded code using pthreads)
##  CFLAGS="-pthread -D__USE_MALLOC $CFLAGS"
#fi

#AC_SUBST(CFLAGS)

# Use TNT library for array types?
AC_ARG_WITH(tnt,
	AS_HELP_STRING([--with-tnt],
	[use TNT library and enable matrix computations @<:@default=no@:>@] ),
	[],
	[ with_tnt=no ] )
AS_IF([test "x$with_tnt" != xno],
		[AC_CHECK_HEADERS([tnt/tnt.h jama/jama_lu.h], [],
	  	[ echo "TNT/JAMA libraries missing" ] )])

# Test for library path. Taken from dune-project.org..
AC_DEFUN([_slu_lib_path],
    [
        my_include_path=include/superlu
        my_lib_path=lib
        my_slu_found=yes
        if test ! -f "$1/$my_include_path/$2" ; then
            #Try to find headers under superlu
            my_include_path=include
            if test ! -f "$with_superlu/$my_include_path/$2" ; then
#                my_include_path=SRC
#                my_lib_path=""
#                if test ! -f "$with_superlu/$my_include_path/$2"; then
                    my_slu_found=no
#               fi
            fi
        fi
    ]
)

AC_DEFUN([_slu_search_versions],
    [
        my_slu_header=slu_ddefs.h
        _slu_lib_path($1, "$my_slu_header")
        if test "$my_slu_found" != "yes"; then 
            my_slu_header="dsp_defs.h"
            _slu_lib_path($1, "$my_slu_header")
        fi
    ]
)

AC_DEFUN([_slu_search_default],
    [
        with_superlu=/usr
        _slu_search_versions($with_superlu)
        
        if test "$my_slu_found" = "no"; then
            with_superlu=/usr/local
            _slu_search_versions($with_superlu)
        fi
    ]
)

AC_DEFUN([SYMPLER_PATH_SUPERLU],[
  #
  # User hints ...
  #
        my_lib_path=""
        my_include_path=""
        AC_ARG_WITH([superlu],
            [AC_HELP_STRING([--with-superlu],[user defined path to SuperLU library])],
            [dnl
                if  test "$withval" != no ; then
                    # get absolute path
                    with_superlu=`eval cd $withval 2>&1 && pwd`
                    if test "$withval" = yes; then
                        # Search in default locations
                        _slu_search_default
                    else
                        # Search for the headers in the specified location
                        with_superlu="$withval"
                        _slu_search_versions("$with_superlu")
                    fi
                fi
                ],
            [dnl
                # Search in default locations
                    _slu_search_default
                ])
  
        AC_ARG_WITH([super_lu_lib],
            [AC_HELP_STRING([--with-superlu-lib],[The name of the static SuperLU library to link to. By default the shared library with the name superlu-mpi is tried])],
            [
                if test "$withval" != no ; then
                    with_superlu_lib=$withval
                fi
            ]
        )
  # store old values
        ac_save_LDFLAGS="$LDFLAGS"
        ac_save_CPPFLAGS="$CPPFLAGS"
        ac_save_LIBS="$LIBS"
        
  ## do nothing if --without-superlu is used
        if test x"$with_superlu" != x"no" ; then
      # defaultpath
            SUPERLU_LIB_PATH="$with_superlu/$my_lib_path"
            SUPERLU_INCLUDE_PATH="$with_superlu/$my_include_path"
            
            SUPERLU_LDFLAGS="-L$SUPERLU_LIB_PATH"
            
      # set variables so that tests can use them
            CPPFLAGS="$CPPFLAGS -I$SUPERLU_INCLUDE_PATH"
            
      # check for central header
            AC_CHECK_HEADER([$my_slu_header],[
                    SUPERLU_CPPFLAGS="$CPPFLAGS"
                    HAVE_SUPERLU="1"],[
                    HAVE_SUPERLU="0"
                    AC_MSG_WARN([$my_slu_header not found in $SUPERLU_INCLUDE_PATH with $CPPFLAGS])]
            )
            
      # if header is found check for the libs
            
            if test x$HAVE_SUPERLU = x1 ; then

                OLDFLAGS="$LDFLAGS"
                LDFLAGS="$LDFLAGS -L$SUPERLU_LIB_PATH"
                LIBS="$LIBS $FLIBS"
                HAVE_SUPERLU=0

                if test x$with_superlu_lib = x; then
                    AC_CHECK_LIB(superlu, [dgssvx],[
                            SUPERLU_LIBS="-lsuperlu $LIBS"
                            SUPERLU_LDFLAGS="$LDFLAGS"
                            HAVE_SUPERLU="1"
                            ],[
                            HAVE_SUPERLU="0"
                            AC_MSG_WARN(libsuperlu not found)])
                fi
                if test "$HAVE_SUPERLU" = 0; then
                    if test x$with_superlu_lib = x ; then
                        with_superlu_lib=superlu.a
                    fi
                    AC_MSG_CHECKING([static superlu library "$with_superlu_lib" in "$SUPERLU_LIB_PATH"])
                    if test -f "$SUPERLU_LIB_PATH/$with_superlu_lib" ; then
                        LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $LIBS"
                        AC_CHECK_FUNC(dgssvx,
                            [
                                SUPERLU_LDFLAGS="$OLDFLAGS"
                                SUPERLU_LIBS="$LIBS"
                                HAVE_SUPERLU="1"
                                AC_MSG_RESULT(yes)
                            ],
                            [
                                HAVE_SUPERLU="0"
                                AC_MSG_RESULT(failed)
                            ]
                        )
                    else
                        HAVE_SUPERLU="0"
                        AC_MSG_RESULT(failed)
                    fi
                fi
            fi
      # pre-set variable for summary
      #with_superlu="no"
      
      # did it work?
            AC_MSG_CHECKING([SuperLU in $with_superlu])
            if test x$HAVE_SUPERLU = x1 ; then
                AC_SUBST(SUPERLU_LDFLAGS, $SUPERLU_LDFLAGS)
                AC_SUBST(SUPERLU_LIBS, $SUPERLU_LIBS)
                AC_SUBST(SUPERLU_CPPFLAGS, $SUPERLU_CPPFLAGS)
                AC_DEFINE(HAVE_SUPERLU, 1, [Define to 1 if SUPERLU is found])
                if test "$my_slu_header" = "slu_ddefs.h"; then
                    AC_DEFINE(SUPERLU_POST_2005_VERSION, 1, [define to 1 if there is  a header slu_ddefs.h in SuperLU])
                fi
                AC_MSG_RESULT(ok)
                
    # re-set variable correctly
                with_superlu="yes"
            else
                with_superlu="no"
                AC_MSG_RESULT(failed)
            fi 
            
  # end of "no --without-superlu"
        else
            with_superlu="no"
        fi
        
  # tell automake       
        AM_CONDITIONAL(SUPERLU, test x$HAVE_SUPERLU = x1)

  # restore variables
        LDFLAGS="$ac_save_LDFLAGS"
        CPPFLAGS="$ac_save_CPPFLAGS"
        LIBS="$ac_save_LIBS"

    ]
)

SYMPLER_PATH_SUPERLU
if test x$HAVE_SUPERLU = x1 ; then
 LDFLAGS="$SUPERLU_LDFLAGS"
 CPPFLAGS="$SUPERLU_CPPFLAGS"
 LIBS="$SUPERLU_LIBS"
fi
# end of check for superlu

if test "${ENDIAN}" = "big"; then
   AC_DEFINE(IS_BIG_ENDIAN, 1, Is the system big endian?)
fi

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([float.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
dnl AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt])

# SVN revision number from good old subversion times
# see http://www.freesoftwaremagazine.com/books/agaal/catalog_of_reusable_solutions
# SVNREV=`svnversion $srcdir | sed 's/:.*//'`
# which svnversion > /dev/null; \
# if test $? -ne 0  || test "x$SVNREV" = "xexported"
#   then SVNREV=`cat $srcdir/SVNREV`
#   else echo -n $SVNREV>$srcdir/SVNREV
# fi
# AC_SUBST(SVNREV)

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 src/basic/Makefile
		 src/basic/head/Makefile
		 src/basic/source/Makefile
		 src/boundary/Makefile
		 src/boundary/head/Makefile
		 src/boundary/source/Makefile
		 src/calculator/Makefile
		 src/calculator/head/Makefile
		 src/calculator/source/Makefile
		 src/force/Makefile
		 src/force/head/Makefile
		 src/force/source/Makefile
		 src/function_parser/Makefile
		 src/function_parser/head/Makefile
		 src/function_parser/source/Makefile
		 src/geometry/Makefile
		 src/geometry/head/Makefile
		 src/geometry/source/Makefile
		 src/geometry/stl/Makefile
		 src/geometry/stl/head/Makefile
		 src/geometry/stl/source/Makefile
		 src/integrator/Makefile
		 src/integrator/head/Makefile
		 src/integrator/source/Makefile
		 src/meter/Makefile
		 src/meter/head/Makefile
		 src/meter/source/Makefile
		 src/meter/grid_meter/Makefile
		 src/meter/grid_meter/head/Makefile
		 src/meter/grid_meter/source/Makefile
		 src/particle_creator/Makefile
		 src/particle_creator/head/Makefile
		 src/particle_creator/source/Makefile
		 src/postprocessor/Makefile
		 src/postprocessor/head/Makefile
		 src/postprocessor/source/Makefile
		 src/reflector/Makefile
		 src/reflector/head/Makefile
		 src/reflector/source/Makefile
		 src/thermostat/Makefile
		 src/thermostat/head/Makefile
		 src/thermostat/source/Makefile
		 src/weighting_function/Makefile
		 src/weighting_function/head/Makefile
		 src/weighting_function/source/Makefile
		 src/tools/Makefile

		 srcUnittest/Makefile
		 srcUnittest/basic/Makefile
		 srcUnittest/basic/head/Makefile
		 srcUnittest/basic/source/Makefile
		 srcUnittest/boundary/Makefile
		 srcUnittest/boundary/head/Makefile
		 srcUnittest/boundary/source/Makefile
		 srcUnittest/calculator/Makefile
		 srcUnittest/calculator/head/Makefile
		 srcUnittest/calculator/source/Makefile
		 srcUnittest/force/Makefile
		 srcUnittest/force/head/Makefile
		 srcUnittest/force/source/Makefile
		 srcUnittest/function_parser/Makefile
		 srcUnittest/function_parser/head/Makefile
		 srcUnittest/function_parser/source/Makefile
		 srcUnittest/geometry/Makefile
		 srcUnittest/geometry/head/Makefile
		 srcUnittest/geometry/source/Makefile
		 srcUnittest/geometry/stl/Makefile
		 srcUnittest/geometry/stl/head/Makefile
		 srcUnittest/geometry/stl/source/Makefile
		 srcUnittest/integrator/Makefile
		 srcUnittest/integrator/head/Makefile
		 srcUnittest/integrator/source/Makefile
		 srcUnittest/meter/Makefile
		 srcUnittest/meter/head/Makefile
		 srcUnittest/meter/source/Makefile
		 srcUnittest/meter/grid_meter/Makefile
		 srcUnittest/meter/grid_meter/head/Makefile
		 srcUnittest/meter/grid_meter/source/Makefile
		 srcUnittest/particle_creator/Makefile
		 srcUnittest/particle_creator/head/Makefile
		 srcUnittest/particle_creator/source/Makefile
		 srcUnittest/postprocessor/Makefile
		 srcUnittest/postprocessor/head/Makefile
		 srcUnittest/postprocessor/source/Makefile
		 srcUnittest/reflector/Makefile
		 srcUnittest/reflector/head/Makefile
		 srcUnittest/reflector/source/Makefile
		 srcUnittest/thermostat/Makefile
		 srcUnittest/thermostat/head/Makefile
		 srcUnittest/thermostat/source/Makefile
		 srcUnittest/weighting_function/Makefile
		 srcUnittest/weighting_function/head/Makefile
		 srcUnittest/weighting_function/source/Makefile
		 srcUnittest/tools/Makefile])
AC_OUTPUT

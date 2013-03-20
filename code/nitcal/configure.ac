dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(nitcal, 0.0.0, pbauman@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)]) # use silent rules where available - automake 1.11

LT_INIT


dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=nitcal

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])

dnl------------------------------
dnl Checks for compilers and MPI
dnl------------------------------

AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

dnl Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

dnl----------------
dnl Libmesh Config
dnl----------------
AX_PATH_LIBMESH_NEW(0.8.0, yes)

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--enable-libmesh-flags],[enable use of libmesh cxxflags for CXXFLAGS]),
               USE_LIBMESH_CXX_FLAGS=1
               AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled]),[])


if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   CXXFLAGS=`$LIBMESH_CONFIG --cxxflags --cppflags`
else
 dnl TODO: WE NEED TO DISTINGUISH BETWEEN COMPILERS IN A BETTER WAY
   CXXFLAGS="-O0 -g -Wall"
fi

dnl---------------------------------------------------------
dnl Add libMesh flags manually if it's not a libtool build
dnl---------------------------------------------------------
if test "x$LIBMESH_LIBTOOL" == "x1"; then
  LIBNITCAL_LIBS='$(top_builddir)/src/libnitcal.la'
else
  LIBNITCAL_LIBS='$(top_builddir)/src/libnitcal.la $(LIBMESH_LDFLAGS)'                                                                                                                                                                       
fi
AC_SUBST(LIBNITCAL_LIBS)

dnl----------------
dnl Check for boost
dnl----------------
BOOST_REQUIRE([1.39.0])


AX_PATH_GRINS(0.2.0,yes)

#----------------
# Check for QUESO
#----------------
AX_PATH_QUESO_NEW(0.45.0, no)


dnl---------------------------------------------------------
dnl Checks for  library prerequisites for other libraries...
dnl---------------------------------------------------------

dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE

dnl---------------------------------
dnl Query configuration environment
dnl Must be done *before* Doxygen to
dnl get some happy variables defined.
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(VAR, doxygen/nitcal.dox, docs/doxygen )


dnl-------------------
dnl Generate makefiles
dnl-------------------
AC_OUTPUT( Makefile \
           src/Makefile \
           test/Makefile \
           doxygen/Makefile \
           doxygen/txt_common/about.page \
           doxygen/nitcal.dox )

#           doxygen/Makefile \
#           doxygen/txt_common/about.page )


dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG

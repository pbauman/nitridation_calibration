# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2010-03-24
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo SVN revision number........... : $BUILD_VERSION
echo
echo Library Dependencies:
echo libMesh....................... : $LIBMESH_PREFIX
echo libMesh CXXFLAGS.............. : $LIBMESH_CXXFLAGS
echo libMesh INCLUDE............... : $LIBMESH_INCLUDE
#echo QUESO......................... : $QUESO_PREFIX
#echo Trilinos...................... : $TRILINOS_PREFIX
#echo HDF5.......................... : $HDF5_PREFIX
#echo GSL........................... : $GSL_PREFIX
#echo GLPK.......................... : $GLPK_PREFIX
echo Boost......................... : $BOOST_ROOT
echo Antioch....................... : $ANTIOCH_PREFIX

# masa optional check:
echo
echo Optional Features:
if test "x$HAVE_GSL" = "x1"; then
  echo '   'Link with GSL................. : $GSL_PREFIX
else
  echo '   'Link with GSL................. : no
fi
if test "x$HAVE_GLPK" = "x1"; then
  echo '   'Link with GLPK................ : $GLPK_PREFIX
else
  echo '   'Link with GLPK................ : no
fi
if test "$HAVE_QUESO" = "1"; then
  echo '   'Link with QUESO............... : $QUESO_PREFIX
else
  echo '   'Link with QUESO............... : no
fi
if test "x$HAVE_GRVY" = "x1"; then
  echo '   'Link with GRVY................ : $GRVY_PREFIX
else
  echo '   'Link with GRVY................ : no
fi
if test "x$USE_GRVY_TIMERS" = "x1"; then
  echo '   'Use GRVY timers............... : yes
else
  echo '   'Use GRVY timers............... : no
fi
if test "$xHAVE_MASA" = "x1"; then
  echo '   'Link with MASA................ : $MASA_PREFIX
else
  echo '   'Link with MASA................ : no
fi
if test "$HAVE_GCOV_TOOLS" = "0"; then
  echo '   'Enable gcov code coverage..... : no
else     
  echo '   'Enable gcov code coverage..... : yes
fi


echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])

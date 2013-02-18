BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = nitfp

lib_LTLIBRARIES = libnitcal.la

#--------------------------------------------------
# List of directories with source and include files
#--------------------------------------------------
boundary_conditions_srcdir = $(top_srcdir)/src/boundary_conditions/src
boundary_conditions_incdir = $(top_srcdir)/src/boundary_conditions/include

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
boundary_conditions_SRC = $(boundary_conditions_srcdir)/nitcal_bc_factory.C \
                          $(boundary_conditions_srcdir)/tube_twall.C

#----------------------------
#INCLUDES we want distributed
#----------------------------
boundary_conditions_INC = $(boundary_conditions_incdir)/nitcal_bc_factory.h \
                          $(boundary_conditions_incdir)/tube_twall.h

include_HEADERS = $(boundary_conditions_INC)

#------------------------------------
# Source for stand alone program 
#------------------------------------

nitfp_SOURCES = $(top_srcdir)/src/apps/nitfp.C

nitfp_LDADD = libnitcal.la

libnitcal_la_SOURCES = $(boundary_conditions_SRC)

libnitcal_la_LIBADD = $(GRINS_PREFIX)/lib/libgrins.la $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(boundary_conditions_incdir)

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRINS_CPPFLAGS)


if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(boundary_conditions_INC) $(boundary_conditions_SRC)
	echo 'updated source license headers' >$@

EXTRA_DIST = common

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@

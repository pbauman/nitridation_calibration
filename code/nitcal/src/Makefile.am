BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = nitfp

lib_LTLIBRARIES = libnitcal.la

#----------------------------------------
# Source files for libnitcal
#----------------------------------------

# boundary_conditions
libnitcal_la_SOURCES  = boundary_conditions/src/nitcal_bc_factory.C
libnitcal_la_SOURCES += boundary_conditions/src/tube_twall.C

# qoi
libnitcal_la_SOURCES += qoi/src/mass_loss.C

#----------------------------------------
# Header files for libnitcal
#----------------------------------------

# boundary_conditions
include_HEADERS  = boundary_conditions/include/nitcal_bc_factory.h
include_HEADERS += boundary_conditions/include/tube_twall.h

# qoi
include_HEADERS += qoi/include/mass_loss.h

libnitcal_la_LIBADD = $(GRINS_PREFIX)/lib/libgrins.la $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la

#------------------------------------
# Source for nitfp
#------------------------------------

nitfp_SOURCES = apps/nitfp.C

nitfp_LDADD = libnitcal.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS  = -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRINS_CPPFLAGS)


if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
        $(top_srcdir)/src/apps/*.C \
	$(top_srcdir)/src/boundary_conditions/include/*.h $(top_srcdir)/src/boundary_conditions/src/*.C \
	$(top_srcdir)/src/qoi/include/*.h $(top_srcdir)/src/qoi/src/*.C
	echo 'updated source license headers' >$@

EXTRA_DIST = common

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@

BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS  = nitfp
bin_PROGRAMS += likelihood_test
bin_PROGRAMS += nitsip
bin_PROGRAMS += reaction_energy_test

lib_LTLIBRARIES = libnitcal.la

#----------------------------------------
# Source files for libnitcal
#----------------------------------------

# boundary_conditions
libnitcal_la_SOURCES  = boundary_conditions/src/nitcal_bc_factory.C
libnitcal_la_SOURCES += boundary_conditions/src/tube_twall.C

#interface
libnitcal_la_SOURCES += interface/src/simulation_interface.C
libnitcal_la_SOURCES += interface/src/nitridation_simulation.C

# inverse problem
libnitcal_la_SOURCES += inverse_problem/src/sip_base.C
libnitcal_la_SOURCES += inverse_problem/src/gamma_cn_sip.C
libnitcal_la_SOURCES += inverse_problem/src/likelihood_base.C
libnitcal_la_SOURCES += inverse_problem/src/gamma_cn_likelihood.C

# qoi
libnitcal_la_SOURCES += qoi/src/qoi_factory.C
libnitcal_la_SOURCES += qoi/src/mass_loss.C

#----------------------------------------
# Header files for libnitcal
#----------------------------------------

# boundary_conditions
include_HEADERS  = boundary_conditions/include/nitcal_bc_factory.h
include_HEADERS += boundary_conditions/include/tube_twall.h

#interface
include_HEADERS += interface/include/simulation_interface.h
include_HEADERS += interface/include/nitridation_simulation.h

# inverse problem
include_HEADERS += inverse_problem/include/sip_base.h
include_HEADERS += inverse_problem/include/gamma_cn_sip.h
include_HEADERS += inverse_problem/include/likelihood_base.h
include_HEADERS += inverse_problem/include/gamma_cn_likelihood.h

# qoi
include_HEADERS += qoi/include/qoi_factory.h
include_HEADERS += qoi/include/mass_loss.h

libnitcal_la_LIBADD = $(GRINS_PREFIX)/lib/libgrins.la $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la

if QUESO_ENABLED
   libnitcal_la_LIBADD += $(QUESO_PREFIX)/lib/libqueso.la
endif

#------------------------------------
# Source for nitfp
#------------------------------------

nitfp_SOURCES = apps/nitfp.C
nitfp_LDADD = libnitcal.la

#------------------------------------
# Source for likelihood_test
#------------------------------------

likelihood_test_SOURCES = apps/likelihood_test.C
likelihood_test_LDADD = libnitcal.la

#------------------------------------
# Source for nitsip
#------------------------------------

nitsip_SOURCES = apps/nitsip.C
nitsip_LDADD = libnitcal.la

#------------------------------------
# Source for reaction_energy_test
#------------------------------------

reaction_energy_test_SOURCES = apps/reaction_energy_test.C
reaction_energy_test_LDADD = libnitcal.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS  = -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/interface/include
AM_CPPFLAGS += -I$(top_srcdir)/src/inverse_problem/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRINS_CPPFLAGS)
AM_CPPFLAGS += $(QUESO_CPPFLAGS)

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

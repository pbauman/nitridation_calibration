BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS  = nitfp
bin_PROGRAMS += likelihood_test
bin_PROGRAMS += nitsip
bin_PROGRAMS += likelihood_sampling

lib_LTLIBRARIES = libnitcal.la

#----------------------------------------
# Source files for libnitcal
#----------------------------------------

# boundary_conditions
libnitcal_la_SOURCES  = boundary_conditions/src/nitcal_bc_factory.C
libnitcal_la_SOURCES += boundary_conditions/src/tube_twall.C
libnitcal_la_SOURCES += boundary_conditions/src/inlet_profile.C

#interface
libnitcal_la_SOURCES += interface/src/simulation_interface.C
libnitcal_la_SOURCES += interface/src/nitridation_simulation.C

# inverse problem
libnitcal_la_SOURCES += inverse_problem/src/sip_base.C
libnitcal_la_SOURCES += inverse_problem/src/gamma_cn_sip.C
libnitcal_la_SOURCES += inverse_problem/src/queso_likelihood_interface.C
libnitcal_la_SOURCES += inverse_problem/src/likelihood_base.C
libnitcal_la_SOURCES += inverse_problem/src/constant_gamma_cn_likelihood.C
libnitcal_la_SOURCES += inverse_problem/src/mass_loss_helper.C
libnitcal_la_SOURCES += inverse_problem/src/average_N_helper.C

# qoi
libnitcal_la_SOURCES += qoi/src/qoi_factory.C
libnitcal_la_SOURCES += qoi/src/mass_loss.C
libnitcal_la_SOURCES += qoi/src/mass_loss_catalytic.C
libnitcal_la_SOURCES += qoi/src/average_N_mole_fraction.C

#----------------------------------------
# Header files for libnitcal
#----------------------------------------

# boundary_conditions
include_HEADERS  = boundary_conditions/include/nitcal_bc_factory.h
include_HEADERS += boundary_conditions/include/tube_twall.h
include_HEADERS += boundary_conditions/include/inlet_profile.h

#interface
include_HEADERS += interface/include/simulation_interface.h
include_HEADERS += interface/include/nitridation_simulation.h

# inverse problem
include_HEADERS += inverse_problem/include/sip_base.h
include_HEADERS += inverse_problem/include/gamma_cn_sip.h
include_HEADERS += inverse_problem/include/queso_likelihood_interface.h
include_HEADERS += inverse_problem/include/likelihood_base.h
include_HEADERS += inverse_problem/include/constant_gamma_cn_likelihood.h
include_HEADERS += inverse_problem/include/mass_loss_helper.h
include_HEADERS += inverse_problem/include/average_N_helper.h

# qoi
include_HEADERS += qoi/include/qoi_factory.h
include_HEADERS += qoi/include/mass_loss.h
include_HEADERS += qoi/include/mass_loss_catalytic.h
include_HEADERS += qoi/include/average_N_mole_fraction.h

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
# Source for likelihood_sampling
#------------------------------------

likelihood_sampling_SOURCES = apps/likelihood_sampling.C
likelihood_sampling_LDADD = libnitcal.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS  = -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/interface/include
AM_CPPFLAGS += -I$(top_srcdir)/src/inverse_problem/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(ANTIOCH_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRINS_CPPFLAGS)
AM_CPPFLAGS += $(QUESO_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CPPFLAGS)
AM_CPPFLAGS += $(GLPK_CPPFLAGS)

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

BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS  = nitfp
bin_PROGRAMS += nitsip
bin_PROGRAMS += nitmap
bin_PROGRAMS += nit_interp
bin_PROGRAMS += nit_data

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

# inverse_problem - SIP stuff
libnitcal_la_SOURCES += inverse_problem/src/queso_sip_interface.C
libnitcal_la_SOURCES += inverse_problem/src/likelihood_comm_handler.C
libnitcal_la_SOURCES += inverse_problem/src/model_interpolation_builder.C
libnitcal_la_SOURCES += inverse_problem/src/constant_gamma_n_constant_gamma_cn_model.C
libnitcal_la_SOURCES += inverse_problem/src/arrhenius_gamma_n_constant_gamma_cn_model.C
libnitcal_la_SOURCES += inverse_problem/src/arrhenius_gamma_n_arrhenius_gamma_cn_model.C
libnitcal_la_SOURCES += inverse_problem/src/full_model_evaluator.C
libnitcal_la_SOURCES += inverse_problem/src/full_model_composition.C
libnitcal_la_SOURCES += inverse_problem/src/model_builder.C
libnitcal_la_SOURCES += inverse_problem/src/queso_gaussian_likelihood_diagonal_covariance_interface.C
libnitcal_la_SOURCES += inverse_problem/src/full_model_likelihood.C
libnitcal_la_SOURCES += inverse_problem/src/surrogate_model_composition.C
libnitcal_la_SOURCES += inverse_problem/src/surrogate_model_likelihood.C

# qoi
libnitcal_la_SOURCES += qoi/src/qoi_factory.C
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

# inverse_problem - SIP stuff
include_HEADERS += inverse_problem/include/queso_sip_interface.h
include_HEADERS += inverse_problem/include/likelihood_comm_handler.h
include_HEADERS += inverse_problem/include/model_interpolation_builder.h
include_HEADERS += inverse_problem/include/model_interface_base.h
include_HEADERS += inverse_problem/include/constant_gamma_n_constant_gamma_cn_model.h
include_HEADERS += inverse_problem/include/arrhenius_gamma_n_constant_gamma_cn_model.h
include_HEADERS += inverse_problem/include/arrhenius_gamma_n_arrhenius_gamma_cn_model.h
include_HEADERS += inverse_problem/include/model_evaluator_base.h
include_HEADERS += inverse_problem/include/full_model_evaluator.h
include_HEADERS += inverse_problem/include/full_model_composition.h
include_HEADERS += inverse_problem/include/model_builder.h
include_HEADERS += inverse_problem/include/queso_gaussian_likelihood_diagonal_covariance_interface.h
include_HEADERS += inverse_problem/include/full_model_likelihood.h
include_HEADERS += inverse_problem/include/surrogate_model_composition.h
include_HEADERS += inverse_problem/include/surrogate_model_likelihood.h

# qoi
include_HEADERS += qoi/include/qoi_names.h
include_HEADERS += qoi/include/qoi_factory.h
include_HEADERS += qoi/include/mass_loss_catalytic.h
include_HEADERS += qoi/include/average_N_mole_fraction.h

libnitcal_la_LIBADD = $(GRINS_PREFIX)/lib/libgrins.la $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la

if HDF5_ENABLED
   libnitcal_la_LIBADD += $(HDF5_PREFIX)/lib/libhdf5.la
endif

if QUESO_ENABLED
   libnitcal_la_LIBADD += $(QUESO_PREFIX)/lib/libqueso.la
endif

if GRVY_ENABLED
   libnitcal_la_LIBADD += $(GRVY_PREFIX)/lib/libgrvy.la
endif

#------------------------------------
# Source for nitfp
#------------------------------------

nitfp_SOURCES = apps/nitfp.C
nitfp_LDADD = libnitcal.la

#------------------------------------
# Source for nitsip
#------------------------------------

nitsip_SOURCES = apps/nitsip.C
nitsip_LDADD = libnitcal.la

#------------------------------------
# Source for nitmap
#------------------------------------
nitmap_SOURCES = apps/nitmap.C
nitmap_LDADD = libnitcal.la

#------------------------------------
# Source for nitmap
#------------------------------------
nit_interp_SOURCES = apps/nit_interp.C
nit_interp_LDADD = libnitcal.la

#------------------------------------
# Source for nitmap
#------------------------------------
nit_data_SOURCES = apps/nit_data.C
nit_data_LDADD = libnitcal.la

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

if HDF5_ENABLED
  AM_CPPFLAGS += $(HDF5_CFLAGS)
endif

AM_CPPFLAGS += $(GRVY_CFLAGS)

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

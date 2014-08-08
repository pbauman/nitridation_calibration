BUILT_SOURCES   = .license.stamp

check_PROGRAMS   = tube_twall_test
check_PROGRAMS  += fp_regression
check_PROGRAMS  += likelihood_regression
check_PROGRAMS  += two_dataset_likelihood_regression

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

AM_LDFLAGS = $(LIBNITCAL_LIBS)

tube_twall_test_SOURCES = tube_twall_test.C
fp_regression_SOURCES = fp_regression.C
likelihood_regression_SOURCES = likelihood_regression.C
two_dataset_likelihood_regression_SOURCES = two_dataset_likelihood_regression.C

#Define tests to actually be run
TESTS  = tube_twall_test.sh
TESTS += K_2_63.sh
TESTS += O_2_28.sh
TESTS += K_2_63_constant_gamma_n_constant_gamma_cn_likelihood_reg.sh
TESTS += O_2_28_constant_gamma_n_constant_gamma_cn_likelihood_reg.sh
TESTS += K_2_63_O_2_28_constant_gamma_n_constant_gamma_cn_likelihood_reg.sh
TESTS += K_2_63_O_2_28_constant_gamma_n_constant_gamma_cn_likelihood_reg2.sh

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

BUILT_SOURCES   = .license.stamp

check_PROGRAMS  = cubic_phase_function_test \
                  phi_step_conductivity_test

AM_CPPFLAGS = -I$(top_srcdir)/src/solver/include \
              -I$(top_srcdir)/src/physics/include \
              -I$(top_srcdir)/src/boundary_conditions/include \
              -I$(top_srcdir)/src/mat_props/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_INCLUDE) \
               $(GRVY_CFLAGS) \
               $(GRINS_CPPFLAGS)

AM_LDFLAGS = $(top_builddir)/src/libvar.la

cubic_phase_function_test_SOURCES = cubic_phase_function_test.C
phi_step_conductivity_SOURCES = phi_step_conductivity_test.C

# List of source files for license stamping
testsrc = $(cubic_phase_function_test_SOURCES) \
          $(phi_step_conductivity_SOURCES)

#Define tests to actually be run
TESTS =  cubic_phase_function_test \
         phi_step_conductivity_test

shellfiles_src = phase_function_test.sh

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src)

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(testsrc)
	echo 'updated source license headers' >$@
BUILT_SOURCES   = .license.stamp

check_PROGRAMS  = tube_twall_test

AM_CPPFLAGS = -I$(top_srcdir)/src/boundary_conditions/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_INCLUDE) \
               $(GRVY_CFLAGS) \
               $(GRINS_CPPFLAGS)

AM_LDFLAGS = $(LIBNITCAL_LIBS)

tube_twall_test_SOURCES = tube_twall_test.C

# List of source files for license stamping
testsrc = $(tube_twall_test_SOURCES)

#Define tests to actually be run
TESTS = tube_twall_test

shellfiles_src =

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src)

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(testsrc)
	echo 'updated source license headers' >$@

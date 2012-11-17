AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = var.dox var.page txt_common fig_common

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`


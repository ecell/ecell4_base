
SUBDIRS = boost.python #src tests

EXTRA_DIST = src tests

dist-hook:
	rm -rf `find $(distdir)/src -name .svn`
	rm -rf `find $(distdir)/boost.python -name .svn`
	rm -rf `find $(distdir)/tests -name .svn`
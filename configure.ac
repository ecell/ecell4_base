AC_REVISION([$Id: configure.in 2563 2006-05-28 03:00:28Z shafi $])
dnl
dnl
AC_INIT
AC_CONFIG_SRCDIR([autogen.sh])
dnl
AM_PROG_LIBTOOL
dnl AC_CONFIG_AUX_DIR()
dnl
dnl 
dnl 
dnl
AC_CANONICAL_TARGET([])
AM_INIT_AUTOMAKE([ecell],${ECELL_VERSION_NUMBER})
AC_PROG_MAKE_SET
dnl
#PYTHON_INCLUDE_DIR=$(shell python -c "import sys; print sys.prefix + '/include/python' + sys.version[:3]")
AM_CHECK_PYTHON_HEADERS()
AM_PATH_PYTHON()
dnl
rm -f conftest
./libtool --config > conftest
. ./conftest
rm -f conftest
dnl
dnl compatibility for libtool 1.5.6
LTDL_SHLIB_EXT=""
if test -n "$shrext_cmds"; then
    LTDL_SHLIB_EXT=$shrext_cmds
    AC_SUBST(LTDL_SHLIB_EXT)
dnl compatibility for libtool 1.5.0
elif test -n "$shrext"; then
    LTDL_SHLIB_EXT=$shrext
    AC_SUBST(LTDL_SHLIB_EXT)
fi
dnl
AC_CHECK_LIB(gslcblas,cblas_dgemm,,[AC_MSG_ERROR([could not find libgslcblas.])])
AC_CHECK_LIB(gsl,gsl_block_alloc,,[AC_MSG_ERROR([could not find libgsl.])])
dnl
AM_CHECK_PYMOD(numpy,,,[AC_MSG_ERROR([could not find Python numpy module.])])
dnl
dnl AC_LIBLTDL_INSTALLABLE
AM_CONFIG_HEADER([config.h])
dnl
dnl
dnl Create the top-level Makefile.
dnl
AC_CONFIG_FILES([Makefile])
dnl
AC_OUTPUT([])


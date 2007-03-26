AC_REVISION([$Id$])
dnl
dnl
AC_INIT
AC_CONFIG_SRCDIR([autogen.sh])
dnl
AM_PROG_LIBTOOL
dnl AC_CONFIG_AUX_DIR()
dnl
dnl 
AH_TEMPLATE(HAVE_SINCOS)
AH_TEMPLATE(HAVE_INLINE)
dnl 
dnl
AC_CANONICAL_TARGET([])
AM_INIT_AUTOMAKE([ecell],${ECELL_VERSION_NUMBER})
AC_PROG_MAKE_SET
AC_LANG(C++)
dnl
dnl checks for programs
dnl
AM_PATH_PYTHON()
dnl
dnl
dnl checks for libraries
dnl
dnl
AM_PATH_GSL(1.0,,AC_MSG_ERROR([could not find required version of GSL.]))
dnl
dnl
AM_CHECK_PYMOD(numpy,,,[AC_MSG_ERROR([could not find Python numpy module.])])
AM_CHECK_PYMOD(scipy,,,[AC_MSG_ERROR([could not find Python scipy module.])])
dnl
dnl
dnl
dnl checks for header files
dnl
dnl
dnl
#PYTHON_INCLUDE_DIR=$(shell python -c "import sys; print sys.prefix + '/include/python' + sys.version[:3]")
AM_CHECK_PYTHON_HEADERS()
dnl
dnl checks for types
dnl
dnl
dnl checks for structures
dnl
dnl
dnl checks for compiler characteristics
dnl
dnl
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
   AC_DEFINE(HAVE_INLINE,1)
   AC_SUBST(HAVE_INLINE)
fi
dnl
dnl extract LTDL_SHLIB_EXT
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
dnl checks for library functions
dnl
AC_CHECK_FUNC(sincos,AC_DEFINE(HAVE_SINCOS))
dnl
dnl
dnl no checks for system services
dnl
AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile test/Makefile])
AC_OUTPUT([])

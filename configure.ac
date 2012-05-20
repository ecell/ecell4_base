AC_REVISION([$Id$])
dnl
dnl
AC_INIT
AC_CONFIG_SRCDIR([autogen.sh])
dnl
AC_DISABLE_STATIC
AM_PROG_LIBTOOL
dnl AC_CONFIG_AUX_DIR()
dnl
dnl 
AH_TEMPLATE(HAVE_SINCOS)
AH_TEMPLATE(HAVE_INLINE)
dnl 
dnl
AC_CANONICAL_TARGET([])
AM_INIT_AUTOMAKE([epdp],[0.3b])

dnl Check pdflatex
AC_PROG_PDFLATEX()

dnl
AC_PROG_MAKE_SET
dnl
dnl checks for programs
dnl
AM_PATH_PYTHON(2.4)

dnl
dnl
dnl checks for libraries
dnl
dnl
AX_PATH_GSL([1.11],,AC_MSG_ERROR([could not find required version of GSL.]))
dnl
AC_CHECK_LIB(m,exp,,AC_MSG_ERROR([could not find libm.]))
dnl
AC_CHECK_LIB(python${PYTHON_VERSION},PySys_GetObject,,AC_MSG_ERROR([could not find libpython.]))
dnl
AC_PROG_CXX
AC_LANG_CPLUSPLUS
dnl
dnl

DEBUG=
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug],
                 [inhibit any optimization options]),
  [ DEBUG=1 ]
)

if test -n "$DEBUG"; then
  CFLAGS=""
  CXXFLAGS=""
  CPPFLAGS="-DDEBUG=1"
fi

if test "$GXX" = "yes"; then
  CXXFLAGS="$CXXFLAGS -Wall -Wstrict-aliasing=0 -Wno-invalid-offsetof"
  CFLAGS="$CFLAGS -Wall"
fi

AC_SUBST(DEBUG)

CXXFLAGS="$CXXFLAGS -g"
CFLAGS="$CFLAGS -g"

AX_BOOST_BASE([1.37],,AC_MSG_ERROR([could not find required version of BOOST.]))

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

BOOST_PYTHON_LIBNAME=boost_python

AC_ARG_WITH(boost-python-libname,
  AC_HELP_STRING([--with-boost-python-libname=LIBNAME],
                 [specify the library name of Boost.Python [[boost_python]]]
                ),
  [BOOST_PYTHON_LIBNAME=$withval]
)
AC_SUBST(BOOST_PYTHON_LIBNAME)

BOOST_REGEX_LIBNAME=boost_regex

AC_ARG_WITH(boost-regex-libname,
  AC_HELP_STRING([--with-boost-date-time-libname=LIBNAME],
                 [specify the library name of Boost.Python [[boost_regex]]]
                ),
  [BOOST_REGEX_LIBNAME=$withval]
)
AC_SUBST(BOOST_REGEX_LIBNAME)

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
ECELL_CHECK_NUMPY
ECELL_CHECK_NUMPY_ARRAY_DESCR
ECELL_CHECK_LOGGING_MODULE
dnl
AM_CHECK_PYTHON_HEADERS()
ac_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_CHECK_TYPES([PyBaseExceptionObject],[],[],[
#include <Python.h>
#include <pyerrors.h>
])
CPPFLAGS="$ac_save_CPPFLAGS"
AC_CHECK_FUNCS([PyInt_FromSize_t],[],[],[
#include <Python.h>
])
dnl
dnl
AC_CHECK_HEADERS([unordered_map boost/unordered_map.hpp boost/functional/hash.hpp])
AC_CHECK_HEADER([tr1/unordered_map], [
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_MSG_CHECKING([for a stupid bug in tr1/unordered_map])
  AC_TRY_COMPILE([
#include <tr1/unordered_map>
  ], [
std::tr1::unordered_map<int, int> a, b(a);
  ], [
    AC_MSG_RESULT([ok])
    AC_DEFINE([HAVE_TR1_UNORDERED_MAP], 1, [Define to 1 if you have usable <tr1/unordered_map> header file])
  ], [
    AC_MSG_RESULT([buggy])
  ])
  AC_LANG_RESTORE
], [
  AC_MSG_RESULT([buggy])
])

AC_CHECK_TYPES([std::hash<int>],[],[],[
#include <functional>
])
AC_CHECK_HEADERS([tr1/functional], [
  AC_CHECK_TYPES([std::tr1::hash<int>],[],[],[
#include <tr1/functional>
  ])
])
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
AC_CHECK_FUNCS([sincos isfinite])
AC_CHECK_DECLS([INFINITY, HUGE_VAL],,,[
#include <math.h>
])
dnl
dnl
dnl no checks for system services
dnl
AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile binding/Makefile test/Makefile doc/Makefile samples/benchmark/Makefile])
AC_OUTPUT([])

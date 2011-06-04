
## this one is commonly used with AM_PATH_PYTHON ...
dnl AM_CHECK_PYMOD(MODNAME [,SYMBOL [,ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]]])
dnl Check if a module containing a given symbol is visible to python.
AC_DEFUN([AM_CHECK_PYMOD],
[AC_REQUIRE([AM_PATH_PYTHON])
py_mod_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_MSG_CHECKING(for ifelse([$2],[],,[$2 in ])python module $1)
AC_CACHE_VAL(py_cv_mod_$py_mod_var, [
ifelse([$2],[], [prog="
import sys
try:
        import $1
except ImportError:
        sys.exit(1)
except:
        sys.exit(0)
sys.exit(0)"], [prog="
import $1
$1.$2"])
if $PYTHON -c "$prog" 1>&AC_FD_CC 2>&AC_FD_CC
  then
    eval "py_cv_mod_$py_mod_var=yes"
  else
    eval "py_cv_mod_$py_mod_var=no"
  fi
])
py_val=`eval "echo \`echo '$py_cv_mod_'$py_mod_var\`"`
if test "x$py_val" != xno; then
  AC_MSG_RESULT(yes)
  ifelse([$3], [],, [$3
])dnl
else
  AC_MSG_RESULT(no)
  ifelse([$4], [],, [$4
])dnl
fi
])

dnl a macro to check for ability to create python extensions
dnl  AM_CHECK_PYTHON_HEADERS([ACTION-IF-POSSIBLE], [ACTION-IF-NOT-POSSIBLE])
dnl function also defines PYTHON_INCLUDES
AC_DEFUN([AM_CHECK_PYTHON_HEADERS],
[AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING(for headers required to compile python extensions)
dnl deduce PYTHON_INCLUDES
py_prefix=`$PYTHON -c "import sys; print sys.prefix"`
py_exec_prefix=`$PYTHON -c "import sys; print sys.exec_prefix"`
PYTHON_INCLUDES="-I${py_prefix}/include/python${PYTHON_VERSION}"
if test "$py_prefix" != "$py_exec_prefix"; then
  PYTHON_INCLUDES="$PYTHON_INCLUDES -I${py_exec_prefix}/include/python${PYTHON_VERSION}"
fi
AC_SUBST(PYTHON_INCLUDES)
dnl check if the headers exist:
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_TRY_CPP([#include <Python.h>],dnl
[AC_MSG_RESULT(found)
$1],dnl
[AC_MSG_RESULT(not found)
$2])
CPPFLAGS="$save_CPPFLAGS"
])

dnl numpy package.
dnl find arrayobject.h.
dnl
AC_DEFUN([ECELL_CHECK_NUMPY], [
  AC_REQUIRE([AM_CHECK_PYTHON_HEADERS])

  AC_ARG_WITH(numpy-includes,
    AC_HELP_STRING([--with-numpy-includes=DIR],
                   [specify the numpy header location]),
    [NUMPY_INCLUDE_DIR=$withval],
    [NUMPY_INCLUDE_DIR=]
  )

  AC_MSG_CHECKING([for numpy include directory])
  if test -z "$NUMPY_INCLUDE_DIR"; then
    if ! NUMPY_INCLUDE_DIR=`$PYTHON -c "import numpy; print numpy.get_include();"`; then
      py_prefix=`$PYTHON -c "import sys; print sys.prefix"`
      pydir=python${PYTHON_VERSION}
      numpy_include="site-packages/numpy/core/include"
      EXT_GUESS= \
        "${py_prefix}/Lib/${numpy_include}" \
        "${py_prefix}/lib/${pydir}/${numpy_include}" \
        "${py_prefix}/lib64/${pydir}/${numpy_include}" \
        "/usr/lib/${pydir}/${numpy_include}" \
        "/usr/lib64/${pydir}/${numpy_include}" \
        "/usr/local/lib/${pydir}/${numpy_include}" \
        "/usr/local/lib64/${pydir}/${numpy_include}" \
        "${prefix}/include" \
        "/usr/include/${pydir}" \
        "/usr/local/include" \
        "/opt/numpy/include"
      NUMPY_INCLUDE_DIR=""
      for ac_dir in $EXT_GUESS ; do
        if test -f ${ac_dir}/numpy/arrayobject.h ; then
           NUMPY_INCLUDE_DIR=`(cd $ac_dir ; pwd)`
        fi
      done
    fi
  fi
  if test -z "${NUMPY_INCLUDE_DIR}"; then        
    AC_MSG_RESULT([not found in ${EXT_GUESS}.])
  else
    AC_MSG_RESULT(${NUMPY_INCLUDE_DIR})
  fi
  ac_save_CPPFLAGS="${CPPFLAGS}"
  CPPFLAGS="-I${NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDES}"
  AC_CHECK_HEADERS([numpy/arrayobject.h], [], [
    AC_MSG_ERROR([no usable NumPy headers were found. please check the installation of NumPy package.])
  ], [
#include <Python.h>
  ])
  CPPFLAGS="${ac_save_CPPFLAGS}"
  AC_SUBST(NUMPY_INCLUDE_DIR)
])

AC_DEFUN([ECELL_CHECK_NUMPY_ARRAY_DESCR], [
  AC_MSG_CHECKING([PyArray_Descr has hasobject])
  ac_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="-I${NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDES} $CPPFLAGS"
  AC_CHECK_MEMBER([PyArray_Descr.hasobject], [
    AC_DEFINE([PYARRAY_DESCR_HAS_HASOBJECT], [1], [Define to 1 if PyArray_Descr has hasobject field])
  ], [], [
#include <Python.h>
#include <numpy/arrayobject.h>
])
  CPPFLAGS="$ac_save_CPPFLAGS"                 
])

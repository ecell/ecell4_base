dnl @synopsis AC_PROG_PDFLATEX
dnl
dnl This macro test if pdflatex is installed. If pdflatex is installed,
dnl it set $pdflatex to the right value
dnl
dnl @category LaTeX
dnl @author Mathieu Boretti <boretti@bss-network.com>
dnl @version 2005-01-21
dnl @license GPLWithACException
dnl
dnl http://ac-archive.sourceforge.net/latex/ac_prog_pdflatex.html

AC_DEFUN([AC_PROG_PDFLATEX],[
AC_CHECK_PROGS(pdflatex,[pdflatex],no)
export pdflatex;
dnl if test $pdflatex = "no" ;
dnl then
dnl AC_MSG_ERROR([Unable to find a PDFLaTeX application]);
dnl fi
AM_CONDITIONAL([HAVE_PDFLATEX], [test "$pdflatex" != "no"])
AC_SUBST(pdflatex)
])

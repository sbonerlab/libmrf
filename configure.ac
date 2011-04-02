#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([libmrf], [1.0.0], [david.chen@yale.edu])
AC_CONFIG_SRCDIR([lib/mrf/mrf.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror])
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot find standard math library])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Cannot find cblas library])])
AC_CHECK_LIB([gsl], [gsl_ran_hypergeometric_pdf], [], [AC_MSG_ERROR([Cannot find gsl library])])
AC_CHECK_LIB([bios], [needMem], [], [AC_MSG_ERROR([Cannot find bios library])])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([strchr])

AC_CONFIG_FILES([Makefile lib/Makefile])
AC_OUTPUT
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([libmrf], [0.7.0], [david.z.chen@outlook.com])
AC_CONFIG_SRCDIR([mrf])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])

#------------------------------------------------------------------------------
# Initialize automake
#------------------------------------------------------------------------------
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#------------------------------------------------------------------------------
# Initialize libtool
#------------------------------------------------------------------------------
LT_INIT
AC_SUBST([MRF_SO_VESION], [0:7:0])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------
AC_PROG_SED
AC_PROG_MKDIR_P
AC_PROG_CC
AC_PROG_LIBTOOL

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------
AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot find standard math library])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Cannot find cblas library])])
AC_CHECK_LIB([gsl], [gsl_ran_hypergeometric_pdf], [], [AC_MSG_ERROR([Cannot find gsl library])])
AC_CHECK_LIB([bios], [needMem], [], [AC_MSG_ERROR([Cannot find bios library])])

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------
AC_CHECK_FUNCS([strchr])

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

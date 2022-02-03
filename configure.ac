AC_INIT([crackMesher], [0.1], [brian.phung@utah.edu])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CXX
AC_PROG_CC
AC_OPENMP

AX_CHECK_COMPILE_FLAG([-std=c++0x], [
                        CXXFLAGS="$CXXFLAGS -std=c++0x"])

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADER_STDBOOL
AC_CHECK_TYPES([ptrdiff_t])
AC_C_INLINE
AC_FUNC_MALLOC
AC_TYPE_SIZE_T

AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

# Checks for libraries.
AX_BOOST_BASE([1.48],, [AC_MSG_ERROR([This program needs Boost, but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS

# Python
AM_PATH_PYTHON([3.5], ,[AC_MSG_ERROR([Need at least python 3.5.])])
AC_ARG_VAR([PYTHON_INCLUDE], [Include flags for python, bypassing python-config])
AC_ARG_VAR([PYTHON_CONFIG], [Path to python-config])
AS_IF([test -z "$PYTHON_INCLUDE"], [
  AS_IF([test -z "$PYTHON_CONFIG"], [
    AC_PATH_PROGS([PYTHON_CONFIG],
                  [python$PYTHON_VERSION-config python-config],
                  [no],
                  [`dirname $PYTHON`])
    AS_IF([test "$PYTHON_CONFIG" = no], [AC_MSG_ERROR([cannot find python-config for $PYTHON.])])
  ])
  AC_MSG_CHECKING([python include flags])
  PYTHON_INCLUDE=`$PYTHON_CONFIG --includes`
  AC_MSG_RESULT([$PYTHON_INCLUDE])
])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

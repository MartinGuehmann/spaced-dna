AC_PREREQ([2.69])
AC_INIT([spaced], [1.0.2+master])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_CONFIG_SRCDIR([extkey.cxx])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_LANG(C++)
AC_OPENMP

# Check available standards.
AC_CONFIG_MACRO_DIR([m4])
AX_CXX_COMPILE_STDCXX_11([],[mandatory])


# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_INIT([main], [0.1], [chenyu@ices.utexas.edu])

# Checks for programs
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL

# Checks for librariesl
AC_CONFIG_MACRO_DIR([m4])

AX_PATH_PETSC([no])

AM_INIT_AUTOMAKE([gnu subdir-objects])

ACLOCAL_AMFLAGS = -I m4

# Checks for library functions
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
                 test/Makefile])

AC_OUTPUT

AC_INIT([phshape], [1.0], [Kenneth.Jansen@colorado.edu], 'phshape')
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([config/config.h])
AM_INIT_AUTOMAKE( [subdir-objects -Wall -Werror foreign ] )

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP
AM_PROG_CC_C_O

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

AX_DEBUG
AX_OPTIMIZATION
AX_PROFILING

#Simmetrix 
AX_SIMMODS
AS_IF([test "x$SIMMODS_MISSING" = xtrue], [echo "Simmetrix is required, specify --with-simmods=<dir>"; exit 1])
CXXFLAGS+="-DSIM "

# Checks for inline support
AC_C_INLINE

# Checks for programs.
AC_PROG_LIBTOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

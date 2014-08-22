AC_INIT([phParAdapt], [1.0], [Kenneth.Jansen@colorado.edu], 'phParAdapt')
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([config/config.h])
AM_INIT_AUTOMAKE( [subdir-objects -Wall -Werror foreign ] )

AX_MPI

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
OLD_FFLAGS=$FFLAGS
OLD_FCFLAGS=$FCFLAGS

AC_PROG_CC( [$CC_LIST] )
AC_PROG_CPP
AC_PROG_CXX( [$CXX_LIST] )
AC_PROG_CXXCPP
AC_PROG_FC( [$FC_LIST] )
AC_PROG_F77( [$F77_LIST] )
AM_PROG_CC_C_O

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS
FFLAGS=$OLD_FFLAGS
FCFLAGS=$OLD_FCFLAGS

# Checks for Fortran compiler
if test "x" != "x$FC"; then
  AC_FC_WRAPPERS
fi

AM_CONDITIONAL([PHPARA_PARALLEL], [test "x$MPI_MISSING" = "xfalse"])

AC_SUBST(MA_CFLAGS)
AC_SUBST(MA_LIBS)
if `pkg-config --atleast-version=1.0 libmeshadapt`; then
    MA_CFLAGS=`pkg-config --cflags libmeshadapt`
    MA_LIBS=`pkg-config --libs libmeshadapt`
    CXXFLAGS+=' -DFMDB '
else
    AC_MSG_ERROR([MeshAdapt package not found])
fi

AC_SUBST(PHSHAPE_CFLAGS)
AC_SUBST(PHSHAPE_LIBS)
if `pkg-config --atleast-version=1.0 libPhShape`; then
    PHSHAPE_CFLAGS=`pkg-config --cflags libPhShape`
    PHSHAPE_LIBS=`pkg-config --libs libPhShape`
else
    AC_MSG_ERROR([PhShape package not found])
fi

AC_SUBST(PHASTAIO_CFLAGS)
AC_SUBST(PHASTAIO_LIBS)
if `pkg-config --atleast-version=1.0 libPhastaIO`; then
    PHASTAIO_CFLAGS=`pkg-config --cflags libPhastaIO`
    PHASTAIO_LIBS=`pkg-config --libs libPhastaIO`
else
    AC_MSG_ERROR([PhastaIO package not found])
fi

AC_SUBST(LU_CFLAGS)
AC_SUBST(LU_LIBS)
if `pkg-config --atleast-version=1.0 libPhUtilLU`; then
    LU_CFLAGS=`pkg-config --cflags libPhUtilLU`
    LU_LIBS=`pkg-config --libs libPhUtilLU`
else
    AC_MSG_ERROR([PhUtilLU package not found])
fi

AC_SUBST(HPS_CFLAGS)
AC_SUBST(HPS_LIBS)
if `pkg-config --atleast-version=0.1.2 libHPS`; then
    HPS_CFLAGS=`pkg-config --cflags libHPS`
    HPS_LIBS=`pkg-config --libs libHPS`
else 
    AC_MSG_ERROR([HPS package not found])
fi

AC_SUBST(DIFFMC_CFLAGS)
AC_SUBST(DIFFMC_LIBS)
if `pkg-config --atleast-version=0.1.275 libDiffMC`; then
    DIFFMC_CFLAGS=`pkg-config --cflags libDiffMC`
    DIFFMC_LIBS=`pkg-config --libs libDiffMC`
else 
    AC_MSG_ERROR([ParMA DiffMC package not found])
fi

# Checks for inline support
AC_C_INLINE

# Checks for programs.
AC_PROG_RANLIB

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

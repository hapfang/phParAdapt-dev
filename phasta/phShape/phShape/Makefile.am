#bring in the macros defined in the m4 directory.. these include AX_GMI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libphshape.la
libphshape_la_SOURCES = \
  src/blend.cc \
  src/blendDrv2.cc \
  src/blendDrv.cc \
  src/Bn.c \
  src/BnDrv2.c \
  src/BnDrv.c \
  src/cequal.cc \
  src/En.c \
  src/EnDrv2.c \
  src/EnDrv.c \
  src/Fn.c \
  src/FnDrv2.c \
  src/FnDrv.c \
  src/localIndex.cc \
  src/modeShape.cc \
  src/modeShapeDrv.cc \
  src/newshape.cc \
  src/parDrv.cc \
  src/setQueryP.cc \
  src/shapeFunction.cc \
  src/shpDrv2.cc 

libphshape_la_include_HEADERS = \
  src/shapeFuncInternals.h \
  src/shapeFunction.h

libphshape_la_includedir = $(includedir)

AM_CXXFLAGS =  \
  $(DEBUG_FLAGS) \
  $(OPTIMIZATION_FLAGS) \
  $(PROFILING_FLAGS) \
  $(SIMMODS_INCLUDES) \
  -I$(top_srcdir) 

AM_CFLAGS = $(AM_CXXFLAGS)

SUBSYSNAME = phasta
MODULENAME = phShape
NVERS = 1.0

dirs := src

ifndef ARCHOS
  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifeq ($(FMDB),1)
  ifndef FMDBVERS
    FMDBVERS := FMDB/trunk
  endif
else
  SIM := 1
  ifndef MESHSIM
#    MESHSIM := /users/zhoum/meshSim/7.0-100403
    MESHSIM := /usr/local/meshSim/7.2-120125
  endif
endif

ifeq ($(SIM),1)
  DEFS := $(DEFS) -DSIM
  ifeq ($(ARCHOS),x86_cygwin)
    INCLUDES := $(INCLUDES) -I"$(MESHSIM)\include"
  else
    INCLUDES := $(INCLUDES) -I$(MESHSIM)/include
  endif
else
  DEFS := $(DEFS) -DFMDB
  DEPS := $(DEPS) $(FMDBVERS)/FMDB $(FMDBVERS)/SCORECUtil $(FMDBVERS)/SCORECModel
endif

include $(DEVROOT)/Util/buildUtil/make.common

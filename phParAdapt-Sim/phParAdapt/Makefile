SUBSYSNAME := phasta
MODULENAME := phParAdapt
ISEXECUTABLE := 1
OUTNAME := phParAdapt

ifeq ($(FMDB),1)
  DEFS := $(DEFS) -DFMDB
  ifndef FMDBVERS
    FMDBVERS := FMDB
  endif
  ifndef MESHADAPTVERS
    MESHADAPTVERS := meshAdapt
  endif
  ifndef MODELER
    MODELER := meshmodel
  endif
  ifeq ($(MODELER),parasolid)
    ifndef MESHSIM
      MESHSIM := /net/common/meshSim/latest
    endif
    ifndef PARASOLID
      PARASOLID := /usr/local/parasolid/latest
    endif
    MCTKVERS := mctk
  endif
  ifeq ($(MODELER),acis)
    ifndef ACIS
      ACIS := /usr/local/acis/R13
    endif
    MCTKVERS := mctk
  endif
else
  SIM := 1
  DEFS := $(DEFS) -DSIM
  ifndef MODELER
    MODELER := Discrete
  endif
  ifndef MESHSIM
	  MESHSIM := /usr/local/meshSim/7.2-120125
#    MESHSIM := /usr/local/meshSim/7.1-100901
  endif
endif

OUTNAME:=$(OUTNAME)-$(MODELER)

ifeq ($(HOSTNAME),master.medusa.scorec.rpi.edu)
  DEVROOT := /import/users/$(USER)/develop
endif

ifndef ARCHOS
  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifeq ($(ARCHOS),ibm-bg)
  SIMARCH := ppc64_rhel6_xlc
  DEFS := $(DEFS) -DPARALLEL -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
  DEFS := $(DEFS) -Dibm -DB64
  parallel := 1
endif

ifeq ($(ARCHOS),x86_64_linux-pgi)
  DEFS := $(DEFS) -DB64 -DPARALLEL -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
  LDFLAGS := $(LDFLAGS) -pgf90libs
  parallel := 1
endif

ifeq ($(ARCHOS),x86_64_linux-icc)
   SIMARCH := x64_rhel7_gcc48
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
   DEFS := $(DEFS) -DB64
endif

ifeq ($(ARCHOS),x86_64_linux)
   SIMARCH := x64_rhel7_gcc48
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
   DEFS := $(DEFS) -DB64
endif

ifeq ($(ARCHOS),x86_64_linux-gcc)
   SIMARCH := x64_rhel7_gcc48
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
   DEFS := $(DEFS) -DB64
endif
dirs := .

ifeq ($(ARCHOS),sun4_5)
  SIMARCH := sun4_5_sc5
  ETCFLAGS := $(ETCFLAGS) -lnsl -lm -lsocket -lsunperf -lF77 -lsunmath -lthread -lfsu
endif

ifeq ($(ARCHOS),sun4_5-64)
  SIMARCH := sun4_5_sc5-64
  ETCFLAGS := $(ETCFLAGS) -lnsl -lm -lsocket -lsunperf -lF77 -lsunmath -lthread -lfsu
endif

ifeq ($(ARCHOS),x86_linux-icc)
   SIMARCH := x86_rhel3_gcc32
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
endif

ifeq ($(ARCHOS),x86_linux)
   SIMARCH := x86_rhel3_gcc32
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
endif

ifeq ($(HOSTNAME),master.medusa.scorec.rpi.edu)
    SIMARCH := x86_rhel3_gcc32
    MESHSIM := /import/users/SCOREC/public/meshSim/latest/
    ETCFLAGS := $(ETCFLAGS)  -lnsl -lpthread
endif

ifeq ($(SIM),1)
 ifeq ($(PARALLEL),mpich1)
    SIMPAR := -mpich
 endif
 ifeq ($(PARALLEL),mvapich-1.2-gcc41)
    SIMPAR := -mpich
 endif
 ifeq ($(PARALLEL),openmpi)
   SIMPAR := -openmpi
 endif
 ifeq ($(PARALLEL),sgimpi)
  SIMPAR := -mpisgi
 endif
endif

ifeq ($(FMDB),1)
  IPCOMMAN := $(DEVROOT)/ipcomman
  INCLUDES := $(INCLUDES) -I$(IPCOMMAN)/include
  DEPS := $(DEPS) ipcomman
  LIBS := $(LIBS) ipcomman

 ifneq ($(ARCHOS),ibm-bg)
  ifndef PARMETIS
    PARMETIS := /usr/local/parmetis/latest
  endif
  ifndef ZOLTAN
    ZOLTAN := /usr/local/zoltan/latest
  endif
 else
  ifndef PARMETIS
    PARMETIS := $(DEVROOT)/parmetis/latest
  endif
  ifndef ZOLTAN
    ZOLTAN := $(DEVROOT)/zoltan
  endif
  endif  
endif

DEPS := $(DEPS) \
	phasta/phShape \
	phasta/phUtil/LU \
	phasta/phastaIO

LIBS := $(LIBS) phShape LU phastaIO

ifeq ($(SIM),1)
  ifeq ($(MODELER),parasolid)
    PARASOLID = $(MESHSIM)/lib/$(SIMARCH)/psKrnl 
    PSKRNL = -lpskernel
    MESHSIMMODELVERS = -lSimParasolid280 
  else
    MODELLIB = -lSim$(MODELER)
  endif 

  INCLUDES := $(INCLUDES) -I$(MESHSIM)/include
  DEFS := $(DEFS) -DSIM -DSIM_$(shell echo $(MODELER)| tr '[a-z]' '[A-Z]')

  MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH)
ifeq ($(MODELER),parasolid) 
  MESHSIMFLAGS := $(MESHSIMFLAGS)  -L$(PARASOLID) 
endif
  MESHSIMFLAGS := $(MESHSIMFLAGS) \
		-lSimPartitionedMesh-mpi -lSimAdvMeshing -lSimMeshing \
		-lSimPartitionWrapper-$(PARALLEL) -lSimPartitionedMesh-mpi $(MESHSIMMODELVERS) -lSimMeshTools -lSimModel\
                $(PSKRNL) 
#		-lSimPartitionWrapper-$(PARALLEL)14 -lSimPartitionedMesh-mpi $(MESHSIMMODELVERS) -lSimMeshTools -lSimModel\
                $(PSKRNL) 
endif

ifeq ($(ENTITY_GROUP),1)
  DEFS := $(DEFS) -DENTITY_GROUP
endif

ifeq ($(FMDB),1)
  ifeq ($(MODELER),meshmodel)
    DEFS := $(DEFS) -DMESHMODEL
    ADDLIB = $(call libtmpl,$(PARMETIS)/lib) -lparmetis -lmetis \
            $(ZOLTAN)/lib/libzoltan.a
  endif
  ifeq ($(MODELER),parasolid)
     DEFS := $(DEFS) -DPARASOLID
     INCLUDES := $(INCLUDES) -I$(PARASOLID)
     DEPS := $(DEPS) $(MCTKVERS)/parasolidModel
     LIBS := $(LIBS) parasolidModel
     PARASOLIDPATH := $(MESHSIM)/lib/$(SIMARCH)/psKrnl
     ADDLIB = $(call libtmpl,$(PARASOLIDPATH)) -lpskernel -lnsl -lm \
              $(call libtmpl,$(PARMETIS)/lib) -lparmetis -lmetis \
              $(ZOLTAN)/lib/libzoltan.a
  endif
  ifeq ($(MODELER),acis)
     DEFS:= $(DEFS) -DACIS
     INCLUDES := $(INCLUDES) -I$(ACIS)/include
     DEPS := $(DEPS) $(MCTKVERS)/acisModel
     LIBS := $(LIBS) acisModel
     ifeq ($(VERS),opt)
       ACISPATH := $(ACIS)/bin/linux_so
     else
       ACISPATH := $(ACIS)/bin/linux_so_debug
     endif
     ADDLIB = $(call libtmpl,$(ACISPATH))  -lSpaAVis  -lSpaBase -lSpaACIS \
              $(call libtmpl,$(PARMETIS)/lib) -lparmetis -lmetis \
              $(ZOLTAN)/lib/libzoltan.a
  endif
  DEPS := $(DEPS) $(FMDBVERS)/SCORECUtil $(FMDBVERS)/SCORECModel $(FMDBVERS)/FMDB 
  LIBS := $(LIBS) SCORECUtil SCORECModel FMDB

  DEPS := $(DEPS) $(MESHADAPTVERS)/meshAdapt \
          $(MESHADAPTVERS)/meshTools
  LIBS := $(LIBS) meshAdapt meshTools
endif

ifeq ($(SIM),1)
LDFLAGS := $(LDFLAGS) $(MESHSIMFLAGS) $(ETCFLAGS) -lnsl -lm -lpthread
else
LDFLAGS := $(LDFLAGS) $(ETCFLAGS) -lnsl -lm
endif

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common

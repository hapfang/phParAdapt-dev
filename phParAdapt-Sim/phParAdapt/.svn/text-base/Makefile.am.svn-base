#bring in the macros defined in the m4 directory.. these include AX_GMI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4


lib_LIBRARIES = libphParAdapt.a
libphParAdapt_a_SOURCES = \
./src/adapt.cc \
./src/allpart.cc \
./src/assignGlobalVars.cc \
./src/attachData.cc \
./src/attachIBC.cc \
./src/attachPeriodicBC.cc \
./src/BCInflowFaceInfo.cc \
./src/BoundaryCondition.cc \
./src/buildSystem.cc \
./src/callback.cc \
./src/commuBCInflowFaceInfo.cc \
./src/commuFix4SolutionTransfer.cc \
./src/commuGradientsFromPatch.cc \
./src/commuHessiansFromPatch.cc \
./src/commuMaxLocalPartLocError.cc \
./src/commuOldMeshSize.cc \
./src/commuSmoothDir.cc \
./src/commuSmoothHessians.cc \
./src/commuSmoothSize.cc \
./src/Eigen.cc \
./src/elementGradient.cc \
./src/elementHessian.cc \
./src/EnsaArrays.cc \
./src/EnsaParameters.cc \
./src/EssentialBC.cc \
./src/face_extraction.cc \
./src/fix4SolutionTransfer.cc \
./src/genblock.cc \
./src/GEN_data.cc \
./src/genILWork.cc \
./src/getData.cc \
./src/getNSHGTOT.cc \
./src/globalInfo.cc \
./src/gradientsFromPatch.cc \
./src/hessiansFromPatch.cc \
./src/hessiansFromSolution.cc \
./src/InitBCinput.cc \
./src/InitialCondition.cc \
./src/Input.cc \
./src/InterpolateSolution.cc \
./src/InverseMap.cc \
./src/lin2quad.cc \
./src/localInfo.cc \
./src/main.cc \
./src/markingStrategy.cc \
./src/maxLocalPartLocError.cc \
./src/mdb2phasta.cc \
./src/M_writeVTKFile.cc \
./src/myAttribute.cc \
./src/NaturalBC.cc \
./src/partitionMesh.cc \
./src/partitionMeshToLBForAdaptivity.cc \
./src/partMeshBypass.cc \
./src/PeriodicBC.cc \
./src/phPartitionCallbacks.cc \
./src/phrase.cc \
./src/phReadWrite.cc \
./src/printBC.cc \
./src/printInfo.cc \
./src/procArgs.cc \
./src/processErrorAG.cc \
./src/processError.cc \
./src/procSize.cc \
./src/readBC.cc \
./src/readErrorFiles.cc \
./src/readErrorFromRestart.cc \
./src/R_edges.cc \
./src/R_entitiesAdapt.cc \
./src/R_entitiesBdry.cc \
./src/R_entities.cc \
./src/reorder.cc \
./src/restart.cc \
./src/return_vertices.cc \
./src/setIsotropicSizeField.cc \
./src/setManualSizeField.cc \
./src/setPeriodic.cc \
./src/setSizeFieldUsingHessians.cc \
./src/setup.cc \
./src/setupGlobalTasks.cc \
./src/SizeLimit.cc \
./src/SmoothErrorIndicators.cc \
./src/SmoothHessians.cc \
./src/SmoothSize.cc \
./src/SmoothStress.cc \
./src/SPEBC.cc \
./src/switchAdapt_Preproc.cc \
./src/tagEntitiesForRefinement.cc \
./src/ThickAdapt.cc \
./src/transformToScalarErrorVal.cc \
./src/UserWeight.cc \
./src/V_Hessian.cc \
./src/Vorticity.cc \
./src/writeEnsaFiles.cc \
./src/writeMEDITHessians.cc \
./src/writeNCorp.cc \
./src/write_restart.cc \
./src/writeRestartFiles.cc \
./src/writeRestartHessians.cc \
./src/bits.c

libphParAdapt_a_include_HEADERS = \
./src/attachData.h \
./src/bits.h \
./src/BoundaryCondition.h \
./src/ccfunc.h \
./src/Eigen.h \
./src/EnsaArrays.h \
./src/EnsaParameters.h \
./src/EssentialBC.h \
./src/func.h \
./src/InitialCondition.h \
./src/Input.h \
./src/mesh_interface.h \
./src/MeshSimInternal.h \
./src/myAttribute.h \
./src/NaturalBC.h \
./src/parallel.h \
./src/PeriodicBC.h \
./src/phParAdapt.h \
./src/phPartitionCallbacks.h \
./src/phReadWrite.h \
./src/Solution.h \
./src/SPEBC.h


libphParAdapt_a_includedir = $(includedir)

AM_CXXFLAGS = -DB64

if PHPARA_PARALLEL
  AM_CXXFLAGS += -DPHPARA_PARALLEL
endif 

AM_CXXFLAGS +=  \
  -DFMDB \
  $(MA_CFLAGS) \
  $(PHASTAIO_INCLUDES) \
  $(PHSHAPE_INCLUDES) \
  -I$(top_srcdir) 

AM_CFLAGS = $(AM_CXXFLAGS)

bin_PROGRAMS = phParAdapt
phParAdapt_SOURCES = ./src/main.cc
phParAdapt_CXXFLAGS = $(HPS_CFLAGS) $(DIFFMC_CFLAGS) 
phParAdapt_LDADD = \
  $(top_builddir)/libphParAdapt.a \
  $(HPS_LIBS) $(DIFFMC_LIBS) \
  $(LU_LIBS) \
  $(PHASTAIO_LIBS) \
  $(PHSHAPE_LIBS) \
  $(MA_LIBS)


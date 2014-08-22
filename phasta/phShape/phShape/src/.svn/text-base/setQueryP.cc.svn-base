/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
            set (query) the polynomial order of interpolation associated 
            with a mesh entity for a given "field"
-------------------------------------------------------------------------*/
#include <stdio.h>
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

int Entity_setP(pEntity entity, pMeshDataId dataID, int p) {
  if( EN_type(entity) == Tvertex ) 
    return 1 ;

  if( EN_getDataInt(entity,dataID,NULL) )
    EN_modifyDataInt(entity,dataID,p);
  else
    EN_attachDataInt(entity,dataID,p);
  return 1 ;
}

int Entity_queryP(pEntity entity, pMeshDataId dataID) {
  int p;
  if( EN_type(entity) == Tvertex )
    return 1;
  else
    if( EN_getDataInt(entity,dataID,&p) )
      return p;
    else
      return 0;
}


#ifdef __cplusplus
}
#endif

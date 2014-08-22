/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             determine the local vertex index of vertices defining the given
             mesh entity within the mesh database.
-------------------------------------------------------------------------*/

#ifdef SIM
#include "MeshSim.h"
pPList R_verticesLeft( pRegion region );
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif
#ifdef __cplusplus
extern "C" {
#ifdef FMDB
pPList R_verticesLeft( pRegion region );
#endif
#endif

int V_index(pVertex v, pEntity ent, int *index) {
  int retval = 1;
  void *temp ;
  int nverts,type = EN_type(ent);
  pVertex vert;
  pPList verts;

  if( type == Tvertex )
    retval = 0;
  else if( type == Tedge ) {
    if( E_vertex((pEdge)ent,0) == v )
      index[0] = 0;
    else if( E_vertex((pEdge)ent,1) == v )
      index[0] = 1;
    else
      retval = 0;
  } else {
    if( type == Tface )
      verts = F_vertices((pFace)ent,1);
    else
      verts = R_vertices((pRegion)ent, 1);

    nverts = PList_size(verts);
    *index = retval = 0 ;
    temp = 0;
    while((vert=(pVertex)PList_next(verts,&temp))) {
      if( vert == v ) {
        retval = 1;
        break;
      }
      else
        (*index)++;
    }
    PList_delete(verts);
  }
  return retval;
}

int E_index(pEdge e, pEntity ent, int *index) {
  pVertex ev,vert;
  pPList verts ;
  void *temp;
  int i,count,retval=1,nverts,type;
  static pEntity lastElement,lastEntity;
  static int lastIndex[2],lastReturn;

  /* if this is the same element entity-element pair dont do the 
     search again 
  */
  if( lastEntity == e && lastElement == ent ) {
    index[0] = lastIndex[0];
    index[1] = lastIndex[1];
    return lastReturn;
  }

  type = EN_type(ent);

  if( type < Tface ) 
    retval = 0;
  else {
    if(type == Tface)
      verts = F_vertices((pFace)ent,1);
    else
      verts = R_vertices((pRegion)ent, 1);
    nverts = PList_size(verts);
    for(i=0; i<2; i++) {
      ev = E_vertex(e,i);
      index[i] = retval = 0;
      temp = 0;
      while((vert=(pVertex)PList_next(verts,&temp))) {
        if( vert == ev) {
          retval = 1;
          break;
	} else
          index[i]++;
      }
      if(!retval)
        break ;
    } 
    PList_delete(verts);
  }
  lastElement = ent;
  lastEntity  = e;
  lastIndex[0] = index[0];
  lastIndex[1] = index[1];
  lastReturn  = retval;
  
  return retval;
}

int F_index(pFace face, pEntity ent, int *index) {
  int type,nverts,count,fcount,retval=1 ;
  pPList fverts,verts;
  pVertex v,fv;
  void *temp,*ftemp;
  static pEntity lastElement,lastEntity;
  static int lastIndex[3],lastReturn;

  /* if this is the same element entity-element pair dont do the 
     search again 
  */
  if( lastEntity == face && lastElement == ent ) {
    index[0] = lastIndex[0];
    index[1] = lastIndex[1];
    index[2] = lastIndex[2];
    return lastReturn;
  }

  type = EN_type(ent);
  if( type == Tregion ) {
    fverts = F_vertices(face,1);
    verts = R_vertices((pRegion)ent, 1);
    nverts = PList_size(verts);
    fcount =0;
    ftemp=0;
    while((fv=(pVertex)PList_next(fverts,&ftemp))) {
      temp = 0;
      count =0;
      while((v=(pVertex)PList_next(verts,&temp))) {
        if( v == fv )
          break;
        else
          count++;
      }
      if( count < nverts )
        index[fcount++] = count;
      else {
        retval = 0;
        break;
      }
    }
    PList_delete(fverts);
    PList_delete(verts);
  } else
    retval = 0;

  lastElement = ent;
  lastEntity  = face;
  lastIndex[0] = index[0];
  lastIndex[1] = index[1];
  lastIndex[2] = index[2];
  lastReturn  = retval;

  return retval ;
}

#ifdef __cplusplus
}
#endif

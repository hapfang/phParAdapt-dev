/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             evaluate the blend funtion for a given mesh entity over another
             mesh entity
-------------------------------------------------------------------------*/

#include <math.h>
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

double V_blendOnEntity(pVertex v, pEntity e, double *L) {
  /* blend a vertex mode on a mesh edge */
  int index ;

  if( V_index(v,e,&index) ) {
    if( EN_type(e) == Tedge )
      return V_blendIndexedOnEdge(index,L);
    else
      return V_blendIndexed(index,L);
  }
  else
    return 0.0;
}

double V_blendIndexed(int i, double *L) {
  return L[i] ;
}

double V_blendIndexedOnEdge(int i, double *L) {
  if( i == 0 )
    return 0.5*(1.0-L[0]);
  else
    return 0.5*(1.0+L[0]);
}

double E_blendOnFace(pEdge edge, pFace face, double *L) {
  /* blend an edge mode on a tri. face */
  int index[2],nedges = F_numEdges(face);

  if(!E_index(edge, face, index) )
    return 0.0;

  /* find out the local edge index in the face */
  if( nedges == 3 ){ 
    return F_edgeBlendTri(index, L);
  } else if( nedges == 4 ){
    return F_edgeBlendQuad(index, L);
  } else
    return 0.0;
}

double F_edgeBlendTri(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_edgeBlendQuad(int *index, double *L) {
  return 0.0;
}

double E_blendOnRegion(pEdge edge, pRegion region, double *L) {
  /* blend a mesh edge mode on a tetra. region */
  int index[2],nfaces= R_numFaces(region);
  
  /* figure our which local edge we are dealing with */
  if(!E_index(edge, region, index) )
    return 0.0;
 if( nfaces == 4 ) 
    return R_edgeBlendTet(index, L);   
 else
  return 0.0;
}

double R_edgeBlendTet(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_blendOnRegion(pFace face, pRegion region, double *L) {
  int index[4] ;

  /* blend a face mode on a tet. region */
  int nfaces = R_numFaces(region);

  if(!F_index(face, region, index)) 
   return 0.0;

  if( nfaces == 4 ) {
    return L[index[0]]*L[index[1]]*L[index[2]] ;
  } else
    return 0.0;
}

#ifdef __cplusplus
}
#endif

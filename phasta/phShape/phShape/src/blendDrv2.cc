/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : May., 97
  Modifi.  : 
  Function :
             evaluate 2nd derivatives of entity blend functions over other 
             entities.
             return array contents:
             Edge: [d2()/dxidxi]
             Face: [d2()/dxidxi,d2()/dxideta,d2()/detadxi,d2()/detadeta]
	     Region:[d2()/d1d1,d2()/d1d2,d2()/d1d3,
	             d2()/d2d1,d2()/d2d2,d2()/d2d3,
                     d2()/d3d1,d2()/d3d2,d2()/d3d3]
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

int V_blendOnEntityDrv2(pVertex v, pEntity e, double *L, double mdrv[]) {
  /* 2nd derivative of blend a vertex mode on a mesh entity 
     since vertex blend is linear for all entities, its zero
  */
  int etype = EN_type(e);
  if(etype == Tedge)
    mdrv[0] = 0.0e0;
  else if(etype == Tface)
    mdrv[0] = mdrv[1] = mdrv[2] = mdrv[3] = 0.0e0;
  else if(etype == Tregion)
    mdrv[0]=mdrv[1]=mdrv[2]=mdrv[3]=mdrv[4]=mdrv[5]=mdrv[6]=mdrv[7]=mdrv[8]=0.0e0;
  return 1;
}

int E_blendOnFaceDrv2(pEdge edge, pFace face, double *L, double bdrv[2]) {
  int index[2], nedges = F_numEdges(face) ;
  
  /* find out the local edge index in the face */
  if( E_index(edge,face,index) ) {
    if( nedges == 3 ) 
      return F_edgeBlendTriDrv2(index, L, bdrv);
  }
  return 0;
}

int F_edgeBlendTriDrv2(int index[2], double *L, double drv[]) {
  double r,s;

  drv[0] = drv[1] = 0.0;
  if( C_equal1(L[0],1.0) || C_equal1(L[1],1.0) || C_equal1(L[2],1.0) )
    return 1;

  r = L[0] ; s = L[1] ;

  /* figure out which edge we are dealing with */
  if( (index[0]==0 && index[1]==1) || (index[0]==1 && index[1]==0) ) {
    /* v0=0, V1=1 */
    drv[0]=drv[3]=0.0e0;
    drv[1]=drv[2]=-2.0e0;
  } else if( (index[0]==1 && index[1]==2) || (index[0]==2 && index[1]==1) ) {
    /* v0=1, V1=2 */
    drv[0]=0.0e0;
    drv[1]=drv[2]=2.0e0;
    drv[3]=4.0e0; 
  } else if( (index[0]==2 && index[1]==0) || (index[0]==0 && index[1]==2) ) {
    /* v0=2, V1=0 */
    drv[0]=4.0e0;
    drv[1]=drv[2]=2.0e0;
    drv[3]=0.0e0; 
  } else
    return 0;
  
  return 1 ;
}

int E_blendOnRegionDrv2(pEdge edge, pRegion region, double *L, double bdrv[3]) {
  int index[2], nfaces=R_numFaces(region);
  
  /* figure our which local edge we are dealing with */
  if( E_index(edge,region,index)) {
    if( nfaces == 4 )
      return R_edgeBlendTetDrv(index, L, bdrv);
  }
   return 0;
}

int F_blendOnRegionDrv2(pFace face, pRegion region, double *L, double drv[3]) {
  /* blend a face mode on a tet. region */
  int nfaces = R_numFaces(region);
  drv[0] = drv[1] = drv[2] = 0.0;
  if( nfaces == 4 )
    return 1;
  else
    return 0;
}

int R_edgeBlendTetDrv2(int *index, double *L, double drv[]) {
  double r=L[0],s=L[1],t=L[2];

  /* the blend is given by -2*L[i]*L[j] */
  if( (index[0]==0 && index[1]==1) || (index[0]==1 && index[1]==0) ) {
    /* v0=0, V1=1 */
    drv[0]=0.0e0;
    drv[1]=-2.0e0;
    drv[2]=0.0e0;
    drv[3]=-2.0e0;
    drv[4]=0.0;
    drv[5]=0.0;
    drv[6]=0.0;
    drv[7]=0.0;
    drv[8]=0.0;
  } else if( (index[0]==0 && index[1]==3) || (index[0]==3 && index[1]==0) ) {
    /* v0=0, V1=3 */
    drv[0] = 4.0e0;
    drv[1] = 2.0e0;
    drv[2] = 2.0e0;
    drv[3]=2.0e0;
    drv[4]=0.0e0;
    drv[5]=0.0e0;
    drv[6]=2.0e0;
    drv[7]=0.0e0;
    drv[8]=0.0e0;
  } else if( (index[0]==1 && index[1]==2) || (index[0]==2 && index[1]==1) ) {
    /* v0=1, V1=2 */
    drv[0]=0.0e0;
    drv[1]=0.0e0;
    drv[2]=0.0e0;
    drv[3]=0.0e0;
    drv[4]=0.0e0;
    drv[5]=-2.0e0;
    drv[6]=0.0e0;
    drv[7]=-2.0e0;
    drv[8]=0.0e0;
  } else if( (index[0]==1 && index[1]==3) || (index[0]==3 && index[1]==1) ) {
    /* v0=1, V1=3 */
    drv[0]=0.0e0;
    drv[1]=2.0e0;
    drv[2]=0.0e0;
    drv[3]=2.0e0;
    drv[4]=4.0e0;
    drv[5]=2.0e0;
    drv[6]=0.0e0;
    drv[7]=2.0e0;
    drv[8]=0.0e0;
  } else if( (index[0]==2 && index[1]==0) || (index[0]==0 && index[1]==2) ) {
    /* v0=2, V1=0 */
    drv[0]=0.0e0;
    drv[1]=0.0e0;
    drv[2]=-2.0e0;
    drv[3]=0.0e0;
    drv[4]=0.0e0;
    drv[5]=0.0e0;
    drv[6]=-2.0e0;
    drv[7]=0.0e0;
    drv[8]=0.0e0;
  } else if( (index[0]==2 && index[1]==3) || (index[0]==3 && index[1]==2) ) {
    /* v0=2, V1=3 */
    drv[0]=0.0e0;
    drv[1]=0.0e0;
    drv[2]=2.0e0;
    drv[3]=0.0e0;
    drv[4]=0.0e0;
    drv[5]=2.0e0;
    drv[6]=2.0e0;
    drv[7]=2.0e0;
    drv[8]=4.0e0;
  } else
    return 0;

  return 1;
}

int R_faceBlendTetDrv2(int *index, double *L, double drv[]) {
  int isum = index[0]+index[1]+index[2];
  double r = L[0];
  double s = L[1];
  double t = L[2];
  if( isum == 3) /* r*s*t */ {
    drv[0]=0.0e0;
    drv[1]=t;
    drv[2]=s;
    drv[3]=t;
    drv[4]=0.0e0;
    drv[5]=r;
    drv[6]=s;
    drv[7]=r;
    drv[8]=0.0e0;
  } else if( isum == 5) /* r*t*(1-r-s-t) */ {
    drv[0]=-2.0e0*t;
    drv[1]=-t;
    drv[2]=1.0e0-2.0e0*(r+t)-s;
    drv[3]=drv[1];
    drv[4]=0.0e0;
    drv[5]=-r;
    drv[6]=drv[2];
    drv[7]=drv[5];
    drv[8]=-2.0e0*r;
  } else if( isum == 4) /* r*s*(1-r-s-t) */ {
    drv[0]=-2.0e0*s;
    drv[1]=1.0e0-2.0e0*(r+s)-t;
    drv[2]=-s;
    drv[3]=drv[1];
    drv[4]=-2.0e0*r;
    drv[5]=-r;
    drv[6]=drv[2];
    drv[7]=drv[5];
    drv[8]=0.0e0;
  } else if( isum == 6) /* s*t*(1-r-s-t) */ {
    drv[0]=0.0e0;
    drv[1]=-t;
    drv[2]=-s;
    drv[3]=drv[1];
    drv[4]=-2.0e0*t;
    drv[5]=1.0e0-2.0e0*(t+s)-r;
    drv[6]=drv[2];
    drv[7]=drv[5];
    drv[8]=-2.0e0*s;
  } else
    return 0;
  return 1;
}
#ifdef __cplusplus
}
#endif

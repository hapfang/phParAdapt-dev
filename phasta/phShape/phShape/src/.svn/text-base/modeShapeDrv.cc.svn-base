/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             return derivative hierarchic mode shapes associated with 
             a given mesh entity of a specified polynomial order.
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

int E_modeShapeDrv(int p, double *L, int order, double drv[]) {
/* Return derivative edge mode shape function evaluated along an edge at 
a point for spectral interpolation order p L[0,1] in [0,1] r+s <= 1
*/
  if(order ==1)
    return EnDrv(p-2,L[0],L[1],drv);
  else if(order ==2)
    return EnDrv2(p-2,L[0],L[1],drv);
  else 
    return 0;
}

int F_modeShapeTriDrv(int p, int i, double *L, int order, double mdrv[]) {
  int alpha,beta,found,count;
  double rs,rst2,P1P2,t,P1,P2,dP1,dP2,t1,t2;
  /* return i-th triangular face mode derivative of polynomial order p 
     note: there are p-2 modes with polynomial order p */

  if( p < 3 || i < 0 || i > p-3 )
    return 0.0 ;

  count = found = 0;
  for(alpha=0; alpha <= p-3; alpha++) {
    for(beta=0; beta <= p-3; beta++) {
      if( alpha+beta == p-3 ) {
        if( count == i )
          found=1;
        else
          count++;
      } 
      if(found)
        break ;   
    }
    if(found)
      break;
  }
  if( found ) {
    if(order==1)
      return FnDrv(alpha,beta,L[0],L[1],mdrv);
    else if(order==2)
      return FnDrv2(alpha,beta,L[0],L[1],mdrv);
    else
      return 0;
  } else
    return 0;
}

int F_modeShapeQuadDrv(int p, int i, double *L, int order, double mdrv[]) {
  return 0;
}

int R_modeShapeTetDrv(int p, int i, double *L, int order, double mdrv[]) {
  int alpha,beta,gamma,count,found ;
  double Pr,Ps,Pt,dPr,dPs,dPt,w,rst,rstw2,PrPsPt,t1;

  /* return the i-th mode shape of polynomial order p for this region , there are
     (p-2)*(p-3)/2 mode shapes of polynomial order p */
  if( p < 4 || i < 0 || i > (((p-2)*(p-3)/2)-1) )
    return 0.0;

  count = 0;
  found = 0;
  for( alpha=0; alpha <= p-4; alpha++) {
    for( beta=0; beta <= p-4; beta++) {
      for( gamma=0; gamma <= p-4; gamma++) {
        if( alpha+beta+gamma == p-4 ) {
          if( count == i )
            found = 1;
          else
            count++;
	}
        if(found) break;
      }
      if(found) break;
    }
    if(found) break;
  }
  if( found ) {
    if(order==1)
      return BnDrv(alpha,beta,gamma,L[0],L[1],L[2],mdrv);
    else if(order==2)
      return BnDrv2(alpha,beta,gamma,L[0],L[1],L[2],mdrv);
    else
      return 0;
  } else
    return 0;
}

int R_modeShapeHexDrv(int p, int i, double *L, int order, double mdrv[]) {
  return 0;
}

#ifdef __cplusplus
}
#endif







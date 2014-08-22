/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             return the derivative of the parametric transformation from
             element to entity coordinate systems.
-------------------------------------------------------------------------*/

#include <math.h>
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif


#ifdef __cplusplus
extern "C" {
#endif

static int SF_par_TABLE[][3][2] = {
   {{0,0},{0,0},{0,0}},
   {{1,0},{0,1},{0,0}},
   {{1,0},{0,0},{0,1}},
   {{1,-1},{0,-1},{0,-1}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,1},{1,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{1,0},{0,1}},
   {{0,-1},{1,-1},{0,-1}}, 
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,1},{0,0},{1,0}},
   {{0,0},{0,1},{1,0}},
   {{0,0},{0,0},{0,0}},
   {{0,-1},{0,-1},{1,-1}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{-1,1},{-1,0},{-1,0}},
   {{-1,0},{-1,1},{-1,0}},
   {{-1,0},{-1,0},{-1,1}}
};

int E_parDrv(int i, int j, pEntity elem, double drv[][2]) {
  int type = EN_type(elem);  
  int retval = 1,ind[2],ii;

  /* find out type of element */
  if( type == Tedge )
    return 0;
  else if( type == Tface ) {
    if( F_numEdges((pFace)elem) == 3) {
      drv[0][0] = drv[0][1] = drv[1][0] = drv[1][1] = 0.0;
      if( i == 0 && j == 1 ) { /* edge +0 , xi = s-r */
        drv[0][0] = 1.0;
        drv[1][1] = 1.0;
      } else if( i == 1 && j == 0 ) { /* edge -0 , xi = r-s */
        drv[0][1] = 1.0;
        drv[1][0] = 1.0;
      } else if( i == 1 && j == 2 ) { /* edge +1 , xi = t-s */
        drv[0][1] = -1.0;
        drv[1][0] = 1.0;
        drv[1][1] = -1.0;
      } else if( i == 2 && j == 1 ) { /* edge -1 , xi = s-t */
        drv[0][0] = -1.0;
        drv[1][0] = -1.0;
        drv[1][1] = 1.0;
      } else if( i == 2 && j == 0 ) { /* edge +2 , xi = r-t */
        drv[0][0] = -1.0;
        drv[0][1] = 1.0;
        drv[1][0] = -1.0;
      } else if( i == 0 && j == 2 ) { /* edge -2 , xi = t-r */
        drv[0][0] = 1.0;
        drv[0][1] = -1.0;
        drv[1][1] = -1.0;
      } else
        retval = 0;
    } else
      retval = 0;
  } else if( type == Tregion ) {
    if( R_numFaces((pRegion)elem) == 4 ) {
      drv[0][0]=drv[0][1]=drv[1][0]=drv[1][1]=drv[2][0]=drv[2][1]=0.0;
      ind[0] = i; ind[1] = j;
      for(ii=0; ii<2; ii++) {
        if(ind[ii] == 0) /* r'= r */
          drv[0][ii] = 1.0;
        else if(ind[ii] == 1) /* r'= s */
          drv[1][ii] = 1.0;
        else if(ind[ii] == 2) /* r'= t */
          drv[2][ii] = 1.0;
        else if(ind[ii] == 3) /* r'= 1-r-s-t */
          drv[0][ii]=drv[1][ii]=drv[2][ii]=-1.0;
      }
    } else
      retval = 0;
  }
  return retval ;
}

int F_parDrv(int i, int j, int k, pEntity elem, int (**drv)[2]) {
  int index = 10*i + j;
  *drv = SF_par_TABLE[index] ;
  return 1;
}

#ifdef __cplusplus
}
#endif

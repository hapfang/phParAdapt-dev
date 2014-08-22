/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : May., 97
  Modifi.  : 
  Function :
             return shape function 2nd derivatives
-------------------------------------------------------------------------*/
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

int Entity_shapeFuncDrv2(pEntity entity, pEntity element, int p, int ith, 
			 double *L, double *shpDrv) { 
  int nent,elType = EN_type(element), enType = EN_type(entity);
  int index[3],(*fpdrv)[2];
  double tmp,mode,blend,modeDrv,bdrv[3],epdrv[3][2],P,dP,tmp1;
  double mdrv2[9],mdrv[3],tmdrv[3],tmdrv2[9],bdrv2[9];
  double rsw,rtw,rst,stw,rstw,Lf[3],rs,rt,st;

  if( enType > elType )
    return 0;

  if( enType == elType ) {
    /* entity and element are same order */
    if( element == entity ) { 
      if( elType == Tedge ) {
        if(!E_modeShapeDrv(p,L,1,mdrv))
          return 0;
        if(!E_modeShapeDrv(p,L,2,mdrv2))
          return 0;
        P = E_modeShape(p,L);
	/* d2/drdr */
        shpDrv[0] = -4.0*L[1]*mdrv[0] -2.0*L[0]*L[1]*mdrv2[0];
	/* d2/drds */
	shpDrv[1] = -2.0*(P+L[1]*mdrv[1]+L[0]*mdrv[0]+L[0]*L[1]*mdrv2[1]);
	/* d2/dsdr */
	shpDrv[2] = shpDrv[1];
	/* d2/dsds */
	shpDrv[3] = -4.0*L[0]*mdrv[0] -2.0*L[0]*L[1]*mdrv2[3];
        return 1 ;
      }
      else if ( elType == Tface ) {
        nent = F_numEdges((pFace)element);
        if( nent == 3 ) {
          if(!F_modeShapeTriDrv(p,ith,L,1,mdrv))
            return 0;
          if(!F_modeShapeTriDrv(p,ith,L,2,mdrv2))
            return 0;
          rs = L[0]*L[1];
	  rt = L[0]*L[2];
	  st = L[1]*L[2];
          rst = rs*L[2];
          mode = F_modeShapeTri(p,ith,L);
	  /* d2/drdr */
          shpDrv[0] = -2.0*L[1]*(mode-L[2]*mdrv[0]) +
	               rs*(L[2]*mdrv2[0]-2.0*mdrv[0]);
	  /* d2/drds */
	  shpDrv[1] = mode*(L[2]-L[1]-L[0]) + mdrv[0]*(rt-rs) + 
	              mdrv[1]*(st-rs) + rst*mdrv2[1];
	  /* d2/dsdr */
	  shpDrv[2] = shpDrv[1];
	  /* d2/dsds */
	  shpDrv[3] = -2.0*L[0]*(mode-L[2]*mdrv[1]) +
	               rs*(L[2]*mdrv2[3]-2.0*mdrv[1]);
          return 1 ;
	}
        else {
          return 0 ;
	}
      }
      else if ( elType == Tregion ) {
        nent = R_numFaces((pRegion)element);
        if( nent == 4 ) {
          if( !R_modeShapeTetDrv(p, ith, L, 1, mdrv) ) 
            return 0;
          if( !R_modeShapeTetDrv(p, ith, L, 2, mdrv2) ) 
            return 0;
          mode = R_modeShapeTet(p,ith,L);
          tmp = 1.0-L[0]-L[1]-L[2] ;
          rsw = L[0]*L[1]*tmp;
          stw = L[1]*L[2]*tmp ;
          rtw = L[0]*L[2]*tmp ;
          rst = L[0]*L[1]*L[2] ;
          rstw = rst*L[3] ;
	  /* d2/drdr */	  
          shpDrv[0] = 2.0*L[1]*L[2]*(L[3]*mdrv[0]-mode) +
	              rst*(L[3]*mdrv2[0]-2.0*mdrv[0]);
	  /* d2/drds */	  
	  shpDrv[1] = mode*(L[2]*L[3]-L[1]*L[2]-L[0]*L[2]) +
	              mdrv[0]*(rtw-rst) + mdrv[1]*(stw-rst) + rstw*mdrv2[1];
	  /* d2/drdt */
	  shpDrv[2] = mode*(L[1]*L[3]-L[1]*L[2]-L[0]*L[1]) +
	              mdrv[0]*(rsw-rst) + mdrv[2]*(stw-rst) + rstw*mdrv2[2];
	  /* d2/dsdr */
	  shpDrv[3] = shpDrv[1];
	  /* d2/dsds */	  
	  shpDrv[4] = 2.0*L[0]*L[2]*(L[3]*mdrv[1]-mode) +
	              rst*(L[3]*mdrv2[4]-2.0*mdrv[1]);
	  /* d2/dsdt */	  
	  shpDrv[5] = mode*(L[0]*L[3]-L[0]*L[2]-L[0]*L[1]) +
	              mdrv[1]*(rsw-rst) + mdrv[2]*(rtw-rst) + rstw*mdrv2[5];
	  /* d2/dtdr */	  
	  shpDrv[6] = shpDrv[2];
	  /* d2/dtds */	  
	  shpDrv[7] = shpDrv[5];
	  /* d2/dtdt */	  
	  shpDrv[8] = 2.0*L[0]*L[1]*(L[3]*mdrv[2]-mode) +
	              rst*(L[3]*mdrv2[8]-2.0*mdrv[2]);
          return 1 ;
	}
        else if ( nent == 6 )
          return 0;
        else 
          return 0;
      }
    } else
      return 0;
  } else {
    /* entity and element are NOT same order, need blending */
    if( enType == Tvertex )
      return V_blendOnEntityDrv2((pVertex)entity, element, L, shpDrv);
    else if( enType == Tedge ) {
      if(!E_index((pEdge)entity,element,index) )
        return 0;
      /* get the mode shape */
      Lf[0] = L[index[0]];
      Lf[1] = L[index[1]];
      mode = E_modeShape(p,Lf);
    
      /* get the mode shape derivative */
      if(!E_modeShapeDrv(p,Lf,1,mdrv))
        return 0;
      if(!E_modeShapeDrv(p,Lf,2,mdrv2))
	return 0;
      
      /* now get the derivative of the edge parameter xi wrt the 
         element coordinate and transform mdrv and mdrv2*/
      if(!E_parDrv(index[0],index[1],element,epdrv))
        return 0;

      if( elType == Tface ) {
	tmdrv[0] = mdrv[0]*epdrv[0][0]+mdrv[1]*epdrv[0][1];
	tmdrv[1] = mdrv[0]*epdrv[1][0]+mdrv[1]*epdrv[1][1];
	
	tmdrv2[0] = mdrv2[0]*epdrv[0][0]*epdrv[0][0] +
	            mdrv2[1]*epdrv[0][0]*epdrv[0][1] +
	            mdrv2[2]*epdrv[0][1]*epdrv[0][0] +
	            mdrv2[3]*epdrv[0][1]*epdrv[0][1] ;
	
	tmdrv2[1] = mdrv2[0]*epdrv[0][0]*epdrv[1][0] +
	            mdrv2[1]*epdrv[0][0]*epdrv[1][1] +
	            mdrv2[2]*epdrv[0][1]*epdrv[1][0] +
	            mdrv2[3]*epdrv[0][1]*epdrv[1][1] ;
	
	tmdrv2[2] = mdrv2[0]*epdrv[1][0]*epdrv[0][0] +
	            mdrv2[1]*epdrv[1][0]*epdrv[0][1] +
	            mdrv2[2]*epdrv[1][1]*epdrv[0][0] +
	            mdrv2[3]*epdrv[1][1]*epdrv[0][1] ;
	
	tmdrv2[3] = mdrv2[0]*epdrv[1][0]*epdrv[1][0] +
	            mdrv2[1]*epdrv[1][0]*epdrv[1][1] +
	            mdrv2[2]*epdrv[1][1]*epdrv[1][0] +
	            mdrv2[3]*epdrv[1][1]*epdrv[1][1] ;
	
        /* get the blend, its derivative */
        if( F_numEdges((pFace)element) == 3 ) {
          blend = F_edgeBlendTri(index,L);
          if(!F_edgeBlendTriDrv(index,L,bdrv))
            return 0;
          if(!F_edgeBlendTriDrv2(index,L,bdrv2))
            return 0;
	  /* d2/drdr */
          shpDrv[0] = blend*mdrv2[0] + 2.0e0*bdrv[0]*tmdrv[0] + mode*tmdrv2[0];
	  /* d2/drds */
	  shpDrv[1] = blend*mdrv2[1] + tmdrv[0]*bdrv[1] + tmdrv[1]*bdrv[0] +
	              mode*tmdrv2[1];
	  /* d2/dsdr */
	  shpDrv[2] = blend*mdrv2[2] + tmdrv[0]*bdrv[1] + tmdrv[1]*bdrv[0] +
	              mode*tmdrv2[2];
	  /* d2/dsds */
	  shpDrv[3] = blend*mdrv2[3] + 2.0e0*bdrv[1]*tmdrv[1] +
	              mode*tmdrv2[3];
          return 1 ;
	}
        return 0;
      } else if ( elType == Tregion ) {
	tmdrv[0] = mdrv[0]*epdrv[0][0]+mdrv[1]*epdrv[0][1];
	tmdrv[1] = mdrv[0]*epdrv[1][0]+mdrv[1]*epdrv[1][1];
	tmdrv[2] = mdrv[0]*epdrv[2][0]+mdrv[1]*epdrv[2][1];

	tmdrv2[0] = mdrv2[0]*epdrv[0][0]*epdrv[0][0] +
	            mdrv2[1]*epdrv[0][0]*epdrv[0][1] +
	            mdrv2[2]*epdrv[0][1]*epdrv[0][0] +
	            mdrv2[3]*epdrv[0][1]*epdrv[0][1] ;

	tmdrv2[1] = mdrv2[0]*epdrv[0][0]*epdrv[1][0] +
	            mdrv2[1]*epdrv[0][0]*epdrv[1][1] +
	            mdrv2[2]*epdrv[0][1]*epdrv[1][0] +
	            mdrv2[3]*epdrv[0][1]*epdrv[1][1] ;

	tmdrv2[2] = mdrv2[0]*epdrv[0][0]*epdrv[2][0] +
	            mdrv2[1]*epdrv[0][0]*epdrv[2][1] +
	            mdrv2[2]*epdrv[0][1]*epdrv[2][0] +
	            mdrv2[3]*epdrv[0][1]*epdrv[2][1] ;
	
	tmdrv2[3] = mdrv2[0]*epdrv[1][0]*epdrv[0][0] +
	            mdrv2[1]*epdrv[1][0]*epdrv[0][1] +
	            mdrv2[2]*epdrv[1][1]*epdrv[0][0] +
	            mdrv2[3]*epdrv[1][1]*epdrv[0][1] ;
	
	tmdrv2[4] = mdrv2[0]*epdrv[1][0]*epdrv[1][0] +
	            mdrv2[1]*epdrv[1][0]*epdrv[1][1] +
	            mdrv2[2]*epdrv[1][1]*epdrv[1][0] +
	            mdrv2[3]*epdrv[1][1]*epdrv[1][1] ;
	
	tmdrv2[5] = mdrv2[0]*epdrv[1][0]*epdrv[2][0] +
	            mdrv2[1]*epdrv[1][0]*epdrv[2][1] +
	            mdrv2[2]*epdrv[1][1]*epdrv[2][0] +
	            mdrv2[3]*epdrv[1][1]*epdrv[2][1] ;
	
	tmdrv2[6] = mdrv2[0]*epdrv[2][0]*epdrv[0][0] +
	            mdrv2[1]*epdrv[2][0]*epdrv[0][1] +
	            mdrv2[2]*epdrv[2][1]*epdrv[0][0] +
	            mdrv2[3]*epdrv[2][1]*epdrv[0][1] ;
	
	tmdrv2[7] = mdrv2[0]*epdrv[2][0]*epdrv[1][0] +
	            mdrv2[1]*epdrv[2][0]*epdrv[1][1] +
	            mdrv2[2]*epdrv[2][1]*epdrv[1][0] +
	            mdrv2[3]*epdrv[2][1]*epdrv[1][1] ;

	tmdrv2[8] = mdrv2[0]*epdrv[2][0]*epdrv[2][0] +
	            mdrv2[1]*epdrv[2][0]*epdrv[2][1] +
	            mdrv2[2]*epdrv[2][1]*epdrv[2][0] +
	            mdrv2[3]*epdrv[2][1]*epdrv[2][1] ;
		
        if( R_numFaces((pRegion)element) == 4 ) {
          blend = R_edgeBlendTet(index,L);
          if(!R_edgeBlendTetDrv(index,L,bdrv))
            return 0;
          if(!R_edgeBlendTetDrv2(index,L,bdrv2))
            return 0;
	  /* d2/drdr */
          shpDrv[0] = blend*tmdrv2[0] + 2.0e0*bdrv[0]*tmdrv[0] + mode*bdrv2[0];

	  /* d2/drds */
          shpDrv[1] = blend*tmdrv2[1] + tmdrv[0]*bdrv[1] + tmdrv[1]*bdrv[0] +
	              mode*bdrv2[1];

	  /* d2/drdt */
          shpDrv[2] = blend*tmdrv2[2] + tmdrv[0]*bdrv[2] + tmdrv[2]*bdrv[0] +
	              mode*bdrv2[2];

	  /* d2/dsdr */
          shpDrv[3] = blend*tmdrv2[3] + tmdrv[1]*bdrv[0] + tmdrv[0]*bdrv[1] +
	              mode*bdrv2[3];
	              
	  /* d2/dsds */
          shpDrv[4] = blend*tmdrv2[4] + 2.0e0*bdrv[1]*tmdrv[1] + mode*bdrv2[4];

	  /* d2/dsdt */
          shpDrv[5] = blend*tmdrv2[5] + tmdrv[1]*bdrv[2] + tmdrv[2]*bdrv[1] +
	              mode*bdrv2[5];

	  /* d2/dtdr */
          shpDrv[6] = blend*tmdrv2[6] + tmdrv[2]*bdrv[0] + tmdrv[0]*bdrv[2] +
	              mode*bdrv2[6];

	  /* d2/dtds */
          shpDrv[7] = blend*tmdrv2[7] + tmdrv[2]*bdrv[1] + tmdrv[1]*bdrv[2] +
	              mode*bdrv2[7];

	  /* d2/dtdt */
          shpDrv[8] = blend*tmdrv2[8] + 2.0e0*bdrv[2]*tmdrv[2] + mode*bdrv2[8];

          return 1 ;
	}
        return 0;
      }
    } else if ( enType == Tface ) {
      /* only works for triangular faces of a te. right now */
      if( F_numEdges((pFace)entity) != 3  || R_numFaces((pRegion)element) != 4 ) 
        return 0;  
      /* determine the local indices of vertices of the face in this "element" */
      if( !F_index((pFace)entity,element,index) )
        return 0.0;

      Lf[0] = L[index[0]];
      Lf[1] = L[index[1]];
      Lf[2] = L[index[2]];
    
      /* get the mode shape of the face, and its derivatives */
      mode = F_modeShapeTri(p,ith,Lf);
      if(!F_modeShapeTriDrv(p,ith,Lf,1,mdrv))
        return 0;

      if(!F_modeShapeTriDrv(p,ith,Lf,2,mdrv2))
        return 0;

      if(!F_parDrv(index[0],index[1],index[2],element,&fpdrv))
        return 0;

      tmdrv[0] = mdrv[0]*(double)fpdrv[0][0] + mdrv[1]*(double)fpdrv[0][1] ;
      tmdrv[1] = mdrv[0]*(double)fpdrv[1][0] + mdrv[1]*(double)fpdrv[1][1] ;
      tmdrv[2] = mdrv[0]*(double)fpdrv[2][0] + mdrv[1]*(double)fpdrv[2][1] ;

      tmdrv2[0] = mdrv2[0]*fpdrv[0][0]*fpdrv[0][0] +
	mdrv2[1]*fpdrv[0][0]*fpdrv[0][1] +
	mdrv2[2]*fpdrv[0][1]*fpdrv[0][0] +
	mdrv2[3]*fpdrv[0][1]*fpdrv[0][1] ;
      
      tmdrv2[1] = mdrv2[0]*fpdrv[0][0]*fpdrv[1][0] +
	mdrv2[1]*fpdrv[0][0]*fpdrv[1][1] +
	mdrv2[2]*fpdrv[0][1]*fpdrv[1][0] +
	mdrv2[3]*fpdrv[0][1]*fpdrv[1][1] ;
      
      tmdrv2[2] = mdrv2[0]*fpdrv[0][0]*fpdrv[2][0] +
	mdrv2[1]*fpdrv[0][0]*fpdrv[2][1] +
	mdrv2[2]*fpdrv[0][1]*fpdrv[2][0] +
	mdrv2[3]*fpdrv[0][1]*fpdrv[2][1] ;
      
      tmdrv2[3] = mdrv2[0]*fpdrv[1][0]*fpdrv[0][0] +
	mdrv2[1]*fpdrv[1][0]*fpdrv[0][1] +
	mdrv2[2]*fpdrv[1][1]*fpdrv[0][0] +
	mdrv2[3]*fpdrv[1][1]*fpdrv[0][1] ;
      
      tmdrv2[4] = mdrv2[0]*fpdrv[1][0]*fpdrv[1][0] +
	mdrv2[1]*fpdrv[1][0]*fpdrv[1][1] +
	mdrv2[2]*fpdrv[1][1]*fpdrv[1][0] +
	mdrv2[3]*fpdrv[1][1]*fpdrv[1][1] ;
      
      tmdrv2[5] = mdrv2[0]*fpdrv[1][0]*fpdrv[2][0] +
	mdrv2[1]*fpdrv[1][0]*fpdrv[2][1] +
	mdrv2[2]*fpdrv[1][1]*fpdrv[2][0] +
	mdrv2[3]*fpdrv[1][1]*fpdrv[2][1] ;
      
      tmdrv2[6] = mdrv2[0]*fpdrv[2][0]*fpdrv[0][0] +
	mdrv2[1]*fpdrv[2][0]*fpdrv[0][1] +
	mdrv2[2]*fpdrv[2][1]*fpdrv[0][0] +
	mdrv2[3]*fpdrv[2][1]*fpdrv[0][1] ;
      
      tmdrv2[7] = mdrv2[0]*fpdrv[2][0]*fpdrv[1][0] +
	mdrv2[1]*fpdrv[2][0]*fpdrv[1][1] +
	mdrv2[2]*fpdrv[2][1]*fpdrv[1][0] +
	mdrv2[3]*fpdrv[2][1]*fpdrv[1][1] ;
      
      tmdrv2[8] = mdrv2[0]*fpdrv[2][0]*fpdrv[2][0] +
	mdrv2[1]*fpdrv[2][0]*fpdrv[2][1] +
	mdrv2[2]*fpdrv[2][1]*fpdrv[2][0] +
	mdrv2[3]*fpdrv[2][1]*fpdrv[2][1] ;
		
      /* get the blend and its derivative */
      blend = Lf[0]*Lf[1]*Lf[2];
      if(!R_faceBlendTetDrv(index,L,bdrv))
        return 0;
      if(!R_faceBlendTetDrv2(index,L,bdrv2))
        return 0;
      /* d2/drdr */
      shpDrv[0] = blend*tmdrv2[0] + 2.0e0*bdrv[0]*tmdrv[0] + mode*bdrv2[0];
      
      /* d2/drds */
      shpDrv[1] = blend*tmdrv2[1] + tmdrv[0]*bdrv[1] + tmdrv[1]*bdrv[0] +
	mode*bdrv2[1];

      /* d2/drdt */
      shpDrv[2] = blend*tmdrv2[2] + tmdrv[0]*bdrv[2] + tmdrv[2]*bdrv[0] +
	mode*bdrv2[2];
      
      /* d2/dsdr */
      shpDrv[3] = blend*tmdrv2[3] + tmdrv[1]*bdrv[0] + tmdrv[0]*bdrv[1] +
	mode*bdrv2[3];
      
      /* d2/dsds */
      shpDrv[4] = blend*tmdrv2[4] + 2.0e0*bdrv[1]*tmdrv[1] + mode*bdrv2[4];
      
      /* d2/dsdt */
      shpDrv[5] = blend*tmdrv2[5] + tmdrv[1]*bdrv[2] + tmdrv[2]*bdrv[1] +
	mode*bdrv2[5];
      
      /* d2/dtdr */
      shpDrv[6] = blend*tmdrv2[6] + tmdrv[2]*bdrv[0] + tmdrv[0]*bdrv[2] +
	mode*bdrv2[6];
      
      /* d2/dtds */
      shpDrv[7] = blend*tmdrv2[7] + tmdrv[2]*bdrv[1] + tmdrv[1]*bdrv[2] +
	mode*bdrv2[7];
      
      /* d2/dtdt */
      shpDrv[8] = blend*tmdrv2[8] + 2.0e0*bdrv[2]*tmdrv[2] + mode*bdrv2[8];
      
      return 1;     
    }
  }
}

#ifdef __cplusplus
}
#endif


/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             return shape function (derivatives) and other related info
             on an entity based query.
-------------------------------------------------------------------------*/

#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

#include "shapeFuncInternals.h"
#include "shapeFunction.h"

#ifdef __cplusplus
extern "C" {
#endif

static int H_tetEmap[][2] = {{0,2},{1,2},{0,1},{0,3},{2,3},{1,3}};

int Entity_shapeMinPolyOrder(pEntity entity) {
  switch( EN_type(entity) ) {
    case Tvertex :
      return 1 ;
    case Tedge :
      return 2 ;
    case Tface :
      return F_numEdges((pFace)entity);
    case Tregion :
      return R_numFaces((pRegion)entity);
    default : /* this is an error */
      return 0;
  }
}

int Entity_numShapeFunction(pEntity entity, int porder) {
  int nfaces,nedges;

  switch( EN_type(entity) ) {
    case Tvertex :
      return ( porder == 1 ? 1 : 0 ) ;
    case Tedge :
      return ( porder >= 2 ? 1 : 0 ) ;
    case Tface :
      nedges = F_numEdges((pFace)entity);
      if( nedges == 3 ) { 
        if( porder >= 3 )
          return (porder-2) ;        
        else
          return 0;
      } else if ( nedges == 4 ) {
        if( porder >= 4 )
          return (porder-3) ;        
        else
          return 0;
      } else
        return 0;
    case Tregion :
      nfaces = R_numFaces((pRegion)entity);
      if( nfaces == 4 ) {
        if( porder >= 4 ) 
          return (porder-3)*(porder-2)/2 ;
        else 
          return 0;
      } else if ( nfaces == 6 ) {
        if( porder >= 6 ) 
          return (porder-5)*(porder-4)/2 ;
        else 
          return 0;
      } else
        return 0;
    default :
      return 0;
  }
}

double Entity_shapeFunction(pEntity entity, pEntity element, int p, int ith, 
                            double *L) {
  int entType = EN_type(entity);
  int index[4], elType  = EN_type(element),nedges,nfaces;
  double Lf[4],modeShape,blend ;

  if( elType < entType )
    return 0.0;
  else if( entity == element ) {
    if( entType == Tedge ) {
      return -2.0*L[0]*L[1]*E_modeShape(p,L);
    }
    else if( entType == Tface ) {
      if( F_numEdges((pFace)entity) == 3 ) 
        return L[0]*L[1]*L[2]*F_modeShapeTri(p,ith,L);
      else
        return F_modeShapeQuad(p,ith,L);
    } else if( entType == Tregion ) {
      if( R_numFaces((pRegion)entity) == 4 )
        return L[0]*L[1]*L[2]*L[3]*R_modeShapeTet(p,ith,L);
      else if( R_numFaces((pRegion)entity) == 6 )
        return R_modeShapeHex(p,ith,L);
    }
  }

  /* do based on type of entity */
  if( entType == Tvertex ) /* vertex mode */
    return V_blendOnEntity((pVertex)entity, element, L) ;
  else if ( entType == Tedge ) { /* edge mode */
    /* determine the local index of the edge in this "element" */
    if( !E_index((pEdge)entity,element,index) )
      return 0.0;

    /* get the edge-->element blend */
    if( elType == Tface ) {
      nedges = F_numEdges((pFace)element) ;
      if( nedges == 3 )
        blend = F_edgeBlendTri(index,L) ;
      else if( nedges == 4 )
        blend = F_edgeBlendQuad(index,L) ;
    } else if ( elType == Tregion ) {
      nfaces = R_numFaces((pRegion)element) ;
      if( nfaces == 4 )
        blend = R_edgeBlendTet(index,L) ;
    }
    Lf[0] = L[index[0]] ;
    Lf[1] = L[index[1]] ;
    /* get the shape function */
    return blend*E_modeShape(p,Lf) ;
  } else if ( entType == Tface ) {
    /* only works for triangular faces of a tet. right now */
    if( F_numEdges((pFace)entity) != 3  || R_numFaces((pRegion)element) != 4 ) 
      return 0.0;  

    /* determine the local indices of vertices of the face in this "element" */
    if( !F_index((pFace)entity,element,index) )
      return 0.0;

    Lf[0] = L[index[0]];
    Lf[1] = L[index[1]];
    Lf[2] = L[index[2]];
    
    /* get the mode shape of the face */
    modeShape = Lf[0]*Lf[1]*Lf[2]*F_modeShapeTri(p,ith,Lf);
    
    return modeShape ;
  } else
    return 0.0;
}

int Entity_shapeFuncDrv(pEntity entity, pEntity element, int p, int ith, 
                            double *L, double *shpDrv) { 
  int nent,elType = EN_type(element), enType = EN_type(entity),index[3],(*fpdrv)[2];
  double tmp,mode,blend,modeDrv,bdrv[3],epdrv[3][2],P,dP,tmp1,mfdrv[2],mdrv[3];
  double rsw,rtw,rst,stw,modeShape,rstw,Lf[3];

  if( enType > elType )
    return 0;


  if( enType == elType ) {
    /* entity and element are same order */
    if( element == entity ) { 
      if( elType == Tedge ) {
        if(!E_modeShapeDrv(p,L,1,mdrv))
          return 0;
        P = E_modeShape(p,L);
        shpDrv[0] = -2.0*L[1]*(P+L[0]*mdrv[0]);
        shpDrv[1] = -2.0*L[0]*(P+L[1]*mdrv[1]);
        return 1 ;
      }
      else if ( elType == Tface ) {
        nent = F_numEdges((pFace)element);
        if( nent == 3 ) {
          if(!F_modeShapeTriDrv(p,ith,L,1,mdrv))
            return 0;
          tmp = L[0]*L[1] ;
          tmp1 = tmp*L[2] ;
          mode = F_modeShapeTri(p,ith,L);
          shpDrv[0] = mdrv[0]*tmp1 + mode*(L[1]*L[2]-tmp);
          shpDrv[1] = mdrv[1]*tmp1 + mode*(L[0]*L[2]-tmp);
          return 1 ;
	}
        else {
          return 0 ;
	}
      }
      else if ( elType == Tregion ) {
        nent = R_numFaces((pRegion)element);
        if( nent == 4 ) {
          if( !R_modeShapeTetDrv(p, ith, L, 1,mdrv) ) 
            return 0;
          mode = R_modeShapeTet(p,ith,L);
          tmp = 1.0-L[0]-L[1]-L[2] ;
          rsw = L[0]*L[1]*tmp;
          stw = L[1]*L[2]*tmp ;
          rtw = L[0]*L[2]*tmp ;
          rst = L[0]*L[1]*L[2] ;
          rstw = rst*L[3] ;
          shpDrv[0] = mode*(stw-rst) + rstw*mdrv[0] ;
          shpDrv[1] = mode*(rtw-rst) + rstw*mdrv[1] ;
          shpDrv[2] = mode*(rsw-rst) + rstw*mdrv[2] ;
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
      return V_blendOnEntityDrv((pVertex)entity, element, L, shpDrv);
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
      
      /* now get the derivative of the edge parameter xi wrt the 
         element coordinate */
      if(!E_parDrv(index[0],index[1],element,epdrv))
        return 0;

      if( elType == Tface ) {
        /* get the blend, its derivative */
        if( F_numEdges((pFace)element) == 3 ) {
          blend = F_edgeBlendTri(index,L);
          if(!F_edgeBlendTriDrv(index,L,bdrv))
            return 0;
          shpDrv[0] = bdrv[0]*mode + 
                      blend*(mdrv[0]*epdrv[0][0]+mdrv[1]*epdrv[0][1]);
          shpDrv[1] = bdrv[1]*mode + 
                      blend*(mdrv[0]*epdrv[1][0]+mdrv[1]*epdrv[1][1]);
          return 1 ;
	}
        return 0;
      } else if ( elType == Tregion ) {
        if( R_numFaces((pRegion)element) == 4 ) {
          blend = R_edgeBlendTet(index,L);
          if(!R_edgeBlendTetDrv(index,L,bdrv))
            return 0;
          shpDrv[0] = bdrv[0]*mode + 
                      blend*(mdrv[0]*epdrv[0][0]+mdrv[1]*epdrv[0][1]);
          shpDrv[1] = bdrv[1]*mode + 
                      blend*(mdrv[0]*epdrv[1][0]+mdrv[1]*epdrv[1][1]);
          shpDrv[2] = bdrv[2]*mode + 
                      blend*(mdrv[0]*epdrv[2][0]+mdrv[1]*epdrv[2][1]);
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
      modeShape = F_modeShapeTri(p,ith,Lf);
      if(!F_modeShapeTriDrv(p,ith,Lf,1,mfdrv))
        return 0;

      if(!F_parDrv(index[0],index[1],index[2],element,&fpdrv))
        return 0;

      mdrv[0] = mfdrv[0]*(double)fpdrv[0][0] + mfdrv[1]*(double)fpdrv[0][1] ;
      mdrv[1] = mfdrv[0]*(double)fpdrv[1][0] + mfdrv[1]*(double)fpdrv[1][1] ;
      mdrv[2] = mfdrv[0]*(double)fpdrv[2][0] + mfdrv[1]*(double)fpdrv[2][1] ;

      /* get the blend and its derivative */
      blend = Lf[0]*Lf[1]*Lf[2];
      if(!R_faceBlendTetDrv(index,L,bdrv))
        return 0;
      
      shpDrv[0] = bdrv[0]*modeShape + blend*mdrv[0] ;      
      shpDrv[1] = bdrv[1]*modeShape + blend*mdrv[1] ;      
      shpDrv[2] = bdrv[2]*modeShape + blend*mdrv[2] ;
      return 0;     
    }
  }
}

#ifdef __cplusplus
}
#endif


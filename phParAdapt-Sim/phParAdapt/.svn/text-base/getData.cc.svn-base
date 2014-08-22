/*

This file contains functions which collect information from a mesh
into arrays needed by PHASTA

*/

#include <stdio.h>
#include<stdlib.h>
#include <map>
#include <iostream>
// #ifndef SIM
// #include "MSopsInternal.h"
// using namespace SCOREC_mesh;
// #else
#include "MeshSimInternal.h"
// #endif
#include "func.h"
#include "parallel.h"


extern forwardblock Iblock;
extern forwardblock Bblock;
extern "C" int topology(pRegion rgn);

extern pMeshDataId MYCTID;
extern pMeshDataId NDOFID;
extern pMeshDataId POLYID;
extern pMeshDataId RNENID;
extern pMeshDataId ReorderR;
extern int isReorder;

static int modeSign(int type, int pdof);
static int Fsign(pFace face,pVertex v[3]);
static int h_modeSign(int type, int pdof);

/**********************************************************************/
/* compute the coordinate array                                       */
/**********************************************************************/
void getX(pMesh mesh, double **x)
{
  pVertex vertex;
  int dof=0;
  void *temp=0;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) { 
      EN_getDataInt(vertex,MYCTID,&dof);
      V_coord(vertex, x[dof]);
  }
  VIter_delete(vIter);
}

/**********************************************************************/
/* compute the interior and boundary element connectivity             */
/**********************************************************************/
void getConnectivity(pMesh mesh,pPList bdry, int ***ien, int ***ienb, globalInfo* info)
{
  /* Connectivity Templates Begin */

  int TetEMap[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
  int TetFMap[4][4]={{0,1,2,-1},{0,3,1,-1},{1,3,2,-1},{0,2,3,-1}};


  int HexEMap[12][2]={{0,1},{1,2},{3,2},{0,3},{4,5},{5,6},{7,6},{4,7},
                      {0,4},{1,5},{2,6},{3,7}};
  int HexFMap[6][4]={{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},
                     {4,5,6,7}};
  int HexfaceDir[6]={0,0,1,0,0,1};

  int WedEMap[9][2]={{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}};
  int WedFMap[5][4]={{0,2,1,-1},{0,1,4,3},{1,2,5,4},{2,0,3,5},{3,4,5,-1}};
  int WedfaceDir[6]={0,0,0,0,0,-1};

  int PyrEMap[8][2]={{0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4}};
  int PyrFMap[5][4]={{0,1,2,3},{1,0,4,-1},{2,1,4,-1},{3,2,4,-1},{0,3,4,-1}};
  int PyrfaceDir[6]={1,1,1,1,1,-1};

  /* Connectivity Templates End */


  int** GenericEdgeMap;
  int** GenericFaceMap;
  int*  GenericfaceDir;

  pVertex *Rvertices;
  pEdge   *Redges;
  pFace   *Rfaces;

  pRegion region;
  pPList ents;
  pVertex v1,v[3];
  pFace face;
  int i,j,nem,nrm;
  int ldof,esign,fsign,msign,k,nsh,is,ip;
  int PolyOrd;
  int nen,nfaces, nedges;
  int nvf ; /* number of vertices on this face */
  int blockid, poly;
  int myctint;
  blockKey BLOCK;
  /* data structures for counting the interior and boundary elements */
  std::map<int,int> iel;
  std::map<int,int> ielb;

  /******************** interior elements ********************/

  GenericEdgeMap = new int* [ 12 ];
  for(i=0; i< 12; i++) GenericEdgeMap[i] = new int [ 2 ];

  GenericFaceMap = new int* [ 6 ];
  for(i=0; i< 6; i++) GenericFaceMap[i] = new int [ 4 ];

  GenericfaceDir = new int [ 6 ];

  /* Using the maximum values possible for verts, edges and faces */

  Rvertices = new pVertex [ 8 ];
  Redges = new pEdge [ 12 ];
  Rfaces = new pFace [ 6 ];

  RIter rIter = M_regionIter(mesh);
  while (region = RIter_next(rIter)){
    /* processor for this region */

    EN_getDataInt((pEntity)region, POLYID, &poly);
    EN_getDataInt((pEntity)region, RNENID, &nen);
    BLOCK.nen = nen;
    BLOCK.maxpoly = poly;
    BLOCK.nenbl = nen == 8 ? 4 : 3;
    BLOCK.lcsyst = topology(region);
    blockid = Iblock[BLOCK]-1;

    nedges = nen+(nen+1)/2;
    nfaces = nen == 8 ? 6: (nen-(nen/6));

    /* R_entites and R_entitiesBdry now passed 0 for edges and faces if linear
     * these fns wont waste time finding higher order ents -- SST */

    /* get all mesh entities connected to this region */
    R_entities(region, Rvertices, info->edgeson?Redges:0, info->faceson?Rfaces:0,nen);
    
    if(isReorder){
        int Rlable;
        EN_getDataInt(region,ReorderR,&Rlable);
        iel[blockid] = Rlable;
    }

    /* collect vertex mode numbers */
    for (i=0; i < nen; i++) {
      EN_getDataInt((pEntity) Rvertices[i], MYCTID, &myctint);
      ien[blockid][iel[blockid]][i] = myctint;
    }

    ldof = nen;

    /* Before we start on the higher order modes we should point the
       topology connectivity templates to the right arrays based on
       the nen (uniquely defines topology) */

    if(info->edgeson) {

      switch(nen){
      case 4:   /* Tetrahedra */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = TetEMap[i][0];
          GenericEdgeMap[i][1] = TetEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = TetFMap[i][j];

        break;
      case 5:   /* Pyramid */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = PyrEMap[i][0];
          GenericEdgeMap[i][1] = PyrEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = PyrFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = PyrfaceDir[i];

        break;

      case 6:  /* Wedge Element */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = WedEMap[i][0];
          GenericEdgeMap[i][1] = WedEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = WedFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = WedfaceDir[i];

        break;

      case 8:  /* Hexahedron */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = HexEMap[i][0];
          GenericEdgeMap[i][1] = HexEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = HexFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = HexfaceDir[i];

        break;
      }

      /* collect edge mode numbers */
      /* find the direction the edge is being used by the region.
         If the first vertex in the edge map for this edge coincides
         with the first vertex of the edge (as it was defined),
         then the sign is 1, otherwise it is 0 */

      nem = poly-1;

      for (i=0; i < nedges; i++){
        v1 = Rvertices[GenericEdgeMap[i][0]];
        esign = (v1 == E_vertex(Redges[i],0));
        for (j=0; j < nem; j++){
          msign = h_modeSign(1,j+2);
          EN_getDataInt((pEntity) Redges[i], MYCTID, &myctint);
          ien[blockid][iel[blockid]][ldof++] = myctint + j;
          if((esign == 0) &&(msign ==0))
            ien[blockid][iel[blockid]][ldof-1] *= -1;
        }
      }
    }

    /* collecting  face mode numbers */
    /* Note : Code Reuse

       Reagardless of element topology only 2 types of faces exist ,
       tri/quad so we can write just 2 cases and use them for all
       topologies. Ofcourse the connectivity template is different for
       each template but that is taken care of using the switch
       statement above .

       In any case it is useless to add edge+triface+quadface
       specifically for each topology we encounter */


    if(info->faceson){

      for(i =0; i < nfaces; i++){

        j =0;
        nvf = F_numEdges(Rfaces[i]);
        PolyOrd =poly;

        /* Triangular Face */
        if (nvf == 3){
          /* find the direction the face is being used */
          for (k=0; k < 3; k++)
            v[k] = Rvertices[GenericFaceMap[i][k]];
          fsign = Fsign(Rfaces[i],v);
          for (ip = 3; ip <= PolyOrd; ip++){
            msign = modeSign(2,ip);
            /* get the number of shape functions of each order */
            nsh = ip-2;
            for (is=0; is < nsh; is++){
              EN_getDataInt((pEntity) Rfaces[i], MYCTID, &myctint);
              ien[blockid][iel[blockid]][ldof++] = myctint + j++;
              if ((fsign == 0) && (msign == 0))
                ien[blockid][iel[blockid]][ldof-1] *= -1;
            }
          }
          /* Quadrilateral Face */
        } else if (nvf == 4) {
          /* find the direction in which the face is being used */
          fsign = R_dirUsingFace(region,Rfaces[0]);
          /* A zero here usually means that the face is defined pointing
             into this particualar region */
          for(ip = 4; ip <= PolyOrd; ip++) {
            msign = h_modeSign(2,ip);
            nsh = ip-3; /* number of shapefunctions of each order */
            for(is =0; is< nsh; is ++) {
              EN_getDataInt((pEntity) Rfaces[0], MYCTID, &myctint);
              ien[blockid][iel[blockid]][ldof++] = myctint + j++;
              if(fsign != GenericfaceDir[0] && msign == 0){
                ien[blockid][iel[blockid]][ldof-1] *= -1;
              }
            }
          }
        }
      }
    }

    /* Region Modes */

    if(info->regnson){
      EN_getDataInt((pEntity) region, NDOFID, &nrm);
      if (nrm > 0)
        for (j=0; j < nrm; j++)  {
          EN_getDataInt((pEntity) region, MYCTID, &myctint);
          ien[blockid][iel[blockid]][ldof++] = myctint + j;
        }
    }

    iel[blockid]++; /* increment the element counter for this block on
                            this processor */

  }
  RIter_delete(rIter);

  /******************** boundary elements ********************/

  void* temp = 0;
  while (face =(pFace) PList_next(bdry,&temp)){
    /* find the list of regions connected to this face. Since this
       is a boundary face, there should only be one such region. */
    ents = F_regions(face);
    region = (pRegion)PList_item(ents,0); /* get the region */
    PList_delete(ents);

    EN_getDataInt((pEntity)region, RNENID, &nen);  
    EN_getDataInt((pEntity)region, POLYID, &poly); 
    BLOCK.nen = nen;
    BLOCK.maxpoly = poly;
    BLOCK.nenbl = F_numEdges(face);
    BLOCK.lcsyst = topology(region);
    if (BLOCK.lcsyst == 3) BLOCK.lcsyst = BLOCK.nenbl;
        if (BLOCK.lcsyst == 5)
                if (BLOCK.nenbl == 3) BLOCK.lcsyst = 6;

    blockid = Bblock[BLOCK]-1;
    nedges = nen+(nen+1)/2;
    nfaces = nen == 8 ? 6: (nen-(nen/6));

    /* get the mesh bounding mesh entities */
    R_entitiesBdry(region,face,Rvertices,info->edgeson?Redges:0,info->faceson?Rfaces:0,nen);

    /*********** gather equation numbers ***********/
    /* vertices */

    for (i=0; i < nen; i++){
      EN_getDataInt((pEntity) Rvertices[i], MYCTID, &myctint);
      ienb[blockid][ielb[blockid]][i] = myctint;
    }

    ldof = nen;

    /* All the stuff below is same as for interior elements , just
       look the same comments... someday should find elegant way of
       reusing code for this */


    if(info->edgeson) {
      switch(nen){
      case 4:   /* Tetrahedra */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = TetEMap[i][0];
          GenericEdgeMap[i][1] = TetEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = TetFMap[i][j];

        break;
      case 5:   /* Pyramid */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = PyrEMap[i][0];
          GenericEdgeMap[i][1] = PyrEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = PyrFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = PyrfaceDir[i];

        break;

      case 6:  /* Wedge Element */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = WedEMap[i][0];
          GenericEdgeMap[i][1] = WedEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = WedFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = WedfaceDir[i];

        break;

      case 8:  /* Hexahedron */

        for(i=0; i< nedges; i++){
          GenericEdgeMap[i][0] = HexEMap[i][0];
          GenericEdgeMap[i][1] = HexEMap[i][1];
        }

        for(i=0; i< nfaces; i++)
          for(j=0; j<4;  j++)
            GenericFaceMap[i][j] = HexFMap[i][j];

        for(i =0; i< nfaces; i++)
          GenericfaceDir[i] = HexfaceDir[i];

        break;
      }

      nem = poly-1;
      for (i=0; i < nedges; i++){
          v1 = Rvertices[GenericEdgeMap[i][0]];
          esign = (v1 == E_vertex(Redges[i],0));
          for (j=0; j < nem; j++){
            msign = h_modeSign(1,j+2);
            EN_getDataInt((pEntity) Redges[i], MYCTID, &myctint);
            ienb[blockid][ielb[blockid]][ldof++] = myctint + j;
            if((esign == 0) &&(msign ==0))
              ienb[blockid][ielb[blockid]][ldof-1] *= -1;
          }
      }
    }

    if(info->faceson){
      /* look at the detailed comment for interior elements */

      for(i =0; i < nfaces; i++){
        j =0;
        nvf = F_numEdges(Rfaces[i]);
        PolyOrd = poly;
        if (nvf == 3){
          for (k=0; k < 3; k++)
          v[k] = Rvertices[GenericFaceMap[i][k]];
          fsign = Fsign(Rfaces[i],v);
          for (ip = 3; ip <= PolyOrd; ip++){
            msign = modeSign(2,ip);
            nsh = ip-2;
            for (is=0; is < nsh; is++){
              EN_getDataInt((pEntity) Rfaces[i], MYCTID, &myctint);
              ienb[blockid][ielb[blockid]][ldof++] = myctint + j++;
              if ((fsign == 0) && (msign == 0))
                ienb[blockid][ielb[blockid]][ldof-1] *= -1;
            }
          }
        } else if(nvf == 4){
          fsign = R_dirUsingFace(region,Rfaces[0]);
          for(ip = 4; ip <= PolyOrd; ip++) {
            msign = h_modeSign(2,ip);
            nsh = ip-3; /* number of shapefunctions of each order */
            for(is =0; is< nsh; is ++) {
              EN_getDataInt((pEntity) Rfaces[0], MYCTID, &myctint);
              ienb[blockid][ielb[blockid]][ldof++] = myctint + j++;
              if(fsign != GenericfaceDir[0] && msign == 0){
                ienb[blockid][ielb[blockid]][ldof-1] *= -1;
              }
            }
          }
        }
      }
    }

    /* Region Modes */

    if(info->regnson){
      EN_getDataInt((pEntity) region, NDOFID, &nrm);
      if (nrm > 0)
        for (j=0; j < nrm; j++) {
          EN_getDataInt((pEntity) region, MYCTID, &myctint);
          ienb[blockid][ielb[blockid]][ldof++] = myctint + j;
        }
    }

    ielb[blockid]++;
  }
  ielb.clear();

  /* free the counter array */
  for(i=0; i< 6; i++) delete [] GenericFaceMap[i];
  for(i=0; i< 12; i++) delete [] GenericEdgeMap[i];
  delete [] GenericfaceDir;
  delete [] GenericEdgeMap;
  delete [] GenericFaceMap;
  delete [] Rvertices;
  delete [] Redges;
  delete [] Rfaces;

}

/* ---------------------------------------------------------------
   return if the shape function associated with the given mode of
   a finite element entity is ODD/EVEN
   Makes sense only for edge/face modes, for regions/vertex its
   always even.

   return 1 for EVEN mode
          0 for ODD  mode

   type = 1 ==> edge mode
          2 ==> tri-face mode

   pdof is the spectral order for the entity.
----------------------------------------------------------------- */
int modeSign(int type, int pdof)
{
  if(type == 1)         /* edge */
    return (pdof%2 ? 0 : 1);
  else if(type == 2)    /* tri face */
    return ((pdof%3)%2 ? 0 : 1) ;
  else
    return 1;
}

int h_modeSign(int type, int pdof)
{
  /* read comments above , pretty much the same thing but for hexes */
  if(type == 1) /* edge */
    return ((pdof%2) ? 0 : 1);
  else if(type == 2)  /* quad face */
    return ((pdof%2) ? 0 : 1) ;
  else
    return 1;
}

/*
return the direction that a face is being used with respect to
how it was defined
*/
int Fsign(pFace face,pVertex v[3])
{
  pPList verts;
  pVertex vl[3];
  int j;

  /* get the vertices in the order they were defined */
  verts = F_vertices(face,1);
  for (j=0; j < 3; j++)
    vl[j] = (pVertex)PList_item(verts,j);
  PList_delete(verts);

  /* compare the vertices */
  if ((v[0]==vl[0] && v[1]==vl[1])
      || (v[0]==vl[1] && v[1]==vl[2])
      || (v[0]==vl[2] && v[1]==vl[0])) return 1;
  else return 0;
}

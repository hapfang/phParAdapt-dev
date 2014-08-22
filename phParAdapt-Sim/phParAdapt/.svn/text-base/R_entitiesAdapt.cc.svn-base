/**********************************************************************/
/* This function gathers all the mesh entities in the closure of a    */
/* tetrahedron.                                                       */
/**********************************************************************/
#include "phParAdapt.h"
#include "MeshSimInternal.h"
#include <stdlib.h>
#include <stdio.h>


void 
R_entitiesAdapt( pRegion region, 
            pVertex *vrts, 
            pEdge *edgs,
            pFace *fcs ) {

    int faceVmap[3][2]={ {0,1}, {1,2}, {2,0} };  
    pPList ents,verts,vts,flist,list2;
    pFace face,face2;
    int i,j,k,idir;
    pPList wverts;
    
    /* vertices on this region */
    ents = R_vertices( region, 0);
    for(i=0; i < 4; i++)
        vrts[i] = (pVertex)PList_item(ents,i);
    PList_delete(ents);
    
    /* find the base face */
    if (!(face = F_exists(Tvertex,( pEntity )vrts[0],
                          ( pEntity )vrts[1],( pEntity )vrts[2],0)) ){
        fprintf(stderr,"Error: face not found in R_entities...");
        exit(-1);
    }
    
    /* six edges on this region */
    for (i=0; i < 3; i++)        /* 3 edges on base face */
        edgs[i] = E_exists( vrts[faceVmap[i][0]],
                            vrts[faceVmap[i][1]] );
    
    edgs[3] = R_gtOppEdg( region, edgs[1]); /* other three edges */
    edgs[4] = R_gtOppEdg( region, edgs[2]);
    edgs[5] = R_gtOppEdg( region, edgs[0]);
  
    /* four faces on this region */
    fcs[0] = face;
    fcs[1] = R_vtOpFc( region, vrts[2]);
    fcs[2] = R_vtOpFc( region, vrts[0]);
    fcs[3] = R_vtOpFc( region, vrts[1]);
}

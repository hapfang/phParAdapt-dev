/**********************************************************************/
/* This function gathers all the mesh entities in the closure of a    */
/* mesh region                                                        */
/**********************************************************************/
#include "func.h"
#include "parallel.h"
#include "phParAdapt.h"
#include "MeshSimInternal.h"


void R_entities(pRegion region, pVertex *vrts, pEdge *edgs, pFace *fcs, int nen)
{
  /* Note: edgs, fcs will be passed as zero if !edgeson, !faceson -- SST */
  int faceVmap[3][2]={{0, 1}, {1, 2}, {2, 0}};
  pPList ents, list, flist, list2;
  pFace face, face2;
  int i, k;
  int direction;

  /* vertices on this region */
  switch (nen) {
  case 4:
    {
      ents = R_vertices(region, 0);
      for(i=0; i < 4; i++)
        vrts[i] = (pVertex)PList_item(ents, i);
      PList_delete(ents);

      if (!(face = F_exists(Tvertex, vrts[0], vrts[1], vrts[2], 0))){
        fprintf(stderr, "Error: face not found in getConnectivity...");
        exit(-1);
      }

      if (edgs)  {
        /* six edges on this region */
        for (i=0; i < 3; i++)        /* 3 edges on base face */
          edgs[i] = E_exists(vrts[faceVmap[i][0]], vrts[faceVmap[i][1]]);

        edgs[3] = R_gtOppEdg(region, edgs[1]); /* other three edges */
        edgs[4] = R_gtOppEdg(region, edgs[2]);
        edgs[5] = R_gtOppEdg(region, edgs[0]);
      }

      if (fcs)  {
        /* four faces on this region */
        fcs[0] = face;
        fcs[1] = R_vtOpFc(region, vrts[2]);
        fcs[2] = R_vtOpFc(region, vrts[0]);
        fcs[3] = R_vtOpFc(region, vrts[1]);
      }
    }
    break;
  case 5:
    {
      face = R_face(region, 0);
      /* get the vertices on the first face of the region
         so that the normal of the face points out*/
      list = F_vertices(face, 1-R_dirUsingFace(region, face));

      /* to get the fifth vertex */
      ents = R_vertices(region,1);
     


      for(i=0; i < 5; i++) {
	if (!PList_inList (list, PList_item(ents, i))) {
          PList_append (list, PList_item (ents, i));
          break;
        }
      }
      PList_delete(ents);

      for (i=0; i < 5; i++)
        vrts[i] = (pVertex)PList_item(list, i);
      if((!(face = F_exists(Tvertex, vrts[0], vrts[1], vrts[2], vrts[3])))) {
        fprintf(stderr, "Error: face not found in getConnectivity...");
        exit(-1);
      }

      if (edgs)  {
        edgs[0] = E_exists(vrts[0], vrts[1]);
        edgs[1] = E_exists(vrts[1], vrts[2]);
        edgs[2] = E_exists(vrts[2], vrts[3]);
        edgs[3] = E_exists(vrts[3], vrts[0]);
        edgs[4] = E_exists(vrts[0], vrts[4]);
        edgs[5] = E_exists(vrts[1], vrts[4]);
        edgs[6] = E_exists(vrts[2], vrts[4]);
        edgs[7] = E_exists(vrts[3], vrts[4]);
      }

      if (fcs)  {
        fcs[0]=face;
        fcs[1]=F_exists(Tvertex, vrts[0], vrts[1], vrts[4], 0);
        fcs[2]=F_exists(Tvertex, vrts[1], vrts[2], vrts[4], 0);
        fcs[3]=F_exists(Tvertex, vrts[2], vrts[3], vrts[4], 0);
        fcs[4]=F_exists(Tvertex, vrts[3], vrts[0], vrts[4], 0);
      }
    }
    break;

  case 6:
    {

      ents = R_vertices(region,1);




      for(i=0; i < 6; i++)
        vrts[i] = (pVertex)PList_item(ents, i);

      PList_delete(ents);
      if((!(face = F_exists(Tvertex, vrts[0], vrts[2], vrts[1], 0)))){
        fprintf(stderr, "Error: face not found in getConnectivity...");
        exit(-1);
      }

      if (edgs)  {
        edgs[0] = E_exists(vrts[0], vrts[1]);
        edgs[1] = E_exists(vrts[1], vrts[2]);
        edgs[2] = E_exists(vrts[2], vrts[0]);
        edgs[3] = E_exists(vrts[3], vrts[4]);
        edgs[4] = E_exists(vrts[4], vrts[5]);
        edgs[5] = E_exists(vrts[5], vrts[3]);
        edgs[6] = E_exists(vrts[0], vrts[3]);
        edgs[7] = E_exists(vrts[1], vrts[4]);
        edgs[8] = E_exists(vrts[2], vrts[5]);
      }

      if (fcs)  {
        fcs[0]=face;
        fcs[1]=F_exists(Tvertex, vrts[0], vrts[1], vrts[4], vrts[3]);
        fcs[2]=F_exists(Tvertex, vrts[1], vrts[2], vrts[5], vrts[4]);
        fcs[3]=F_exists(Tvertex, vrts[2], vrts[0], vrts[3], vrts[5]);
        fcs[4]=F_exists(Tvertex, vrts[3], vrts[4], vrts[5], 0);
      }
    }
    break;
  
  case 8:
    {


        flist = R_faces(region,1);





      face = R_face(region, 0);

      /* get the vertices on the first face of the region */
      direction=R_dirUsingFace(region, face);
      list = F_vertices(face, 1-direction);

      for (i=0; i < 4; i++) {
        vrts[i] = (pVertex)PList_item(list, i);
      }

      /* find the face on the opposite side of the hex */

      face2 = (pFace)PList_item(flist, 5);

      /* find the vertex on the other face matching vertex 0 */

      direction = R_dirUsingFace(region, face2);

      list2 = F_vertices(face2, direction);

      for (k=0; k < 4; k++) {
        for (i=0; i < 4; i++) {
          if (E_exists((pVertex)PList_item(list2, i), vrts[k])) {
            vrts[4+k] = (pVertex)PList_item(list2, i);
            break;
          }
        }
      }

      PList_delete(list);
      PList_delete(list2);
      if (!(face = F_exists(Tvertex, vrts[0], vrts[1], vrts[2], vrts[3]))){
        fprintf(stderr, "Error: face not found in getConnectivity...");
        exit(-1);
      }

      if (edgs)  {
        edgs[0] = E_exists(vrts[0], vrts[1]);
        edgs[1] = E_exists(vrts[1], vrts[2]);
        edgs[2] = E_exists(vrts[2], vrts[3]);
        edgs[3] = E_exists(vrts[3], vrts[0]);
        edgs[4] = E_exists(vrts[4], vrts[5]);
        edgs[5] = E_exists(vrts[5], vrts[6]);
        edgs[6] = E_exists(vrts[6], vrts[7]);
        edgs[7] = E_exists(vrts[7], vrts[4]);
        edgs[8] = E_exists(vrts[0], vrts[4]);
        edgs[9] = E_exists(vrts[1], vrts[5]);
        edgs[10] = E_exists(vrts[2], vrts[6]);
        edgs[11] = E_exists(vrts[3], vrts[7]);
      }

      if (fcs)  {
        fcs[0]=face;
        fcs[1]=F_exists(Tvertex, vrts[0], vrts[1], vrts[5], vrts[4]);
        fcs[2]=F_exists(Tvertex, vrts[1], vrts[2], vrts[6], vrts[5]);
        fcs[3]=F_exists(Tvertex, vrts[2], vrts[3], vrts[7], vrts[6]);
        fcs[4]=F_exists(Tvertex, vrts[3], vrts[0], vrts[4], vrts[7]);
        fcs[5]=F_exists(Tvertex, vrts[4], vrts[5], vrts[6], vrts[7]);
      }
    }
    break;
  }
}

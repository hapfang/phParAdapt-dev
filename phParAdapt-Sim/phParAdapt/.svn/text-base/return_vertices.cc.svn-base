#include "func.h"
#include "ccfunc.h"
#include <iostream>
using namespace std;

extern pMeshDataId RNENID;

pPList return_vertices(pRegion region)
{
  pPList vertexlist;
  int x, n;
  
  EN_getDataInt((pEntity)region, RNENID, &n);

  if (n==0)  {
    switch(topology(region)) {
      case 1: n = 4; break;
      case 5: n = 5; break;
      case 3: n = 6; break;
      case 2: n = 8; break;
    }
  }

  switch(n) {
  case 4: /* tetrahedra */


    vertexlist = R_vertices(region, 0);


    break;

  case 8: /* hexahedra */
    {
      vertexlist = PList_new();


      pPList flist = R_faces(region,1);

     


      pPList list1,list2;
      /* first face of the region , ideally the 0th face */
      pFace face0 =(pFace) PList_item(flist,0);

      // Right now I am just taking the 5th face as the opposite face
      // If any problems come up return to looking for the opp face
      // :)

      pFace face1 =(pFace) PList_item(flist,5); /* the opposite face */

      /* pFace tmpface ; */
      /* pVertex  tvertex; */

      int direction;
      pVertex tvtx[8];

      //      for(int y=0; y < 6 ; y++) {
      //         tmpface=(pFace) PList_item(flist,y);
      //         if((face0 != tmpface) && !F_conToFace(face0,tmpface)) {
      //           face1 =(pFace) PList_item(flist,y);
      //           break;
      //         }
      //       }

      PList_delete(flist);
      direction= R_dirUsingFace(region,face0);

      if(direction){
        list1 = F_vertices(face0,0);
      } else {
        list1 = F_vertices(face0,1);
      }

      for(x=0; x<4;x++) {
        tvtx[x]= (pVertex)PList_item(list1,x);
      }

      direction= R_dirUsingFace(region,face1);

      if(direction) {
        list2 = F_vertices(face1,1);
      } else {
        list2 = F_vertices(face1,0);
      }

      for(x=0;x<4;x++){
        for(int y=0;y<4;y++){
          if(E_exists((pVertex)PList_item(list2,y),tvtx[x])) {
            tvtx[4+x]=(pVertex)PList_item(list2,y);
            break;
          }
        }
      }

      for(x=0;x<8;x++) {
        PList_append(vertexlist,tvtx[x]);
      }

      PList_delete(list1);
      PList_delete(list2);

    } /* end of case loop for hex elements */
    break;

  case 5: /* pyramids */



      vertexlist = R_vertices(region,1);


    break;

  case 6: /* wedges */



      vertexlist = R_vertices(region,1);



    break;

  default:
    cerr << "Unknown element type " << n << " in return_vertices "<<endl;
    exit(-1);
  }
  return vertexlist;
}

#include "parallel.h"
#ifdef FMDB
#include "AOMD.h"
#endif
#ifdef SIM
#include "MeshTypes.h"
#endif

/* The following function retuns a list of edges on a mesh region */
/* not in any particular order */

#ifdef __cplusplus
extern "C" {
#endif

pPList nsR_edges(pRegion region)
{
  pPList eList = PList_new();
  pPList eList2;
  pFace  face;
  int i,nfaces;

  nfaces = R_numFaces(region);

  for(i=0; i< nfaces; i++){  /* iterating over faces */
    face = R_face(region,i);
    eList2 = F_edges(face,1,0);
    PList_appPListUnique(eList,eList2);
        PList_delete(eList2);
  }
  return eList;
}

#ifdef __cplusplus
}
#endif

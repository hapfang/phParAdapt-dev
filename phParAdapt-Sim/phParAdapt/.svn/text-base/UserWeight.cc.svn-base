#include "parallel.h"
#include "func.h"
// #ifndef SIM
// using namespace SCOREC_mesh;
// #endif

#ifdef __cplusplus
extern "C" {
#endif
extern pMeshDataId NDOFID;

void UserWeight(pMesh mesh)
{
/* commenting this out right now (that means, all regions will have weight of 1
 * by default). this is since if we do this here, initLocalInfo will have to be
 * before partitioning, that means inefficiency since one proc will likely do
 * most of the stuff. if this is really needed, move initlocalinfo before
 * partitioning, and make an AttachDataCommu object like for inodes in readBC.cc
 * so it's migrated. then refer to ndof thru the AttachDataId object. else
 * simply calculate a rough weight here based on region type, globalP etc. */
#if 0
  RIter rIter = M_regionIter(mesh);
  pRegion rgn;
  pFace face;
  pEdge edge;
  pVertex vertex;

  int ndof;
  pPList vertices, edges, faces;
  void* temp=0;

/*    User comments
    cout << " Doing weighted Partitioning based on the total number ";
         << " of DOFs per element \n";
         << " If you want it some other way , please change the function ";
         << " UserWeight "<<endl;  */

  while(rgn = RIter_next(rIter)){
    /* the region modes */
    EN_getDataInt((pEntity)rgn, NDOFID, &ndof);

    /* the face modes */
    faces = R_faces(rgn);
    temp =0;
    while(face = (pFace)PList_next(faces,&temp))
      int ndofint;
      EN_getDataInt((pEntity)face, NDOFID, &ndofint);
      ndof += ndofint;
    PList_delete(faces);

    /* the edge modes */
    edges = nsR_edges(rgn);
    temp =0;
    while(edge = (pEdge)PList_next(edges,&temp))
      int ndofint;
      EN_getDataInt((pEntity)edge, NDOFID, &ndofint);
      ndof += ndofint;
    PList_delete(edges);

    /* the vertex modes */
    vertices = R_vertices(rgn);
    temp=0;
    while(vertex = (pVertex)PList_next(vertices,&temp)) ndof++;
    PList_delete(vertices);

    R_setWeight(rgn, ndof);
  }
  RIter_delete(rIter);
#endif
}

#ifdef __cplusplus
}
#endif

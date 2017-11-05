//
//
// This function reads the internal boundary condition information and
// assigns it to the appropriate mesh vertices
//
// Internal boundary conditions are necessary when the user wants to
// set boundary condtions on mesh entities which aren't classified on
// a geometric boundary. For example, a trip strip at the inflow of
// a turbulent cavity is produced by setting boundary conditions on
// a row of nodes near the base of the cavity. The user must specify
// these nodes in a separate file.
//
#include <map>
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"
#include "phParAdapt.h"

using namespace std;
extern map<pair<string,pGEntity>,void *> GEntityDataContainerPtr;
extern string bc; 



#ifdef PARALLEL
// Currently works only if a single (unpartitioned) mesh was initially loaded
void readBC(pGModel model, pParMesh pmesh)
#else
void readBC(pGModel model, pMesh mesh)
#endif
{
  int      nn,ibc,num,i;
  pVertex  vertex;
  pGEntity  g_ent;
  double   *list; // for BC values
  void     *temp;
  FILE     *fb;

  // create an array of all vertices in the mesh
  temp=0;
  i=0;

  /******************** internal BC's ********************/
  /* To specify boundary conditions on internal nodes,   */
  /* such as a trip-strip, a file must be provided that  */
  /* contains:                                           */
  /*   nn                                                */
  /*   11 BC values                                      */
  /*   iBC-code node# (for each node)                    */
  /* which should be named: inodes.dat                   */
  /*******************************************************/
  fb = fopen("inodes.dat","r");
  list = new double [ 11 ];
  fscanf(fb,"%d",&nn);
  for (i=0; i < 11; i++)
    fscanf(fb,"%lf",&list[i]);

  GRIter griter = GM_regionIter(model);
  g_ent = (pGEntity) GRIter_next(griter);
  GRIter_delete(griter);

  // attach the BC values to the first geometric region
  GEN_attachDataP(g_ent,"bc  ",(void *)list);

  pMeshDataId inod = MD_newMeshDataId("ibc");

  /* This data will be transmitted upon partitioning. Maybe later we can write
   * inodes.dat.<n> type files along with the mesh partition files, and read
   * those while reading back an already partitioned mesh. */

#ifdef PARALLEL
  pAttachDataCommu adc = AttachDataCommu_new(1,0,1);
  if(adc){
      MD_setMeshCallback(inod, CBmigrateOut, pm_sendAnInt, adc);
      MD_setMeshCallback(inod, CBmigrateIn, pm_recvAnInt, adc);
      PM_setMigrId(pmesh, inod);
  }
  pMesh mesh = PM_mesh(pmesh, 0);
#endif

  if (PMU_rank() == 0)  {
    pVertex  *vArray = new pVertex [M_numVertices(mesh)];
    VIter vIter = M_vertexIter(mesh);
    while( vertex = VIter_next(vIter) ) vArray[i++]=vertex;
    VIter_delete(vIter);
    // read the internal nodes, and attach the new ibc code
    for (i=0; i < nn; i++){
      fscanf(fb,"%d %d",&ibc,&num);
      vertex = vArray[num-1];
      EN_attachDataInt((pEntity)vertex, inod, ibc);
    }
    delete vArray;
  }
  fclose(fb);
}

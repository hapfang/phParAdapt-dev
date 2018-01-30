#include "SimPartitionedMesh.h"
#include "SimParasolidKrnl.h"
#include "SimAdvMeshing.h"
#include "MeshSimAdapt.h"
#include "SimModel.h"
#include "SimUtil.h"
#include <stdio.h>
#include <set>

int main(int argc, char **argv)
{
  SimPartitionedMesh_start(&argc, &argv);
  SimParasolid_start(1);
  SimAdvMeshing_start();
  Sim_readLicenseFile(NULL);     // Fix this

  pProgress prog = Progress_new();
  Progress_setDefaultCallback(prog);
  pNativeModel nmodel = ParasolidNM_createFromFile("geom_nat.x_t", 0);
  pGModel model = GM_load("geom.smd", nmodel, prog);
  pParMesh pmesh = PM_load("geom.sms", model, prog);
  pMesh mesh = PM_mesh(pmesh, 0);

  std::set<int> edgeids;
  char fname[80];
  sprintf(fname, "edgeids.%d", PMU_rank());
  FILE *fin = fopen(fname, "r");
  int id;
  while (fscanf(fin, "%d", &id) == 1)
    edgeids.insert(id);

// two lines above are from Saurabh when he thought we were still putting one int per line and this probably reads until eof
// but Cameron convinced us to go with c++ vectors in the interim and that is one line below
/*  fseek (fin , 0 , SEEK_END);
  long lSize = ftell (fin);
  rewind (fin);

  fread(&edgeids[0],sizeof(int),edgeids.size(),fin);
*/

  fclose(fin);

  pMSAdapt msa = MSA_new(pmesh, 0);
  EIter eiter = M_edgeIter(mesh);
  while (pEdge e = EIter_next(eiter)) {
    if (edgeids.count(EN_id(e)) == 1)
      MSA_setRefineLevel(msa, e, 1);
  }
  MSA_adapt(msa, prog);
  MSA_delete(msa);

  PM_write(pmesh, "adapted.sms", prog);

  M_release(pmesh);
  GM_release(model);
  NM_release(nmodel);
  Progress_delete(prog);

  Sim_unregisterAllKeys();
  SimAdvMeshing_stop();
  SimParasolid_stop(1);
  SimPartitionedMesh_stop();
  return 0;
}

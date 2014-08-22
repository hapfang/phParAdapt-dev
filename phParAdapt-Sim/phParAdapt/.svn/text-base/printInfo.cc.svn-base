#include "func.h"
#include "parallel.h"

#ifdef __cplusplus
extern "C" {
#endif
extern int numParts;
extern int numTotParts;

#ifdef PARALLEL
void printInfo(pParMesh pmesh, globalInfo *info, int ipart)
#else
void printInfo(pMesh mesh, globalInfo *info)
#endif
{
#ifdef DEBUG_PER
  pEntity ent;
  int count=0;
  int nem, nfm, nrm, nsh, nsb;
  pVertex vertex;
  pEdge   edge;
  pFace   face;
  pRegion region;
  void *temp = 0;
  char fn[20];
  FILE *fp;

#ifdef PARALLEL
  pMesh mesh = PM_mesh(pmesh, 0);
  sprintf(fn, "info.out.%d", PMU_rank()*numParts+ipart+1);
#else
  strcpy(fn, "info.out");
#endif

  fp = fopen(fn, "w");

  fprintf(fp, "\n==========> ENSA-5C preprocessing information <==========\n\n"
             "\t%d (times n_flow) total degrees of freedom (dofs) \n", info->nshgTot);

  fprintf(fp, "Global mesh entity information \n"
             "------------------------------ \n"
             "\t%6d Vertices\n\t%6d Edges\n\t%6d Faces\n\t%6d Regions\n\n",
   M_numVertices(mesh), M_numEdges(mesh), M_numFaces(mesh), M_numRegions(mesh));

  fprintf(fp, "\n\nIndividual partition information \n"
                 "--------------------------------\n"
             "Total number of processors: %d\n\n", numTotParts);
  fprintf(fp, "  Processor %d \n", PMU_rank()*numParts+ipart);

  fprintf(fp, "  ------------ \n"
             "\t%6d elements\n\t%6d vertices\n\t%6d boundary elements\n\t%6d "
             "prescribed essential BCs\n\t%6d is the dimension of the local "
             "work array\n\n",
             info->numel, info->numnp, info->numelb, info->numpbc,
             info->nlwork);

  fprintf(fp, "\nLocal mesh information \n");
  fprintf(fp, "----------------------- \n\n");

  fprintf(fp, "Vertices: (id, master)\n");
  vIter = M_vertexIter(mesh);
  while (vertex = VIter_next(vIter)) 
    fprintf(fp, "V%-2d ---------> %d\n", count++, EN_ownerProc((pEntity)vertex));
  VIter_delete(vIter);

  temp = 0;
  count = 0;
  fprintf(fp, "\nEdges: (id, master)\n");
  eIter = M_edgeIter(mesh);
  while (edge = EIter_next(eIter)) 
    fprintf(fp, "V%-2d-V%-2d ---------> %-6d\n", EN_id(E_vertex(edge, 0)),
            EN_id(E_vertex(edge, 1)), EN_ownerProc((Pentity)ent))
  EIter_delete(eIter);

  temp = 0;
  count = 0;
  fprintf(fp, "\nFaces: (id, master)\n");
  fIter = M_faceIter(mesh);
  while (face = FIter_next(fIter))
    fprintf(fp, "%6d \t%6d\n", count++, EN_ownerProc((pEntity)face));
  FIter_delete(fIter);

  temp = 0;
  count = 0;
  fprintf(fp, "\nRegions: (id, master)\n");
  rIter = M_regionIter(mesh);
  while (region =RIter_next(rIter))  
    fprintf(fp, "%6d \t%6d\n", count++, EN_ownerProc((pEntity)region));
  RIter_delete(rIter);

  fclose(fp);
#endif
}

#ifdef __cplusplus
}
#endif

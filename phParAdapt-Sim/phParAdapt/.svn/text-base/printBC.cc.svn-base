#include <stdio.h>
#include "func.h"



void printPeriodicBC(pMesh mesh)
{
  pVertex vertex;
  FILE *fout;
  void *temp = 0;
  pMatch pm;
#ifdef PARALLEL
  char s[20];
  sprintf(s, "perBC.out.%d", PMU_rank()+1);
  fout = fopen(s, "w");
#else
  fout = fopen("perBC.out", "w");
#endif

  fprintf(fout, "Vertices: \n");
  fprintf(fout, "--------- \n");

  VIter vIter = M_vertexIter(mesh);
  while( vertex = VIter_next(vIter) ) {
    if (isPeriodic((pEntity)vertex, &pm))  {
      pPoint pt = V_point(vertex);
      fprintf(fout, "Slave:  (%8.5e %8.5e %8.5e)  ID = %4d\n",
                     P_x(pt), P_y(pt), P_z(pt), EN_id((pEntity)vertex));
      if (isOnThisProc(pm))  {
//old        pVertex vertm = (pVertex) PMatch_ent(pm);
          pVertex vertm = (pVertex) Match_ent(pm);  //new
        pPoint ptm = V_point(vertm);
        fprintf(fout, "Master: (%8.5e %8.5e %8.5e)  ID = %4d\n\n",
                      P_x(ptm), P_y(ptm), P_z(ptm), EN_id((pEntity)vertm));
      }
      else
        fprintf(fout, "Master: off-proc\n\n");
    }
  }
  VIter_delete(vIter);
  fclose(fout);
}

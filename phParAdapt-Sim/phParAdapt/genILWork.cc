#ifdef PARALLEL

/**********************************************************************/
/* this function generates the local work array for parallel */
/* communication for the given processor */
/**********************************************************************/
// called in writeEnsaFiles.cc: writeEnsaFiles
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "func.h"
#include "parallel.h"
#include "assert.h"

using namespace std;
#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId MYCTID;
extern pMeshDataId NDOFID;
extern int numParts;
extern int numTotParts;

void ilworkTask(Task *t, int *ilwork, int *count, FILE *fp)
{
  int i;
  int myctint, ndofint;
  pEntity ent;
  if (t)  {
    ilwork[0]++;
    ilwork[*count+1] = t->tag;
    ilwork[*count+2] = t->type;
    ilwork[*count+3] = t->other_pid + 1;
    ilwork[*count+4] = t->numSeg;
   
#ifdef DEBUG
    fprintf(fp, "%d %d %d %d\n", t->tag, t->type, t->other_pid, t->numSeg);
#endif
    for (i=0; i < t->numSeg; i++)  {
      ent = (pEntity) PList_item(t->ents, i);
      if(!EN_getDataInt(ent, MYCTID, &myctint)) cout << "Error ilwork" << endl;
      EN_getDataInt(ent, NDOFID, &ndofint);
      ilwork[*count + 5 + 2*i] = myctint + 1;  /* Fortran idx */
      ilwork[*count + 6 + 2*i] = ndofint;
    }
    *count += 4 + 2*t->numSeg;
  }
}


void genILWork(globalInfo *info, int *ilwork, int ipart)
{
  char fname[20];
  FILE *ilwf;
  int i, count;

  sprintf(fname, "graph.out.%d", PMU_rank()*numParts+ipart+1);
#ifdef DEBUG
  ilwf = fopen(fname, "w");
  fprintf(ilwf, "tag, type, other, numSeg\n");
#endif
  ilwork[0] = 0;
  count = 0;

  /* Exchange ordering of ilwork from all-recvs-follow-all-sends to vice-versa
   * by swapping rtask and stask below */
#ifdef DEBUG
  for (i=0; i < numTotParts; i++)
    ilworkTask(info->rtask[i], ilwork, &count, ilwf);
  for (i=0; i < numTotParts; i++)
    ilworkTask(info->stask[i], ilwork, &count, ilwf);
  assert(count+1 == info->nlwork);  /* +1 for ilwork[0] which has numtasks */
  fclose(ilwf);

#else
  for (i=0; i < numTotParts; i++)
    ilworkTask(info->rtask[i], ilwork, &count, NULL);
  for (i=0; i < numTotParts; i++)
    ilworkTask(info->stask[i], ilwork, &count, NULL);
  assert(count+1 == info->nlwork);  /* +1 for ilwork[0] which has numtasks */
#endif
}

#ifdef __cplusplus
}
#endif

#endif


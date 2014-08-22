/**********************************************************************/
/* functions relating to global multiprocessor information */
/**********************************************************************/
#include <stdio.h>
#include "func.h"
#include "parallel.h"
extern int numParts;
extern int numTotParts;
/**********************************************************************/
/* initialize the global interprocessor communication and information */
/* structure. */
//
// allocation and initialization to 0
// functions is called in mdb2phasta before partitioning
////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif

void initGlobalInfo(globalInfo *info)
{
  int j;

  info->numnp = info->numel = info->numpbc = info->numelb = info->nshgTot
              = info->nshgOwn = info->nlwork = info->nshg = 0;

  info->stask = new Task* [ numTotParts ];
  info->rtask = new Task* [ numTotParts ];
  for (j=0; j < numTotParts; j++) {
    info->stask[j] = 0;
    info->rtask[j] = 0;
  }
  info->perSlv = PList_new();
  info->perMst = PList_new();
}
// called in addSendSegment() below
// called in addRecvSegment() below
int gettag(int slvproc, int mstproc)
{
  /* Let us have unique tags to each master-slave task based on the
   * slave pids.. thus we will have the range [PMU_size()*i .. PMU_size()*(i+1)-1]
   * for tags from proc i to proc 0..PMU_size()-1 when proc i is slave. If proc
   * j is master wrt proc i, for the reverse commu it will have proc i's tag
   * for proc j from the above range, i.e. PMU_size()*i + j. */
  return 0;
  return PMU_size()*slvproc + mstproc;
}

// segment is an entity
// slave is sending partition
//
// function is called in setupGlobalTasks.c
void addSendSegment(globalInfo *info, pEntity ent, pEntity ment, int pid)
{
  Task *t = info->stask[pid];
  if (!t)  {
    t = info->stask[pid] = new Task;
    t->tag = gettag(PMU_rank(), pid);
    t->type = 0;
    t->other_pid = pid;
    t->numSeg = 0;
    t->ents = PList_new();
    t->ments = PList_new();
  }
  t->numSeg++;
  PList_append(t->ents, ent);
  PList_append(t->ments, ment);
}
// master is receiving partition
// function is called in setupGlobalTasks.c
void addRecvSegment(globalInfo *info, pEntity ent, int pid)
{
  Task *t = info->rtask[pid];
  if (!t)  {
    t = info->rtask[pid] = new Task;
    t->tag = gettag(pid, PMU_rank());  /* order revsd from addSendSegment */
    t->type = 1;
    t->other_pid = pid;
    t->numSeg = 0;
    t->ents = PList_new();
  }
  t->numSeg++;
  PList_append(t->ents, ent);
}

void nlworkCalc(globalInfo *info)
{
#ifdef PARALLEL
  int i;
  Task *t;

  info->nlwork = 1;   /* 1 for numtasks */

  for (i=0; i < numTotParts; i++)  {
    /* 4 for tag,type,pid,numseg + begin,length for each numseg */
    if (t = info->stask[i])
      info->nlwork += 4 + 2*t->numSeg;
    if (t = info->rtask[i])
      info->nlwork += 4 + 2*t->numSeg;
  }
#endif
}

void freeGlobalInfo(globalInfo *info)
{
  int i;
  Task *t;

  for (i=0; i < numTotParts; i++)  {
    if (t = info->stask[i])  {
      PList_delete(t->ents);
      PList_delete(t->ments);
      delete t;
    }
    if (t = info->rtask[i])  {
      PList_delete(t->ents);
      delete t;
    }
  }
  delete [] info->stask;
  delete [] info->rtask;
  delete [] info->ncvec;
  PList_delete(info->perSlv);
  PList_delete(info->perMst);
}

#ifdef __cplusplus
}
#endif

#ifdef FMDB
#ifdef PARALLEL
#ifndef _PHPARTITIONCALLBACKS_H
#define _PHPARTITIONCALLBACKS_H

#include "AOMD.h"
#include "pmGraphs.h"
#include "pmZoltanCallbacks.h"

using namespace AOMD;
class phPartitionCallbacks: public AOMD::pmZoltanCallbacks
{
private:
  pMeshDataId phastaSolutionID;
  int nVar;
public :
  phPartitionCallbacks(pMeshDataId phastaSolID, int nV)
  : pmZoltanCallbacks(pmZoltanCallbacks::PARMETIS),phastaSolutionID(phastaSolID), nVar(nV) {}

  /* virtual void partition(AOMD::AOMD_distributed_graph2 &dg,
                         int *partitionVector); */
  virtual bool useWeights() const {return false;}
  virtual float getWeight (pEntity) const {return 1.0;}
  virtual void * getUserData (pEntity, int dest_proc, int &size);
  virtual void recieveUserData (pEntity, int pid, int tag, void *buf);
  virtual void deleteUserData (void *buf);
  virtual void deleteEntityData (pEntity ent);

#ifdef ENTITY_GROUP
   virtual void * getEntGrpUserData (mEntityGroup *, int dest_proc, int &size) { size = 0; return 0; }
   virtual void deleteEntGrpUserData (void *);
   virtual void deleteEntGrpData (mEntityGroup *);
   virtual void recieveEntGrpUserData (mEntityGroup *, int pid, int tag, void *buf);
#endif

};

#endif
#endif

#endif

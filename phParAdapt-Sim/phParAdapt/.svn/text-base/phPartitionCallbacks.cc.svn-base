#ifdef FMDB 
#ifdef PARALLEL

#include "phPartitionCallbacks.h"
#include "AOMD_cint.h"

/* void phPartitionCallbacks::partition(AOMD::AOMD_distributed_graph2 &dg,
                                     int *partitionVector)
{

  for(int i=0; i<dg.theGraph->nn; i++)
    partitionVector[i] = 0;

} */

void* phPartitionCallbacks::getUserData (pEntity ent, int dest_proc, int &size)
{

  if (EN_type(ent)!=0 || nVar ==0)
  { 
    size=0;
    return (0);
  }

  double *nodalData;
  if(nVar) {   
    size = nVar*sizeof(double);
    char *mybuffer = (char*) malloc (size);
    double *dbuf = (double*)malloc(nVar*sizeof(double));

    if(!EN_getDataPtr(ent,phastaSolutionID,(void **)&nodalData)) {
        printf("\nphasta solution not attached in phPartitionCallBacks\n");
        exit(0);
    }
    for (int i=0;i<nVar;i++)
      dbuf[i]= nodalData[i];
    memcpy(&mybuffer[0], dbuf,nVar*sizeof(double));
    
    free(dbuf);
    return mybuffer;
  }

  size = 0;
  return (0);
}

void phPartitionCallbacks::recieveUserData (pEntity ent, int dest_proc, int tag, void *buf)
{
  if (EN_type(ent)!=0 || nVar ==0)
    return;  

  char *mybuffer = (char*) buf;
  
  if(nVar) {
    double *dbuf = (double*)malloc(nVar*sizeof(double));
    memcpy(dbuf,&mybuffer[0], nVar*sizeof(double));
    double *value = new double[nVar];
    for (int i=0;i<nVar;i++)
      value[i]=dbuf[i];
#ifdef DEBUG
    double *nodalData;
    if(EN_getDataPtr(ent,phastaSolutionID,(void **)&nodalData)) {
        printf("\nphasta solution already attached to vertex in phPartitionCallBacks\n");
        exit(0);
    }
#endif
    EN_attachDataPtr(ent,phastaSolutionID,(void *)value);
    free(dbuf);
  }
}

void phPartitionCallbacks::deleteUserData (void *buf)
{
  free(buf);
}

void phPartitionCallbacks::deleteEntityData (pEntity ent)
{
  if(EN_type(ent)!=0)
    return;

  if(nVar) {
    double * nodalData;
    if(!EN_getDataPtr(ent,phastaSolutionID,(void **)&nodalData)) {
        printf("\nphasta solution not attached\n");
        exit(0);
    }
    delete [] nodalData;
    EN_deleteData(ent,phastaSolutionID);
  }
}


#ifdef ENTITY_GROUP
//   void * phPartitionCallbacks::getEntGrpUserData (mEntityGroup *, int dest_proc, int &size){}
   void phPartitionCallbacks::deleteEntGrpUserData (void *){}
   void phPartitionCallbacks::deleteEntGrpData (mEntityGroup *){}
   void phPartitionCallbacks::recieveEntGrpUserData (mEntityGroup *, int pid, int tag, void *buf){}
#endif


#endif
#endif

#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern int BCInflowFaceTag;

extern pMeshDataId vertID;
extern pMeshDataId inflowOwnedVertID;


// communicates contributions for verticies on inflow face and
// on partition bdry.

// owner vertex on inflow face and on part. bdry.
// (and have been tagged with vertID)
// take the id of owner vertex in connectivity

// data are gathered via
// vertID
void
commuBCInflowFaceInfo(pParMesh pmesh, pMesh mesh)
{
    int mult = sizeof(P_int)/sizeof(int);
    int *ns = new int [ PMU_size() ];
    int *nr = new int [ PMU_size() ];
    int *nsNodes = new int [ PMU_size() ];
    int *nrNodes = new int [ PMU_size() ];
    void ***se = new void** [ PMU_size() ];
    void ***re = new void** [ PMU_size() ];

    pEntity ent;
    int** dsend = new int*[PMU_size()];
    int** drecv = new int*[PMU_size()];

    for(int i=0; i<PMU_size(); i++) {
        ns[i] = 0;
	nsNodes[i] = 0;
	nrNodes[i] = 0;
    }

    int inflowOwnedVertTag;
    pVertex vertex;
    VIter vit = M_vertexIter(mesh);
    while(vertex = VIter_next(vit)) {
      if(EN_getDataInt((pEntity)vertex,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
	// need this for accumulatedID
	// ids of owned vertex by higher rank procs would
	// be added by accumulatedID
	for(int i=0; i<PMU_size(); i++)
	  nsNodes[i]++;
      }
    }
    VIter_reset(vit);

    // communicate the NUMBER of verts owned by each proc
    PMU_commuInt(nsNodes, nrNodes);

    int accumulatedID = 0;
    // add for the procs with lower ranks (excluding this proc)
    for(int i=0; i<PMU_rank(); i++)
      accumulatedID += nrNodes[i];

    VIter_reset(vit);
    while(vertex = VIter_next(vit)) {
      if(EN_getDataInt((pEntity)vertex,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
        int ID;
        // get the local data
        if(!EN_getDataInt((pEntity)vertex,vertID,&ID)) {
            cout<<"\nerror in commuBCInflowFaceInfo: no vertex id. attached to vertex (during offset)\n";
            exit(0);
        }

	// offset the IDs of owned vertices to make it global
	ID = ID + accumulatedID;
	EN_modifyDataInt((pEntity)vertex,vertID,ID);
      }
    }
    VIter_delete(vit);

    pEntity en;
    pBdryProcIter myBdryIter = PM_bdryProcIter(pmesh,0);
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt((pEntity)en,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
	pEntCopies vCopies = EN_copies(en); 

	// how many of them on OTHER parts/=procs ?
	int vCopiesSize = EntCopies_size(vCopies);
        
	// fill the send array with ents/vertices
	for(int i=0; i<vCopiesSize; i++) {

	  // global id of partition/proc on which en's n'th copy lies.
	  int vCopyGid = EntCopies_gid(vCopies,i);

	  if (vCopyGid != PMU_rank()) {
	    // increment number of copies shared between this proc and proc
	    // with rank vCopyGid
	    ns[vCopyGid]++;
	  }
	}
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
      }
    }
    BdryProcIter_reset(myBdryIter);
    
    for(int i=0; i < PMU_size(); i++) {
      // ns[i] is   number of ent-copies shared by this proc with
      // proc i

      // each se[i] is of length ns[i] and
      // therefore each se[i] will be filled with (off-proc) copies of all ents
      // that are shared with proc i

      // along with that we communicate an int array of length 2*ns[i]: dsend, drecv
      // -- id. of owned vertex
      // -- owner gid
        se[i] = new void* [ns[i]];
        dsend[i] = new int[ns[i]];
        ns[i] = nr[i] = 0;
    }

    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt((pEntity)en,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
        int ID;
        // get the local data
        if(!EN_getDataInt(en,vertID,&ID)) {
            cout<<"\nerror in commuBCInflowFaceInfo: no vertex id. attached to vertex (during send)\n";
            exit(0);
        }

        // copies of the vertex on other parts/=procs
        // does not include the entity en that vCopies was created with
        pEntCopies vCopies = EN_copies(en); 

        // how many of them on OTHER parts/=procs ?
        int vCopiesSize = EntCopies_size(vCopies );
        
        // fill the send array with ents/vertices
        for(int i=0; i<vCopiesSize; i++) {

	  // global id of partition/proc on which en's n'th copy lies.
	  int vCopyGid = EntCopies_gid(vCopies,i);

	  if (vCopyGid != PMU_rank()) {
	    // fill send array to proc of rank vCopyGid with the
	    // (ns[vCopyGid])-th copy
	    // among all copies of verts shared
	    // between this proc and proc vCopyGid
	    se[vCopyGid][ ns[vCopyGid] ] = EntCopies_ent(vCopies,i);
		
	    dsend[vCopyGid][ ns[vCopyGid] ] = ID;

	    // increment number of copies shared between this proc and proc
	    // with rank vCopyGid
	    ns[vCopyGid]++;
	  }
        }
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
      }
    }
    BdryProcIter_reset(myBdryIter);
    BdryProcIter_delete(myBdryIter);

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for(int i=0; i<PMU_size(); i++) {
        re[i] = new void* [nr[i]];
        drecv[i] = new int[nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_INT, 1);

    // retrieve communicated data and globalize
    for(int i=0; i<PMU_size(); i++) {
      // nr[i] is number of ent-copies shared by this proc with
      // proc i
      for(int j=0; j<nr[i]; j++) {            
	// the vertex itself
	ent = (pEntity) re[i][j];

	int ID;
	if(EN_getDataInt(ent,vertID,&ID)) {
	  cout<<"\nerror in commuBCInflowFaceInfo globalize: already have vertex id. data attached to vertex (during receive)\n";
	  exit(0);
	}
	
	ID = drecv[i][j];
	EN_attachDataInt(ent,vertID,ID);
      }
    }

    for(int i=0; i<PMU_size(); i++)  {
        delete [] se[i];
        delete [] re[i];
        delete [] dsend[i];
        delete [] drecv[i];
    }

    delete [] dsend;
    delete [] drecv;

    delete [] se;
    delete [] re;
    delete [] ns;
    delete [] nr;
    delete [] nsNodes;
    delete [] nrNodes;
}

#ifdef __cplusplus
}
#endif

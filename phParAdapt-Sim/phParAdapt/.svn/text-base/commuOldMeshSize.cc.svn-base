#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId oldMeshSizeID;
extern pMeshDataId numSurroundNodesID;


// communicates missing contributions for verticies
// on partition bdry. 

// each vertex on part. bdry. has accumulated size 
// of entities owned by the proc. (like edge lengths) around it
// AND also has number of entities that contributed in the accumulated size

// data are gathered via
// oldMeshSizeID
// numSurroundNodesID
void
commuOldMeshSize(pParMesh pmesh, pMesh mesh)
{
    int mult = sizeof(P_int)/sizeof(int);
    int *ns = new int [ PMU_size() ];
    int *nr = new int [ PMU_size() ];
    void ***se = new void** [ PMU_size() ];
    void ***re = new void** [ PMU_size() ];

    pEntity ent;
    double** dsend = new double*[PMU_size()];
    double** drecv = new double*[PMU_size()];

    for(int i=0; i<PMU_size(); i++) {
        ns[i] = 0;
    }

    pBdryProcIter myBdryIter = PM_bdryProcIter(pmesh,0);
    pEntity en;
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {
         pEntCopies vCopies = EN_copies(en); 

        // how many of them on OTHER parts/=procs ?
        int vCopiesSize = EntCopies_size(vCopies);
        
        // fill the send array with ents/vertices
        for(int i=0; i<vCopiesSize; i++) {

            // global id of partition/proc on which en's n'th copy lies.
            int vCopyGid = EntCopies_gid(vCopies,i);

	    if (vCopyGid != PMU_rank()) {
            // if (EntCopies_gid(vCopies,i) != PMU_rank()) {

                // increment number of copies shared between this proc and proc
                // with rank vCopyGid
                ns[vCopyGid]++;
            }
        }
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
    }
    BdryProcIter_reset(myBdryIter);
    
    for (int i=0; i < PMU_size(); i++) {
        // ns[i] is   number of ent-copies shared by this proc with
        // proc i and in turn, proc i shares the same number of bdry verts with
        // this proc

        // each se[i] is of length ns[i] and
        // therefore each se[i] will be filled with (off-proc) copies of all ents
        // that are shared with proc i

        // along with that we communicate a double array of length 2*ns[i]: dsend, drecv
        // -- numSurroundNodes (packed as double data)
        // -- oldMeshSize
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[2*ns[i]];
        ns[i] = nr[i] = 0;
    }

    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
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
            // if (EntCopies_gid(vCopies,i) != PMU_rank()) {

                // fill send array to proc of rank vCopyGid with the
                // (ns[vCopyGid])-th copy
	        // among all copies of verts shared
                // between this proc and proc vCopyGid
                se[vCopyGid][ ns[vCopyGid] ] = EntCopies_ent(vCopies,i);
                
                double* oldMeshSize;
                // get the local data
                if(!EN_getDataPtr(en,oldMeshSizeID,(void**)&oldMeshSize)) {
                    cout<<"\nerror in commuOldMeshSize: no old mesh size attached to vertex (during send)\n";
                    exit(0);
                }

                int* numSurroundNodes;
                if(!EN_getDataPtr(en,numSurroundNodesID,(void**)&numSurroundNodes)) {
                    cout<<"\nerror in commuOldMeshSize: no numSurroundNodes data attached to vertex (during send)\n";
                    exit(0);
                }

                // numSurroundNodes
                dsend[vCopyGid][ 2*ns[vCopyGid] ] = 1.0*(*numSurroundNodes);

		// old mesh size
		dsend[vCopyGid][ 2*ns[vCopyGid]+1 ] = *oldMeshSize;
               
                // increment number of copies shared between this proc and proc
                // with rank vCopyGid
                ns[vCopyGid]++;
            }
        }
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
    }
    BdryProcIter_reset(myBdryIter);
    BdryProcIter_delete(myBdryIter);

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for (int i=0; i<PMU_size(); i++) {
        re[i] = new void* [ nr[i] ];
        drecv[i] = new double[2*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, 2); 

    // retrieve communicated data and globalize
    for (int i=0; i<PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            // numSurroundNodes: add the incoming numSurroundNodes up to what
            // is already there in ent
            int* currentNumSurrNodes;
            if(!EN_getDataPtr(ent,numSurroundNodesID,
			      (void**)&currentNumSurrNodes)) {
	      cout<<"\nerror in commuOldMeshSize globalize: no currentNumSurrNodes data attached to vertex (during receive)\n";
	      exit(0);
            }

	    // num. surrounding nodes was packed as double during communication
	    *currentNumSurrNodes += (int)drecv[i][2*j];

            // get old mesh size (has been commuicated previous to this call)
            // and average it
            double* oldMeshSize;
            if(!EN_getDataPtr(ent,oldMeshSizeID,
			      (void**)&oldMeshSize)) {
                cout<<"\nerror in commuOldMeshSize globalize: no oldMeshSize data attached to vertex (during receive)\n";
                exit(0);
            }

	    *oldMeshSize += drecv[i][2*j+1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i=0; i<PMU_size(); i++)  {
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


    // average the accumulated old mesh size
    // do this for all vertices of the mesh
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;
    while(vertex = VIter_next(vIter)) {

        double* oldMeshSize; 
	if(!EN_getDataPtr((pEntity)vertex,oldMeshSizeID,
			  (void**)&oldMeshSize)) {
	  cout<<"\nerror in commuOldMeshSize: no oldMeshSize data attached to vertex (during average)\n";
	  exit(0);
	}
            
	int*  numSurrVerts;
	if(!EN_getDataPtr((pEntity)vertex,numSurroundNodesID,
			  (void**)&numSurrVerts)) {
	  cout<<"\nerror in commuOldMeshSize: no numSurrVerts data attached to vertex (during average)\n";
	  exit(0);
	}

	if(*numSurrVerts)
	  *oldMeshSize = *oldMeshSize/(*numSurrVerts);
	else {
	  V_info(vertex);
	  cout<<"\n["<<PMU_rank()<<"] Info : error in commuOldMeshSize: no surrounding entities found to compute average\n";
	  exit(0);
	}
    }
    VIter_delete(vIter);
}

#ifdef __cplusplus
}
#endif

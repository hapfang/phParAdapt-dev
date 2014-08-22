#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"


using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId partBdryFix4SolutionID;
extern int numVars;


// communicates contributions for verticies
// on partition bdry. that was fixed by fix4NodesOnSurface()

// each vertex on part. bdry. that was fixed by fix4NodesOnSurface()
// (and have been tagged with partBdryFix4SolutionID)
// has fixed solution attached; take average over all copies 
// to get unique solution value on all partitions

// data are gathered via
// phasta_solution
void
commuFix4SolutionTransfer(pParMesh pmesh, pMesh mesh)
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
    int partBdryFixTag;
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt(en,partBdryFix4SolutionID,&partBdryFixTag)) {

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
    }
    BdryProcIter_reset(myBdryIter);
    
    for (int i=0; i < PMU_size(); i++) {
        // ns[i] is   number of ent-copies shared by this proc with
        // proc i and in turn, proc i shares the same number of bdry verts with
        // this proc

        // each se[i] is of length ns[i] and
        // therefore each se[i] will be filled with (off-proc) copies of all ents
        // that are shared with proc i

        // along with that we communicate a double array of length numVars*ns[i]: dsend, drecv
        // -- nodalSolution (numVars)
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[numVars*ns[i]];
        ns[i] = nr[i] = 0;
    }

    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt(en,partBdryFix4SolutionID,&partBdryFixTag)) {
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
                
                double* nodalSol;
                // get the local data
                if(!EN_getDataPtr(en,phasta_solution,(void**)&nodalSol)) {
                    cout<<"\nerror in commuFix4SolutionTransfer: no solution attached to vertex (during send)\n";
                    exit(0);
                }

		for(int iVar=0; iVar<numVars; iVar++)
		  dsend[vCopyGid][ numVars*ns[vCopyGid]+iVar ] = nodalSol[iVar];

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

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for (int i=0; i<PMU_size(); i++) {
        re[i] = new void* [ nr[i] ];
        drecv[i] = new double[numVars*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, numVars); 

    // retrieve communicated data and globalize
    for (int i=0; i<PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            double* nodalSol;
            if(!EN_getDataPtr(ent,phasta_solution,
			      (void**)&nodalSol)) {
	      cout<<"\nerror in commuFix4SolutionTransfer globalize: no solution data attached to vertex (during receive)\n";
	      exit(0);
            }

	    for(int iVar=0; iVar<numVars; iVar++)
	      nodalSol[iVar] += drecv[i][numVars*j+iVar];
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

    // average the fixed solution on part. bdry. over all procs
    BdryProcIter_reset(myBdryIter);
    while(en =  BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt(en,partBdryFix4SolutionID,&partBdryFixTag)) {
        double* nodalSol; 
	if(!EN_getDataPtr(en,phasta_solution,
			  (void**)&nodalSol)) {
	  cout<<"\nerror in commuFix4SolutionTransfer: no solution attached to vertex (during average)\n";
	  exit(0);
	}
	
	pEntCopies vCopies = EN_copies(en); 
        int vCopiesSize = EntCopies_size(vCopies);

	for(int iVar=0; iVar<numVars; iVar++)
	  nodalSol[iVar] /= vCopiesSize;
	
	EN_deleteData(en,partBdryFix4SolutionID);
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
      }
    }
    BdryProcIter_delete(myBdryIter);
}

#ifdef __cplusplus
}
#endif

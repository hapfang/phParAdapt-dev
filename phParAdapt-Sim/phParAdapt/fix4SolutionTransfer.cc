#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "phParAdapt.h"



#ifdef __cplusplus
extern "C" {
#endif

// fix4NodesOnSurface introduces new vertices for which 
// a solution is needed

extern pMeshDataId partBdryFix4SolutionID;
extern pMeshDataId phasta_solution;
extern int numVars;


void
fix4SolutionTransfer(pMesh mesh)
{

  VIter vit=M_vertexIter(mesh);
  pVertex vertex;
  double* nodalData;

  int nodesCount = 0;
  int nodesCountOnPartBdry = 0;
  pPList fixedNodesWithNoSol = PList_new();
  while(vertex=VIter_next(vit)) { 
    // fix the solution for vertex
    // for vertex on part. bdry. communicate to copies and take average
    if(!EN_getDataPtr((pEntity)vertex, phasta_solution,
		      (void**)&nodalData)) {
      nodesCount++;
      PList_append(fixedNodesWithNoSol,vertex);

      // do not know whether this would happen or not
      // (i.e., new nodes created by fix4NodesOnSurface at partition bdry.)
      // in any case does not hurt
      if(EN_isOnPartBdry((pEntity)vertex)) {
	nodesCountOnPartBdry++;
	EN_attachDataInt((pEntity)vertex,partBdryFix4SolutionID,1);
      }
    }
  }
  VIter_delete(vit);

  printf("\n[%d] Fixing %d nodes in fix4SolutionTransfer()\n",PMU_rank(),nodesCount);
  printf("[%d] Fixing %d part. bdry. nodes in fix4SolutionTransfer()\n",PMU_rank(),nodesCountOnPartBdry);

  int iPass, numPasses = 5, countPasses = 0;
  pPList fixedNodesWithSol = PList_new();
  for(iPass=0; iPass<numPasses; iPass++) {
    if(!PList_size(fixedNodesWithNoSol))
      break;
    
    countPasses++;
    int sizeOfNodesWithNoSol = PList_size(fixedNodesWithNoSol);
    for(int iNode=0; iNode<sizeOfNodesWithNoSol; iNode++) {
      vertex = (pVertex)PList_item(fixedNodesWithNoSol,iNode);

      // if the node has NO solution, i.e. it is a NEW node
      if(!EN_getDataPtr((pEntity)vertex, phasta_solution,
			(void**)&nodalData)) {

	// need to find a neighbor for this one
	int nedges = V_numEdges( vertex);
	for(int j=0; j<nedges; j++) {
                
	  pEdge edge = V_edge(vertex,j);
                
	  pVertex vex = E_otherVertex(edge,vertex);

	  if(EN_getDataPtr((pEntity)vex, phasta_solution,
			   (void**)&nodalData)) {

	    double* sol = new double[numVars];
                    
	    for(int k=0; k<numVars; k++) {
	      sol[k] = nodalData[k];
	    }

	    EN_attachDataPtr( (pEntity)vertex, phasta_solution, (void *)
			      sol );

	    PList_append(fixedNodesWithSol,vertex);
	    break;
	  }
	}
      }
      else {
	printf("\n[%d] Error in fix4SolutionTransfer()... \n",PMU_rank());
	printf("node with solution attached is in list fixedNodesWithNoSol\n");
	exit(0);
      }
    }

    int sizeOfNodesWithSol = PList_size(fixedNodesWithSol);
    for(int iNode=0; iNode<sizeOfNodesWithSol; iNode++) {
      vertex = (pVertex)PList_item(fixedNodesWithSol,iNode);

      PList_remItem(fixedNodesWithNoSol,vertex);
    }
    PList_clear(fixedNodesWithSol);
  }

  if(PList_size(fixedNodesWithNoSol)) {
    printf("\n[%d] Error in fix4SolutionTransfer()... \n",PMU_rank());
    printf("%d have all surrounding vertices with no solution after [%d] passes\n",PList_size(fixedNodesWithNoSol),numPasses);
    exit(0);
  }
  else {
    printf("\n[%d] fix4SolutionTransfer done in %d pass(es)... \n",PMU_rank(),countPasses);
  }

  PList_delete(fixedNodesWithSol);
  PList_delete(fixedNodesWithNoSol);

}
                
#ifdef __cplusplus
}
#endif    



    
    

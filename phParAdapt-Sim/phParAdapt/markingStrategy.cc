#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"

using namespace std;

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;

// this routine tags/marks the mesh entities for refinement (i.e., tag driven)
// as of now only tags the edges (later, may introduce other choices)
// tags entities for refinement which have error values greater than threshold
// as of now do not use hmin and hmax
// can introduce one more factor to compute threshold for coarsening
int
applyMarkingStrategy(pMesh mesh, pMSAdapt simAdapter,
		     double factor, double hmax, double hmin,
		     double totalError, double maxError, double minError,
		     double threshold, int option) 
{

  if (option == 10 || option == 11) {
    threshold = factor;
  }
  else {   
  threshold = getErrorThreshold(mesh,factor,totalError,
				maxError,minError,option);
  }

  pEdge edge;    
  double *tagged = new double;
  pMeshDataId tagged_edges = MD_newMeshDataId("tagged edges");

  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
    double *nodalValues;
    if(!EN_getDataPtr((pEntity)vertex,errorIndicatorID,(void**)&nodalValues)) {
      cout<<"\nerror in applyMarkingStrategy(...) : no solution attached to vertex"<<endl;
      exit(0);
    }

    double errorValue;
    errorValue = getErrorValue(nodalValues,option);
//    errorValue = nodalValues[4];

    // If the error is above threshold
    // tag/mark all the edges of vertex for refinement
    // threshold is fraction of maximum error 
    if(errorValue > threshold) {
      int numEdges = V_numEdges(vertex);
      for (int i=0; i < numEdges; i++) {
	edge = V_edge(vertex,i);
	MSA_setRefineLevel(simAdapter,(pEntity)edge,1);
	if(!EN_getDataPtr((pEntity)edge,tagged_edges,(void**)&tagged)) {
	  EN_attachDataPtr((pEntity)edge,tagged_edges,(void *)tagged);
	}
      }
    }
  }
  VIter_delete(vIter);

  int edgesTagged = 0;
  EIter eIter = M_edgeIter(mesh);
  while(edge = EIter_next(eIter)) {
    if(EN_getDataPtr((pEntity)edge,tagged_edges,(void**)&tagged)) {
      edgesTagged++;
      EN_deleteData((pEntity)edge,tagged_edges);
    }
  }
  EIter_delete(eIter);

  delete tagged;
  MD_deleteMeshDataId(tagged_edges);

  return edgesTagged;
}
  
double
getErrorThreshold(pMesh mesh, double factor,
		  double totalError, 
		  double maxError, double minError,
		  int option)
{
  maxError = 0.;
  minError = 1.e15;
  totalError = 0.;
  double globalMax = 0;

  // std::ofstream errorLog("error.log");
  
  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
    double *nodalValues;
    if(!EN_getDataPtr((pEntity)vertex,errorIndicatorID,(void**)&nodalValues)) {
      cout<<"\nerror in getErrorThreshold(...) : no error data attached to vertex"<<endl;
      exit(0);
    }

    double errorValue = getErrorValue(nodalValues,option);

    if(errorValue > maxError) {
      maxError = errorValue;
    }

    if(errorValue < minError) {
      minError = errorValue;
    }

    totalError+=errorValue;
//    printf("totalError %lf\n",errorValue);

//     errorLog<<errorValue<<endl;
  }
  VIter_delete(vIter);

  MPI_Allreduce(&maxError, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  // errorLog.close();

  totalError = sqrt(totalError);

  return factor*globalMax;
} 

// option is to decide how to compute the error value
// (i.e., use 3 EI for flow problem or use 1 EI for scalar problem)
double 
getErrorValue(double *nodalValues, int option) {

  double errorValue = 0.;

  switch(option) {
  case 0:
    {
      double weight = 1.;
      errorValue = weight*nodalValues[0];  
    }
    break;
  case 1:
    {
      double weight[3]={1.,1.,1.};
      for (int i=1; i<4; i++) {
	errorValue += weight[i-1]*nodalValues[i]*nodalValues[i];
//       errorValue = nodalValues[4];
      }
      errorValue = sqrt(errorValue);
//    printf("error: %lf\n",errorValue);
    }
    break;
  default :
    cout<<"\nSpecify correct `option' to compute error value in getErrorValue(...)"<<endl;
    break;
  }
  return errorValue;
}

#ifdef __cplusplus
}
#endif
	

#include <iostream>
#include <stdlib.h>
#include "phParAdapt.h"
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId errorIndicatorID;

void SmoothErrorIndicators(pMesh mesh, int option)
{
  pVertex v;
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;

  double* averageNodalValue = new double[M_numVertices(mesh)];

  while(v = VIter_next(vIter)) {

    // get each vert's surrounding nodes 
    // if they are off-bdry
    int numEdges = V_numEdges(v);
    pEdge edge;
	
    averageNodalValue[vCount] = 0.0;

    // can be different from V_numEdges
    // will include the actual vertex, too
    int numSurroundingVerts=0;
        
    // exclude non-interor vertex from contributing to
    // its own average
    if(V_whatInType(v)== 3){
      // first the "central" vertex
      double *nodalEIs;
      if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&nodalEIs)){
	cout<<"\nerror in SmoothErrorIndicators: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      averageNodalValue[vCount] = averageNodalValue[vCount] + getErrorValue(nodalEIs,option);
      numSurroundingVerts++;
    }
        
    // now for the surrounding vertices
    for(int i=0; i<numEdges  ;i++){
      
      edge = V_edge(v,i);
      pVertex otherVertex;
      otherVertex = E_otherVertex(edge,v);
      
      // if neighbor vert is NOT on bdry
      if(V_whatInType(otherVertex)== 3){//is in interior
                
	double *nodalEIs;
	if(!EN_getDataPtr((pEntity)otherVertex,errorIndicatorID,(void**)&nodalEIs)){
	  cout<<"\nerror in SmoothErrorIndicators: no data attached to OTHER vertex\n";
	  V_info(otherVertex);
	  exit(0);
	}
	// add values up
	averageNodalValue[vCount] = averageNodalValue[vCount] + getErrorValue(nodalEIs,option);
	numSurroundingVerts++;
      }
    }//for( ..i<numEdges...)
    // if all surrounding verts on bdry that this verts' own value
    if(!numSurroundingVerts){

      double *nodalEIs;
      if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&nodalEIs)){
	cout<<"\nerror in SmoothErrorIndicators: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      averageNodalValue[vCount] =  getErrorValue(nodalEIs,option);
      numSurroundingVerts = 1;  


    }
    
    averageNodalValue[vCount] = averageNodalValue[vCount]/numSurroundingVerts;
    vCount++;

    }// while v = vIter
    VIter_reset(vIter); 

    vCount=0;
    while(v = VIter_next(vIter)) {
    
      // delete the old nodal EI values
      double* oldNodalEIs;
      if(!EN_getDataPtr((pEntity)v,errorIndicatorID,(void**)&oldNodalEIs)){
	cout<<"\nerror in SmoothErrorIndicators: no data attached to OTHER vertex\n";
	V_info(v);
	exit(0);
      }
//      delete [] oldNodalEIs;

      EN_deleteData((pEntity)v,errorIndicatorID);
    
//      double* averageEI = new double;
//      *averageEI = averageNodalValue[vCount];
      oldNodalEIs[4] = averageNodalValue[vCount];

      EN_attachDataPtr( (pEntity)v, errorIndicatorID, (void **)
			oldNodalEIs );

      vCount++;
    }
    VIter_delete(vIter); 

    delete [] averageNodalValue;

}

#ifdef __cplusplus
}
#endif

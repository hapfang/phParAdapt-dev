#include "phParAdapt.h"
//#include "MA.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalSizeID;
extern pMeshDataId numSurroundNodesID;
  
// simple average over a patch surrounding the vertex    
void 
SmoothSize(pMesh mesh) {

  pVertex v;
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  int* numSurroundingVerts = new int[M_numVertices(mesh)];
  // keep the vals in memory before finally setting them
  double** averageSize;
  averageSize = new double*[M_numVertices(mesh)];  
  
#if  ( defined  DEBUG )
//   printf("[%i] in SmoothHessians:\n",PMU_rank());
#endif 

  while(v = VIter_next(vIter)) {
    averageSize[vCount] = new double[6];
    
    for(int i=0;i<3;i++) {
      averageSize[vCount][i]=0.0;
    }
    // can be different from V_numEdges
    // will include the actual vertex, too
    numSurroundingVerts[vCount]=0;
    
    // exclude non-interor vertex from contributing to
    // its own average
    // also make sure this vertex is owned by this proc
    if(V_whatInType(v)== 3 && EN_isOwnerProc((pEntity)v) ) {
      // first the "central" vertex
      double *nodalSize;
      if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)){
	
	cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      for(int i=0; i<3; i++) {
	averageSize[vCount][i]= nodalSize[i];
      }
      numSurroundingVerts[vCount]++;

#if  ( defined  DEBUG )
//       double c[3];
//       V_coord(v,c);
//       printf("\n[%i]contributing to vert: coords %f %f %f\n",PMU_rank(),c[0],c[1],c[2]);
#endif     


    }//if(V_whatInType(v)== 3 && ..

    int numEdges = V_numEdges(v);
    pEdge edge;


    // now for the surrounding vertices
    // now this vertex v also can be non-owned
    // but EDGES to surrounding verts MUST be owned and verts also be off model-bdry
    for(int i=0; i<numEdges ; i++){
        edge = V_edge(v,i);

      if(EN_isOwnerProc((pEntity)edge)){

            pVertex otherVertex;
            otherVertex = E_otherVertex(edge,v);
            
        
#if  ( defined  DEBUG )
//             double cfq[3], cvq[3];
//             V_coord(otherVertex,cfq);
//             V_coord(v,cvq);
//             printf("\n[%i]considering other vertex (ORIG) (ownflag %i) %f %f %f:\n",PMU_rank(),EN_isOwnerProc((pEntity)v),cvq[0],cvq[1],cvq[2]);
//             printf("coords %f %f %f (orig) %f %f %f\n",cfq[0],cfq[1],cfq[2],cvq[0],cvq[1],cvq[2]);
#endif   

            // if neighbor vert is NOT on bdry NOR is it NOT owned (neighbor  HAS to be owned)
            if(V_whatInType(otherVertex)== 3){//is in interior
            
                double *nodalSizeOther;
                if(!EN_getDataPtr((pEntity)otherVertex, nodalSizeID,(void**)&nodalSizeOther)){
                    
                    cout<<"\nerror in SmoothSize: no data attached to OTHER vertex\n";
                    V_info(otherVertex);
                    exit(0);
                }
                // add values up
                for(int i=0;i<3;i++) {
                    averageSize[vCount][i] = averageSize[vCount][i] + nodalSizeOther[i];
                }
                numSurroundingVerts[vCount]++;
                
#if  ( defined  DEBUG )
//                 printf("\n[%i]is also vert (orig):\n",PMU_rank());
//                 double cf[3], cv[3];
//                 V_coord(otherVertex,cf);
//                 V_coord(v,cv);
//                 printf("coords %f %f %f (orig) %f %f %f\n",cf[0],cf[1],cf[2],cv[0],cv[1],cv[2]);
                
#endif  
            }//if(V_whatInType(otherVertex)== 3
        }// if(EN_isOwnerProc(edge){
            
    }//for( numEdges

    // At this moment we cannnot decide whether a vertex' neighbors are
    // exclusively on the bdry --> also look at adjacent partitions
    // for now, we rely that the initial mesh hase been prepared and got 
    // all its bdry-only elements removed

    // if no surronding vertex found 
    // this is not a good point to take care of this for partitioned mesh
    
    if(!numSurroundingVerts[vCount]) {
      double *nodalSize;
      if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)){
	cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      
      // just add its own values up
      // mostly the case when a region has no intr. node/vert.
      for(int i=0;i<3;i++) {
       averageSize[vCount][i] = nodalSize[i];
      }
      numSurroundingVerts[vCount]++;
    }

    //divide by number of Surr vertices
  for(int i=0;i<3;i++) {
     averageSize[vCount][i] = averageSize[vCount][i]/(numSurroundingVerts[vCount]);
  }


//      if(!numSurroundingVerts[vCount]){
//        cout<<"\nerror in SmoothHessians: there is a boundary vertex whose\n"
//  	  <<"neighbors are exclusively classsified on modelfaces/edges/vertices\n"
//  	  <<"and NOT in the interior\n";
      
//        cout<<"For the following vertex : "<<endl;
//        V_info(v);
//        exit(0);
//      }
    vCount++;
  }// while v = vIter
  VIter_reset(vIter); 
  
  vCount=0;
  while(v = VIter_next(vIter)) {
    // delete the old Hessian data
    double* oldSize;
    if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&oldSize)){
      
      cout<<"\nerror in SmoothSize: no data attached to OTHER vertex\n";
      V_info(v);
      exit(0);
    }
//      delete [] oldHessian;
    
//      EN_deleteData((pEntity)v, nodalHessianID);
   
    
    double* smoothSize = new double[3];

    for(int k=0;k<3;k++){
        smoothSize[k] = averageSize[vCount][k];
    }

    delete [] oldSize;
//    EN_deleteData( (pEntity)v, nodalSizeID);
    EN_modifyDataPtr( (pEntity)v, nodalSizeID, (void *)smoothSize);

    // 
    double* nSurrNodes = new double[1];
    nSurrNodes[0] = 1.0*numSurroundingVerts[vCount] ; 
    EN_attachDataPtr( (pEntity)v, numSurroundNodesID, (void *)
                       nSurrNodes);
    
#if  ( defined  DEBUG )
//             printf("\n[%i]attaching in  SmoothHessians:\n %f %f %f\n",PMU_rank(),smoothHess[0],smoothHess[1],smoothHess[2]);
//             printf(" %f %f %f\n",smoothHess[3],smoothHess[4],smoothHess[5]);
//             printf("\nand local numSurrVertex Num %f\n",nSurrNodes[0]);
//             printf("for vertex:");
//             double c[3];
//             V_coord(v,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif     

    vCount++;
  }
  VIter_delete(vIter);

  for(int i=0;i<M_numVertices(mesh);i++){
      delete [] averageSize[i];
  }
  delete [] averageSize;
  delete [] numSurroundingVerts;
//  M_writeVTKFile(mesh, "nodalSizeSmooth", nodalSizeID, 3);

  int numExtraSmoothPassesBdryFaces = 1;
  int numExtraSmoothPassesBdryEdges = 1;
  int numExtraSmoothPassesBdryVerts = 1;

  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryFaces,2);
  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryEdges,1);
  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryVerts,0);


}
  
// performs "numPasses" cycles of simple average over 
// a patch (consisting specific nodes) surrounding the 
// vertex classified on model entities of dimension "genDim"
// nodes contained in the surrounding patch are the ones 
// classified on model ents. with dim. higher than "genDim"
void SmoothSizeOnBdry(pMesh mesh, int numPasses, int genDim) {

  pMeshDataId avgSizeID = MD_newMeshDataId("avg. at bdry.");
  for(int iPass=0; iPass<numPasses; iPass++) {
    pVertex v;
    VIter vit = M_vertexIter(mesh);
    while(v = VIter_next(vit)) {
      // vertex on domain bdry. of dimension "genDim"
      if(V_whatInType(v)==genDim && EN_isOwnerProc((pEntity)v) ) {
	
   double *avgSize = new double[3];
	for(int i=0;i<3;i++)
	  avgSize[i] = 0.;

	// start from 1 due to its own contri.
	int numSurroundingVerts = 1;
	double *nodalSize;
	if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex\n";
	  V_info(v);
	  exit(0);
	}

	for(int i=0;i<3;i++)
	  avgSize[i] += nodalSize[i];

	int numEdges = V_numEdges(v);
	pEdge edge;
	for(int i=0; i<numEdges; i++) {
	  edge = V_edge(v,i);
     if(EN_isOwnerProc((pEntity)edge)){
   	  pVertex otherVertex;
	     otherVertex = E_otherVertex(edge,v);
      
	     // do not account for vertices
	     // classified on model ents. with 
	     // dimension lower than "genDim"
	    if(V_whatInType(otherVertex)<genDim)
       continue;

	    double *nodalSizeOther;
       if(!EN_getDataPtr((pEntity)otherVertex, nodalSizeID,(void**)&nodalSizeOther)){
	     cout<<"\nerror in SmoothSizeOnBdry: no data attached to OTHER vertex of bdry. vertex\n";
	     V_info(otherVertex);
	     exit(0);
	    } 

	    // add values up
	    for(int i=0;i<3;i++)
	      avgSize[i] += nodalSizeOther[i];
	    numSurroundingVerts++;
     } //if owned edge 
	}//for( ..i<numEdges...)

	for(int i=0;i<3;i++)
	  avgSize[i] /= numSurroundingVerts;

	EN_attachDataPtr((pEntity)v,avgSizeID,(void *)avgSize);
      }
    }
    VIter_reset(vit);
    
    VIter_reset(vit);
    while(v = VIter_next(vit)) {
      if(V_whatInType(v)==genDim && EN_isOwnerProc((pEntity)v)) {
	double *nodalSize, *avgSize;

	if(!EN_getDataPtr((pEntity)v, nodalSizeID,
            (void**)&nodalSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex\n";
	  V_info(v);
	  exit(0);
	}
	
	if(!EN_getDataPtr((pEntity)v, avgSizeID,(void**)&avgSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex with avghessianID\n";
	  V_info(v);
	  exit(0);
	}

//	EN_deleteData((pEntity)v, nodalSizeID);
	// attach the new data ptr. (i.e., averaged hessian)
	EN_modifyDataPtr((pEntity)v, nodalSizeID, (void *)avgSize);
	// delete the old data ptr. (i.e., old nodal hessian)
	delete [] nodalSize;
	
	// clean the ptr. attached with avghessianID for other passes 
	EN_deleteData((pEntity)v, avgSizeID);
      }
    }
    VIter_delete(vit);
  }

  MD_deleteMeshDataId(avgSizeID);
//  M_writeVTKFile(mesh, "nodalSizeSmooth2", nodalSizeID, 3);
}

void SizeLimit(pMesh mesh){

  pMeshDataId limitSizeID = MD_newMeshDataId("avg. at bdry.");
  pVertex v;

  char sizefile[256];
  int nshg;
  sprintf(sizefile,"OrgSize%i.dat", PMU_rank());
  FILE *sf;
  sf=fopen(sizefile, "rb");
  fscanf(sf,"%i",&nshg);
  printf("%i\n",nshg);
  
  printf("Limiting the size field request by factor of original size\n");
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  int count=0;
  double value[nshg*3]; 

  while(1){
     fscanf(sf,"%lf",&value[count]);
//     printf("%i %lf\n",count, value[count]);
     count++;
     if(feof(sf)) break;
  }
  
//  printf("%lf %lf %lf\n", value[0], value[1],value[2]);
  fclose(sf);
  while(v = VIter_next(vIter)) {
     
     double *h = new double[3];
     double anisoSize[3][3];
     if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&h)){

      cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
      V_info(v);
      exit(0);
     }
   
//     MA_Vtx_Size (v, anisoSize);

//put solution from value array in anisoSize 3x3 matrix for each vertex and then
//attach it to an     
     int vtxID;
     vtxID = EN_id(v);
     
    double* h2 = new double[3];
    for(int i=0; i<3; i++){
           h2[i]=value[vCount*3+i];
    }

    for (int i=0;i<3;i++){
       if(h[i]<h2[i]/4)
           h[i]=h2[i]/4;
//       if(h[i]>h2[i]/4) 
//           h[i]=h2[i]*4;
    }
//       if(h[0]>h2[0]/16) 
//          h[0]=h2[0]*16;
//       if(h[1]>h2[1]/16) 
//          h[1]=h2[1]*16;
       
//           printf("%f\n",h[i]);

   EN_attachDataPtr( (pEntity)v, limitSizeID, (void *)
                         h2);
   
//   EN_deleteData( (pEntity)v, nodalSizeID);
//   EN_modifyDataPtr( (pEntity)v, nodalSizeID, (void *)h);
   vCount++;
  }//vertex loop
  VIter_delete(vIter);
  M_writeVTKFile(mesh, "OrgSize", limitSizeID, 3);

}
#ifdef __cplusplus
}
#endif

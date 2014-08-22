#include "phParAdapt.h"
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

  
// simple average over a patch surrounding the vertex    
void 
SmoothWallField(pMesh mesh, pMeshDataId fieldID, int ndof) {

  pVertex v;
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  // keep the vals in memory before finally setting them
  double** averageField;
  averageField = new double*[M_numVertices(mesh)];  
  
#if  ( defined  DEBUG )
//   printf("[%i] in SmoothFields:\n",PMU_rank());
#endif 

  while(v = VIter_next(vIter)) {
    averageField[vCount] = new double[ndof];
    
    for(int iField=0;iField<ndof;iField++) {
      averageField[vCount][iField]=0.0;
    }
    // can be different from V_numEdges
    // will include the actual vertex, too
    int numSurroundingVerts=0;
    
    // exclude non-interor vertex from contributing to
    // its own average
    // also make sure this vertex is owned by this proc
    if(V_isOriginatingNode(v)){// && EN_isOwnerProc((pEntity)v) ) {
      // first the "central" vertex
      double *Field;
      if(!EN_getDataPtr((pEntity)v, fieldID,(void**)&Field)){
	
	cout<<"\nerror in SmoothFields: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      for(int iField=0; iField<ndof; iField++) {
	averageField[vCount][iField]= averageField[vCount][iField] + Field[iField];
      }
      numSurroundingVerts++;

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

//        if(EN_isOwnerProc((pEntity)edge)){

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
            if(V_isOriginatingNode(otherVertex)){//is an originating node
            
                double *FieldOther;
                if(!EN_getDataPtr((pEntity)otherVertex, fieldID,(void**)&FieldOther)){
                    
                    cout<<"\nerror in SmoothFields: no data attached to OTHER vertex\n";
                    V_info(otherVertex);
                    exit(0);
                }
                // add values up
                for(int iField=0;iField<ndof;iField++) {
                    averageField[vCount][iField] = averageField[vCount][iField] + FieldOther[iField];
                }
                numSurroundingVerts++;
                
#if  ( defined  DEBUG )
//                 printf("\n[%i]is also vert (orig):\n",PMU_rank());
//                 double cf[3], cv[3];
//                 V_coord(otherVertex,cf);
//                 V_coord(v,cv);
//                 printf("coords %f %f %f (orig) %f %f %f\n",cf[0],cf[1],cf[2],cv[0],cv[1],cv[2]);
                
#endif  
//            }//if(V_whatInType(otherVertex)== 3
        }// if(EN_isOwnerProc(edge){
            
    }//for( numEdges

    // At this moment we cannnot decide whether a vertex' neighbors are
    // exclusively on the bdry --> also look at adjacent partitions
    // for now, we rely that the initial mesh hase been prepared and got 
    // all its bdry-only elements removed

    // if no surronding vertex found 
    // this is not a good point to take care of this for partitioned mesh
    if(!numSurroundingVerts && V_isOriginatingNode(v)) {
      double *Field;
      if(!EN_getDataPtr((pEntity)v, fieldID,(void**)&Field)){
	cout<<"\nerror in SmoothFields: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      
      // just add its own values up
      // mostly the case when a region has no intr. node/vert.
      for(int iField=0;iField<ndof;iField++) {
       averageField[vCount][iField] = averageField[vCount][iField] + Field[iField];
      }
      numSurroundingVerts=1;
    }

//      if(!numSurroundingVerts[vCount]){
//        cout<<"\nerror in SmoothFields: there is a boundary vertex whose\n"
//  	  <<"neighbors are exclusively classsified on modelfaces/edges/vertices\n"
//  	  <<"and NOT in the interior\n";
      
//        cout<<"For the following vertex : "<<endl;
//        V_info(v);
//        exit(0);
//      }

    //divide by number of Surr vertices
    for(int i=0;i<ndof;i++) {
       if(numSurroundingVerts>0)
         averageField[vCount][i] = averageField[vCount][i]/numSurroundingVerts;      else
         averageField[vCount][i]=0.0;
     }


    vCount++;
  }// while v = vIter
  VIter_reset(vIter); 
  
  vCount=0;
  while(v = VIter_next(vIter)) {
    // delete the old Field data
    double* oldField;
    if(!EN_getDataPtr((pEntity)v, fieldID,(void**)&oldField)){
      
      cout<<"\nerror in SmoothFields: no data attached to OTHER vertex\n";
      V_info(v);
      exit(0);
    }
//      delete [] oldField;
    
    double* smoothField = new double[ndof];

    for(int k=0;k<ndof;k++){
        smoothField[k] = averageField[vCount][k];
    }
//    delete [] oldField;
//    EN_deleteData( (pEntity)v, fieldID);
    EN_modifyDataPtr( (pEntity)v, fieldID, (void *)
		       smoothField);

    // 
    
#if  ( defined  DEBUG )
//             printf("\n[%i]attaching in  SmoothFields:\n %f %f %f\n",PMU_rank(),smoothField[0],smoothField[1],smoothField[2]);
//             printf(" %f %f %f\n",smoothField[3],smoothField[4],smoothField[5]);
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
      delete [] averageField[i];
  }
  delete [] averageField;

}
  
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef SIM
#include "MeshSim.h"
#endif
#ifdef FMDB
#include "AOMD.h"
#endif
#include "MeshSimAdapt.h"
#include "phParAdapt.h"i
#include "func.h"


// implementation of the callback fctn
extern int CBcounter;
extern pMeshDataId phasta_solution;
extern pMeshDataId modes;
extern int numVars;
extern int DisplacementMigration;
int isEdgeBasedTransfer = 1;
//for measuring procsize

// arguments are :
// mesh mod. type, look at file MeshSimAdapt.h
// data containing mesh entities created, deleted or reshaped
// userData used to set a pointer to be passed in callback function
#ifdef SIM
void
phastaTransferTopSIM( MeshModType mtype, pMeshChanges mco, void *userData)  {

    CBcounter++;

    pPList fresh;
    pPList old;
    pRegion  rgn;
    pEntity ent;

    pVertex deleteVertex;
    int i,j,k,deletedVerts;
    
//    if(mtype!=3 && mtype!=0 && mtype!=4)    
//        printf("mtype: %d\n", mtype); 

//     sprintf(stringName, "adaptLog.%d",PMU_rank+1);
    //logfile=fopen(stringName,"a");
//    printf("\n entering callback fct\n");

    pPList freshEnts = MCO_createdEntities(mco); 
    pPList oldEnts  = MCO_deletedEntities(mco);

    int oldEntsSize = PList_size(oldEnts);
    int freshEntsSize = PList_size(freshEnts);

//     // this is done as callback function is invoked for
//     // some mesh mod. operations which are not successful
//     // this must be fixed with newer versions of Simmetrix libs.
//     if(!freshEntsSize) {
//       PList_delete(freshEnts);
//       PList_delete(oldEnts);

//       // do not count this call
//       CBcounter--;

//       return;
//     }
    if(mtype==3 || mtype==4)
       return;

    // old - holds old regions (to be deleted)
    old = PList_new();
    // freash - holds new regions (created)
    fresh = PList_new();
    // oldVertices - hold old vertices (to be deleted)
    pPList oldVertices = PList_new();
    pPList newVertices = PList_new();
    pPList oldEdges = PList_new();
    pPList newEdges = PList_new();
    pPList oldFaces = PList_new();
    pPList newFaces = PList_new();
    
    for(i=0; i<oldEntsSize; i++){
        ent = (pEntity)PList_item(oldEnts,i);
        int en_type = EN_type(ent);
        if(en_type==0) {
            PList_append(oldVertices,ent);
        }
        else if(en_type==1) {
            PList_append(oldEdges,ent);
        }
        else if(en_type==2) {
            PList_append(oldFaces,ent);
        }
        else if(en_type==3) {
            PList_append(old,ent);
        }
    }

    for(i=0; i<freshEntsSize; i++){
        ent = (pEntity)PList_item(freshEnts,i);
        int en_type = EN_type(ent);
        if(en_type==0) {
            PList_append(newVertices,ent);
        }
        else if(en_type==1) {
            PList_append(newEdges,ent);
        }
        else if(en_type==2) {
            PList_append(newFaces,ent);
        }
        if(en_type==3) {
            PList_append(fresh,ent);
        }
    }
    
    PList_delete(freshEnts);
    PList_delete(oldEnts);

    int ndof=numVars;
    
    if(mtype==0) {
      phastaTransferSIMBottom(newVertices, newEdges, newFaces, fresh, ndof, mtype);
    } else {
      phastaTransferBottom(old, fresh, oldVertices, mtype);
    }

    PList_delete(fresh);
    PList_delete(old);
    PList_delete(oldVertices);
    PList_delete(newVertices);
    PList_delete(oldEdges);
    PList_delete(newEdges);
    PList_delete(oldFaces);
    PList_delete(newFaces);
}


void
phastaTransferSIMBottom(pPList Vertices, pPList Edges, pPList Faces, pPList Regions, int ndof, modType mtype) {

     pVertex vertex, vert1, vert2;
     if(PList_size(Regions)==0 && PList_size(Faces)>1) {  //F_refine or E_refine
       
        for(int i=0; i<PList_size(Vertices); i++) {
           vertex = (pVertex)PList_item(Vertices, i);
           pPList VertCon;
           pVertex vFound[2];

           VertCon = V_vertices(vertex);
           
           int iFound=0;
           for(int j=0; j<PList_size(VertCon); j++) {
              pVertex otherVert;
              otherVert = (pVertex)PList_item(VertCon,j);
              pPList VertFaces;
              VertFaces = V_faces(otherVert);            
              
              int count = 0;
              for(int k=0; k<PList_size(VertFaces); k++) {
                 pFace otherFace;
                 otherFace = (pFace)PList_item(VertFaces, k);

                 if(PList_inList(Faces, otherFace))
                   count++; 
              }
              
              double *Data;
              if(count==1) { //found the vertex for edge split
                 vFound[iFound] = otherVert;
                 iFound++;
              }  
              PList_delete(VertFaces); 
           }
          
           if(iFound==2) { //edge split

              printf("found edge split. yo!\n");

              double xyz[3], node1[3], node2[3];
              double dist1, dist2, distance, xi1, xi2;
              double *Data1, *Data2;
              double* q = new double[ndof];

              V_coord(vertex, xyz);
              V_coord(vFound[0], node1);
              V_coord(vFound[1], node2);
              dist1 = sqrt(dist(node1, xyz));
              dist2 = sqrt(dist(node2, xyz));
              distance = dist1+dist2;
              xi1=dist2/distance;
              xi2=dist1/distance;

              if(!EN_getDataPtr((pEntity)vFound[0], phasta_solution, (void**)&Data1)) {
                 printf("error in callback\n");
//                 exit(0);
              }
              if(!EN_getDataPtr((pEntity)vFound[1], phasta_solution, (void**)&Data2)) {
                 printf("error in callback\n");
//                 exit(0);
              }

              for(int i=0; i<ndof; i++) q[i]=0.0;

              for(int i=0;i<ndof;i++)
                  q[i] = Data1[i]*xi1+ Data2[i]*xi2;
          }    
          else {
             printf("not an edge split :(\n");
          }

           PList_delete(VertCon);
        }
     }

}
#endif

#ifdef FMDB
void
phastaTransferTopSCOREC(pPList oldEnts, pPList newRgn, void *userData, modType mtype, pEntity ent){

    if(isEdgeBasedTransfer && mtype==F_REFINE)
        return;
    if(!isEdgeBasedTransfer && (mtype==E_REFINE || mtype==F_REFINE))
        return;

    CBcounter++;

    pPList fresh;
    pPList old;
    pRegion  rgn;
    pEntity entity;

    pVertex deleteVertex;
    int i,j,k,deletedVerts;
    
    int oldEntsSize = PList_size(oldEnts);
    int freshEntsSize = PList_size(newRgn);

    // old - holds old regions (to be deleted)
    old = PList_new();
    // freash - holds new regions (created)
    fresh = PList_new();

    // VerticestoHandle - hold vertices to be attached solution or deleted
    pPList VerticestoHandle = PList_new();
    if(isEdgeBasedTransfer)
      PList_append(VerticestoHandle,(void*)ent);
    else
      if(mtype==ECOLAPS || mtype==SPLTCLPS)
       if(EN_type(ent)==0) // just to be safe
          PList_append(VerticestoHandle,(void*)ent);

    for(i=0; i<oldEntsSize; i++){
        entity = (pEntity)PList_item(oldEnts,i);
        PList_append(old,entity);
    }

    for(i=0; i<freshEntsSize; i++){
        entity = (pEntity)PList_item(newRgn,i);
        PList_append(fresh,entity);
    }
    
    if(isEdgeBasedTransfer) //fresh is not used in ESPLIT and E_REFINE
       phastaTransferBottomE(old, fresh, VerticestoHandle, mtype);
     else
       phastaTransferBottom(old, fresh, VerticestoHandle, mtype);    

    PList_delete(fresh);
    PList_delete(old);
    PList_delete(VerticestoHandle);

}

void 
phastaTransferBottomE(pPList parent, pPList fresh, pPList VtxstoHandle, modType mtype) {

    pEntity ent;
    pEdge edge;
    pVertex vtx;
    double xietazeta[2], xyz[3];
    int numVarsTot, ndisp = 0;
    double* solution, *data;

    if(DisplacementMigration)
      ndisp = 3;
    numVarsTot = numVars+ndisp;

  if(mtype==E_REFINE || mtype==ESPLIT){

    if(DisplacementMigration)
      ndisp = 3;
    numVarsTot = numVars+ndisp;

    for(int i=0; i<PList_size(VtxstoHandle);i++) { //loop over the new vertices
       vtx = (pVertex)PList_item(VtxstoHandle, i);
       if(!EN_getDataPtr((pEntity)vtx, phasta_solution, (void**)&data)) { //no solution attached yet
           V_coord(vtx, xyz);
           edge = (pEdge)PList_item(parent, i);
           inverseMapE(edge, xyz, xietazeta);
           solution = (double *)malloc(sizeof(double)*numVarsTot);   
           solution = InterpolateSolutionE(edge, xietazeta, numVarsTot, modes);
          
           EN_attachDataPtr((pEntity)vtx, phasta_solution, (void*) solution);
       }
     }
    }
    
   
   //take care of the centroid point created by R_REFINE
   else if(mtype==R_REFINE){
    vtx=(pVertex)PList_item(VtxstoHandle, 0);
    if(vtx && !EN_getDataPtr((pEntity)vtx, phasta_solution, (void**)&data)) { //no solution attached yet
          solution = (double *)malloc(sizeof(double)*numVarsTot);
          for(int i=0; i<numVarsTot;i++) solution[i] = 0.0;
          pRegion rgn = (pRegion)PList_item(parent,0);
          pPList vertices = R_vertices(rgn, 1);
          int numVtx = PList_size(vertices);
          for(int j=0; j<numVtx; j++) {
             pVertex oldvtx = (pVertex)PList_item(vertices,j);
             if(!EN_getDataPtr((pEntity)oldvtx, phasta_solution, (void**)&data)){
                printf("Error in callback function, no solution attached to an old vertex \n");
                exit(-1);
             }
             for(int i=0; i<numVarsTot; i++)
                solution[i] += data[i];
           }
           PList_delete(vertices);

           for(int i=0; i<numVarsTot; i++)
              solution[i] = solution[i]/numVtx;

           EN_attachDataPtr((pEntity)vtx, phasta_solution, (void*)solution);
      }
   }     
   // clean the solution field attached to the vertices which will be deleted in the mesh adapt
   else if(mtype==ECOLAPS || mtype==SPLTCLPS){
     for(int i=0; i< PList_size(VtxstoHandle);i++) {
       pEntity ent = (pEntity)PList_item(VtxstoHandle, i);
       double *data;
       if(EN_getDataPtr((pEntity)ent, phasta_solution, (void**)&data)) {
          delete [] data;
          EN_deleteData((pEntity)ent, phasta_solution);
      }
    }
   }
} 
#endif

void
phastaTransferBottom(pPList old, pPList fresh, pPList oldVertices, modType mtype){

    pPList verticesO;    
    pPList verticesN;
    pRegion  parent_region=0;
    pEntity ent;
    pVertex v;
    pPList newVertices;

    int count = 0;
    int counter = 0;
    int numFoundOldSolution;
    int i,j,k;
    numFoundOldSolution = 0;
    double vol, volTol = 1.e-14;
    double xyz[3], xietazeta[3];
    double* nodalSolution;
    double coords[3];
//     FILE  *logfile; 
    char stringName[255];
    double* averageSolution;
    double* avSol;
    double* solution;

    int numVarsTot, ndisp = 0;
    if(DisplacementMigration)
        ndisp = 3;
    numVarsTot = numVars+ndisp; 
    ////////////////////////////////////////////////////////////////////////////////////////
    // more than 1 region in old
    ////////////////////////////////////////////////////////////////////////////////////////


    if ( PList_size( old ) > 1 ) {
        
        count = PList_size( old );

//        printf("\nsize of old ent list %d  \n",count);

/*          printf("\nMesh mod type is %d and sizeof(old ent list) >1\n",mtype); */
//        PList_printx(old); 
//           PList_printx(fresh); 
        // get all vertices that have some solution on them and build average
        // set this average to all NEW vertices of the NEW regions
        averageSolution = new double [ numVarsTot ];
        for( i =0; i <numVarsTot ; i++ ){
        averageSolution[i]=0.0;
        }
          
        if(mtype==MeshModType_EColapse) 
            printf("\nhere\n");      

        for( i =0; i < count ; i++ ){
             
            verticesO = R_vertices((pRegion)PList_item(old , i ),1);

/*              printf("\n counter is  %d\n",counter); */
/*              printf("\ngot %d old vertices of region %d\n",PList_size(verticesO),i); */

            for( j =0; j < PList_size(verticesO) ; j++ ){

                v = (pVertex)PList_item(verticesO , j );
//                printf("\n got ONE old vertex num %d of region %d\n",j,i );


                if( EN_getDataPtr( ( pEntity )v, phasta_solution, (void
                                                                    **)&nodalSolution) != 0){
                    for( k =0; k <numVarsTot ; k++ ){
                        averageSolution[k] += nodalSolution[k];
                    }
                    numFoundOldSolution++;
//                    printf("\n found solution on old vertex \n" );
                }
            }

	    PList_delete(verticesO);

            counter++;
        }
        if(numFoundOldSolution != 0){
            for( i =0; i <numVarsTot ; i++ ){
                averageSolution[i] = averageSolution[i]/numFoundOldSolution;
            }
        }
        else{
            printf("\n all vertices in old region list have NO solution \n" );
            exit(1);
        }
        // assign the average to ALL new vertices
        for( i =0; i < PList_size( fresh ) ; i++ ){

//            printf("\n got vertices of a new region \n" );

            verticesN = R_vertices((pRegion)PList_item(fresh , i ),1);



            for( j =0; j < PList_size(verticesN  ) ; j++ ){

                v =  (pVertex)PList_item(verticesN , j );


                if(! EN_getDataPtr( ( pEntity )v, phasta_solution, (void
                                                                    **)&nodalSolution)){
                
                    // provide memory for new solution
		    // avSol= (double *)malloc(sizeof(double)*numVars );
		    avSol = new double[numVarsTot];
                    for( i =0; i <numVarsTot ; i++ ){
                        avSol[i]=averageSolution[i];
                    }
//                    printf("\n found vertex in new regions that got assigned a sol \n" );
                    EN_attachDataPtr( (pEntity)v, phasta_solution, (void *)
                                      avSol ); 
                }
            }

	    PList_delete(verticesN);
        }
        delete [] averageSolution;


        
    }// PList_size(old) > 1 
    ////////////////////////////////////////////////////////////////////////////////////////
    // exactly 1 region in old list ==> it is  a split operation
    ////////////////////////////////////////////////////////////////////////////////////////
    else if( PList_size( old ) == 1){
//        printf("\n  1 region in old, split type is %d\n",mtype );
        
        parent_region = (pRegion)PList_item( old, 0 );
        
        // split operation
        // retrieve new vertices
        // get all vertices on old and new element(s)
        verticesO = R_vertices(parent_region,1);
        
        newVertices = PList_new();

        for( i =0; i < PList_size( fresh ) ; i++ ){

            verticesN = R_vertices((pRegion)PList_item(fresh , i ),1);

            for( j =0; j < PList_size(verticesN  ) ; j++ ){
                
                if(!PList_inList(newVertices,PList_item(verticesN , j ))){
                     PList_append(newVertices,PList_item(verticesN , j ));
                }
            }

	    PList_delete(verticesN);
        }         
        // find the vertices which are not in the old list (-> a new vertex)
        for( i =0; i < PList_size(newVertices  ) ; i++ ){
            
            v =  (pVertex)PList_item(newVertices , i );
            V_coord(v,xyz);
             
            if(! PList_inList(verticesO,v)){
                
                //if there is not a solution attached already
                if(! EN_getDataPtr( ( pEntity )v, phasta_solution, (void
                                                                    **)&nodalSolution)){

                    // provide  memory for new solution
                    solution= (double *)malloc(sizeof(double)*numVarsTot );

                    // do a solution interpolation for that vertex
                    // calculate the parametric co-ordinates of the current vertex in the
                    // space of the parent region. ( also use these to make sure that the
                    // vertex lies inside the given region )
                    // inversemap implemented in InverseMap.cc
                    if ( !inverseMap( parent_region, xyz, xietazeta ) ) {	
                        printf("Could not locate point inside the parent region: \n");
                        R_info(parent_region);
                        V_info(v);
                        exit( -1 );
                    }
                    // interpolate the solution from the original region onto the finer one 
                    // use appropriate shapeFuns
                    solution = InterpolateSolution( parent_region,
                                                    xietazeta,
                                                    numVarsTot,
                                                    modes);
                    
                    // attach it to the vertex (if not done yet)
                    // otherwise do nothing
                    EN_attachDataPtr( (pEntity)v, phasta_solution, (void *)
                                      solution );
                    
//                    printf("\n info: assigning a solution val to a new vertex: \n" );
/*                      V_coord( v ,coords ); */
/*                      for(j=0; j<3; j++){ */
/*                          //fprintf(logfile,"%f\n",coords[j]); */
/*                      } */
                    //fprintf(logfile,"\n");
                }//if(! EN_getDataPtr(  
                else{
                    // a solution already exists at that vertex: do nothing
                }
                
            }// if(! PList_inList(v
        }//for(...
       

        PList_delete(verticesO);
        PList_delete(newVertices);
    }
    else{//  size of old regions list =0
//        printf("\n  in callback.c: size of old regions list =0 \n" );
    }


    for( i =0; i < PList_size(oldVertices) ; i++ ){
        ent = (pEntity)PList_item(oldVertices,i);
        double* data;
        if(!EN_getDataPtr( ( pEntity )ent, phasta_solution, (void
							     **)&data)) {
	  printf("\nError in callback.cc : solution NOT attached to an old vertex\n");
	  exit(1);
	}
        
        delete [] data;
	EN_deleteData((pEntity)ent,phasta_solution);
    }
    
/*      printf("\n  in callback.c:exit callback function \n" ); */
/*      printf("[%d] exiting Callback.cc procSize (MY) : %d \n",PMU_rank(),myProcSize()); */
}

extern int delDblArraycounter;
#ifdef SIM
int delDblArray(void *ent, pAttachDataId id, int cb, void** data, void* c) {
#endif
#ifdef FMDB
int delDblArray(pAttachableData ad, pAttachDataId id, int cb, void** data, void* c) {
#endif
  delDblArraycounter++;

  double *dblData = (double *)(*data);
  delete [] dblData;

  return 1;
}

extern int delDblcounter;
#ifdef SIM
int delDbl(void *ent, pAttachDataId id, int cb, void** data, void* c) {
#endif
#ifdef FMDB
int delDbl(pAttachableData ad, pAttachDataId id, int cb, void** data, void* c) {
#endif
  delDblcounter++;

  double *dblData = (double *)(*data);
  delete dblData;

  return 1;
}


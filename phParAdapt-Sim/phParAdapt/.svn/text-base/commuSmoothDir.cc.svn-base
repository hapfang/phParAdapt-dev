#include "phParAdapt.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"
#include "Eigen.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalDirectionID;
extern pMeshDataId numSurroundNodesID;


// communicates missing contributions over procs/=parts
// each part p's  vertex has
// average nodal Dir via nodalDirectionID
//
// data are gathered via
// nodalDirectionID
// numSurroundNodesID
// The commu attaches the combined data  via nodalDirectionID
// ( delete old ones on part bdry)
void
commuSmoothDir(pParMesh pmesh, pMesh mesh)
{
    int mult = sizeof(P_int)/sizeof(int);
    int *ns = new int [ PMU_size() ];
    int *nr = new int [ PMU_size() ];
    void ***se = new void** [ PMU_size() ];
    void ***re = new void** [ PMU_size() ];

    pEntity ent;
    double** dsend = new double*[PMU_size()];
    double** drecv = new double*[PMU_size()];

    for(int i=0; i< PMU_size(); i++)  {
        ns[i] = 0;
    }
    pBdryProcIter myBdryIter = PM_bdryProcIter(pmesh,0);
    pEntity en;
    pEntOrig eo;
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {

       if(!EN_isOwnerProc(en)) { //only nonowned are going to send
         
         int OwnerRank;
         OwnerRank = EN_ownerProc(en);

          // increment number of copies shared between this proc and proc
          // with rank of owned
          ns[OwnerRank]++;
      } //!EN_isOwnerProc
   }
    BdryProcIter_reset(myBdryIter);

    for (int i=0; i < PMU_size(); i++)  {
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[10*ns[i]];
        ns[i] = nr[i] = 0;
    }

    while(en = BdryProcIter_next(myBdryIter)) {

       if(!EN_isOwnerProc(en)) { //only nonowned are going to send
         
         eo = EN_original(en);
         int OwnerRank;
         OwnerRank = EN_ownerProc(en);
          // fill send array to proc of rank vCopyGid with the
          // (ns[vCopyGid])-th copy (among all copies of verts shared
          // between this proc and proc OwnerRank)
          se[OwnerRank][ ns[OwnerRank] ] = EntOrig_ent(eo);
                
          double* nodalDir;
          // get the local data
          if(!EN_getDataPtr(en, nodalDirectionID, (void**)&nodalDir)){
               cout<<"\nerror in commuSmoothDir: no nodalDir data attached to vertex\n";
               exit(0);
          }

          double* numSurroundNodes;
          if(!EN_getDataPtr(en, numSurroundNodesID,
                (void**)&numSurroundNodes)){
             cout<<"\nerror in commuSmoothDir: no numSurroundNodes data attached to vertex\n";
             exit(0);
          }
          // numSurroundNodes
          dsend[OwnerRank][ 10*ns[OwnerRank]  ] =   numSurroundNodes[0];
                
          // nodal Dir
          for(int k=0;k<9;k++){
             dsend[OwnerRank][ 10*ns[OwnerRank] +1+k ] = nodalDir[k];
          }
               
          // increment number of copies shared between this proc and proc
          // with rank vCopyGid
          ns[OwnerRank]++;
#ifdef FMDB
        EntOrig_delete(eo);
#endif
       } //isOwned
    }
    BdryProcIter_reset(myBdryIter);

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for (int i=0; i < PMU_size(); i++){
        re[i] = new void* [ nr[i] ];
        drecv[i] = new double[10*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);    
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, 10); 

    double neweOwn[3][3];
    double cureOwn[3][3];
    double eNorm[3][3];
    double dprod;

    // retrieve communicated data and globalize
    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i with this proc being owner of those
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];
            if(!EN_isOwnerProc(ent)) {
               printf("Something is wrong, I am not owner, no...\n");
               exit(0);
            }   

            // get nodal Dir (has been commuicated previous to this call)
            // and smooth = average
            double* currentDir;
            if(!EN_getDataPtr(ent, nodalDirectionID,
			      (void**)&currentDir)){
                cout<<"\nerror in commuSmoothDir globalize: no currentDir data attached to vertex\n";
                exit(0);
            }
            double* newDir = new double[9];  
            for(int k=0;k<9;k++){
                newDir[k]= drecv[i][10*j+1+k];
            }

            for(int iRow=0; iRow<3; iRow++) {
               for(int iCol=0; iCol<3; iCol++) {
                  neweOwn[iRow][iCol]=newDir[iRow*3+iCol];
                  cureOwn[iRow][iCol]=currentDir[iRow*3+iCol];
               }
               dprod = dotProd(neweOwn[iRow], cureOwn[iRow]);
               if(dprod < 0.0) {
                  for(int iCol=0; iCol<3; iCol++) {
                     newDir[iRow*3+iCol]=cureOwn[iRow][iCol]-neweOwn[iRow][iCol];
                  }
               }
               else {
                  for(int iCol=0; iCol<3; iCol++) {
                     newDir[iRow*3+iCol]=cureOwn[iRow][iCol]+neweOwn[iRow][iCol];
                  }   
               }   //if else dot prod
            } //iRow  

            delete [] currentDir;
            EN_deleteData(ent, nodalDirectionID);
            EN_attachDataPtr( ent, nodalDirectionID, (void *)
                              newDir );      

        } //j loop
    } // iloop

    MPI_Barrier(MPI_COMM_WORLD);

    //now to copy from owner to remote copies
    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i with this proc being owner of those
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];
           
            if(!EN_isOwnerProc(ent)) {
               printf("Again something is wrong, I am not owner, no...\n");
               exit(0);
            }

           pEntCopies vCopies = EN_copies(ent);
          // how many of them on OTHER apts/=procs ?
          int vCopiesSize = EntCopies_size(vCopies );

          // fill the send array with ents/vertices
           int Index = -1;
          for(int k=0; k< vCopiesSize; k++){

           // global id of partition/proc on which en's n'th copy lies.
           int vCopyGid = EntCopies_gid(vCopies,k);

           if (vCopyGid == i)  {
              Index = k; //kth copy is on proc i
           
           // fill send array to proc of rank vCopyGid with the
           // (ns[vCopyGid])-th copy (among all copies of verts shared
           // between this proc and proc vCopyGid)
            re[i][j] = EntCopies_ent(vCopies,Index);

              double* nodalDir;
            // get the local data
             if(!EN_getDataPtr(ent, nodalDirectionID, (void**)&nodalDir)){
               cout<<"\nerror in commuSmoothDir: no nodalDir data attached to vertex\n";
               exit(0);
             }
            double* numSurroundNodes;
            if(!EN_getDataPtr(ent, numSurroundNodesID,
                  (void**)&numSurroundNodes)){
               cout<<"\nerror in commuSmoothDir: no numSurroundNodes data attached to vertex\n";
               exit(0);
            }
            // numSurroundNodes
             drecv[i][ 10*j  ] =   numSurroundNodes[0];
                
            // nodal Dir
            for(int m=0;m<9;m++){
               drecv[i][ 10*j+1+m ] = nodalDir[m];
            }
               
           } //if found ent copy
           
         } //loop over number of copies
       } //j loop
     } //i loop

    PMU_commuInt(nr, ns);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for (int i=0; i < PMU_size(); i++){
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[10*ns[i]];
    }
    
    //send total value from owner to remote copies
    PMU_commuArr((void**)re, nr, (void**)se, ns, MPI_INT, mult);    
    PMU_commuArr((void**)drecv, nr, (void**)dsend, ns, MPI_DOUBLE, 10); 

    //receive values on this remote copy from the owner
    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i with this proc being owner of those
        for (int j=0; j < ns[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) se[i][j];

            double* currentDir;
            if(!EN_getDataPtr(ent, nodalDirectionID,
			      (void**)&currentDir)){
                cout<<"\nerror in commuSmoothDir globalize: no currentDir data attached to vertex\n";
                exit(0);
            }

            double* newDir = new double[9];  
            for(int k=0;k<9;k++){
                newDir[k]= dsend[i][10*j+1+k];
            }

            delete [] currentDir;
            EN_deleteData(ent, nodalDirectionID);
            EN_attachDataPtr( ent, nodalDirectionID, (void *)
                              newDir );      

         } //j loop
    } //i loop

/*
    //now have to copy the data from owner to the remote copies
    while(en = BdryProcIter_next(myBdryIter)) {

       if(EN_isOwnerProc(en))  { 
       pEntCopies vCopies = EN_copies(en);

   send// how many of them on OTHER apts/=procs ?
       int vCopiesDir = EntCopies_size(vCopies );
       
       // fill the send array with ents/vertices
       for(int i=0 ; i< vCopiesDir; i++){
          
          // global id of partition/proc on which en's n'th copy lies.
          int vCopyGid = EntCopies_gid(vCopies,i);
          
          if (EntCopies_gid(vCopies,i) != PMU_rank())  {
             
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

    for (int i=0; i < PMU_size(); i++)  {
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[10*ns[i]];
        ns[i] = nr[i] = 0;
    }
    while(en = BdryProcIter_next(myBdryIter)) {
       if(EN_isOwnerProc(en) ) {
          pEntCopies vCopies = EN_copies(en);

        // how many of them on OTHER apts/=procs ?
        int vCopiesDir = EntCopies_size(vCopies );

         // fill the send array with ents/vertices
         for(int i=0 ; i< vCopiesDir; i++){

            // global id of partition/proc on which en's n'th copy lies.
            int vCopyGid = EntCopies_gid(vCopies,i);
              
            if (EntCopies_gid(vCopies,i) != PMU_rank())  {
               
               // fill send array to proc of rank vCopyGid with the
               // (ns[vCopyGid])-th copy (among all copies of verts shared
               // between this proc and proc vCopyGid)
               se[vCopyGid][ ns[vCopyGid] ] = EntCopies_ent(vCopies,i);
               
             double* nodalDir;
             // get the local data

             if(!EN_getDataPtr(en, nodalDirectionID, (void**)&nodalDir)){
                
                cout<<"\nerror in commuSmoothDir: no nodalDir data attached\n";
                   exit(0);
             }
             double* numSurroundNodes;
             if(!EN_getDataPtr(en, numSurroundNodesID,
                      (void**)&numSurroundNodes)){
                cout<<"\nerror in commuSmoothDir: no numSurroundNodes data\n";
                   exit(0);
             }
             // numSurroundNodes
             dsend[vCopyGid][ 4*ns[vCopyGid]  ] =   numSurroundNodes[0];
             
             // nodal Dir
 
             for(int k=0;k<3;k++){
                dsend[vCopyGid][ 4*ns[vCopyGid] +1+k ] = nodalDir[k];
             }
             
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
*/
/*      
    // average the data attached via nodalDirectionID and numSurroundNodesID
    // over to data attched via nodalDirectionID again
    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
        

        double* currentNumSurrNodes;
        if(!EN_getDataPtr(en, numSurroundNodesID,
			  (void**)&currentNumSurrNodes)){
            printf("\n[%i]error in commuSmoothDir average: no currentNumSurrNodes data attached to vertex\n",PMU_rank());
            exit(0);
        }
        double* nodalDir;
        if(!EN_getDataPtr(en, nodalDirectionID,
			  (void**)&nodalDir)){
            printf("\n[%i]error in commuSmoothDir average: no nodalDir data attached to vertex\n",PMU_rank());
            double c[3];
            V_coord((pVertex)en,c);
            printf("coords %f %f %f\n",c[0],c[1],c[2]);
            exit(0);
        }

        double* newNodalDir = new double[9];
        for(int k=0;k<9;k++){
	  if(currentNumSurrNodes[0])
            newNodalDir[k] = nodalDir[k]/currentNumSurrNodes[0];
	  else
	    newNodalDir[k] = nodalDir[k];
        }
        delete [] nodalDir;
        EN_deleteData(en, nodalDirectionID);
        EN_attachDataPtr( en, nodalDirectionID, (void *)
                          newNodalDir) ;


#if  ( defined  DEBUG )
//             printf("\nattaching (globBDR) smooth hessian:\n %f %f %f\n",newNodalDir[0],newNodalDir[1],newNodalDir[2]);
//             printf(" %f %f %f\n",newNodalDir[3],newNodalDir[4],newNodalDir[5]);
//             printf("for vertex:");
//             double c[3];
//             V_coord((pVertex)en,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 

    }//while(en = BdryProcIter_n
    BdryProcIter_delete(myBdryIter);
*/
    for (int i=0; i < PMU_size(); i++)  {
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


    // loop over all verts that are NOT on part. bdry.
    // and shovel local data into gobal data
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;
    
    int Count=0;
    while(vertex = VIter_next(vIter)) {

        double* nodalDir; 
//        if(! EN_isOnPartBdry((pEntity)vertex)){

            if(!EN_getDataPtr((pEntity)vertex, nodalDirectionID,
			      (void**)&nodalDir)){
                    cout<<"\nerror in INTR commuSmoothDir: no nodalDir data attached to vertex\n";
                    exit(0);
            }
            
            double*  numSurrVerts;
            double* smoothDir = new double[9];
            for(int i=0; i<9; i++){
                smoothDir[i]= 0.0  ;
            }
            if(!EN_getDataPtr((pEntity)vertex,numSurroundNodesID ,
			      (void**)&numSurrVerts)){
                    cout<<"\nerror in INTR commuSmoothDir: no numSurrVerts data attached to vertex\n";
                    exit(0);
            }

            for(int iRow=0; iRow<3;iRow++) {
               for(int iCol=0; iCol<3; iCol++) {
                  neweOwn[iRow][iCol]=nodalDir[iRow*3+iCol];
               }
               normVt(neweOwn[iRow], eNorm[iRow]);
            }

            //check orthogonality of the resultant vectors            
/*            int factor;
            if(!checkUnitaryOthoganal(eNorm, factor)) {
               dprod=dotProd(eNorm[0], eNorm[2]);
               for(int iCol=0; iCol<3; iCol++) {
                  eNorm[2][iCol]=eNorm[2][iCol]-dprod*eNorm[0][iCol];
               }
               normVt(eNorm[2], eNorm[2]);
               crossProd(eNorm[0], eNorm[2], eNorm[1]);
            }
            if(!checkUnitaryOthoganal(eNorm, factor)) {
               Count++;
               printf("still not orthogonal\n");
            }
*/            
            for(int iRow=0; iRow<3; iRow++) {
              for(int iCol=0; iCol<3; iCol++){
                smoothDir[iRow*3+iCol]=eNorm[iRow][iCol];
              }
            }

            // attach these values
            delete [] nodalDir;
            EN_deleteData((pEntity)vertex, nodalDirectionID);

            EN_attachDataPtr( (pEntity)vertex, nodalDirectionID, (void *)
                              smoothDir);
        }
//    }// while loop over non bdry vertices
    VIter_delete(vIter);
//    printf("number of non orthogonal in commu %i\n", Count);
}

#ifdef __cplusplus
}
#endif

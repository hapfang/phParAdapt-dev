#ifdef PARALLEL

/************************************************************************/
/* This function sets up the interprocessor communication tasks.      */
/**********************************************************************/
#include <stdio.h>
#include "func.h"
#include "parallel.h"
#include <map>
#ifdef SIM
#include "SimPartitionedMesh.h"
#endif
#include "phParAdapt.h"

#ifdef __cplusplus
extern "C" {
#endif

pMeshDataId REIID;
pMeshDataId REJID;
pMeshDataId IPERID;
extern pMeshDataId MYCTID;
extern int numTotParts;
extern int numParts;

// called in mdb2phasta
void setupGlobalTasks(vector<pParMesh> pmesh, globalInfo **ginfo)
{
  pEntity ent;
  pEntity em;
  pVertex vertex;
  pMatch ment;
  pEntOrig eo;
  pEntCopies ec;

  PBEntPartIter ei;
  int ipart;
  pParMesh mesh;
  globalInfo *info;

  
  REIID =  MD_newMeshDataId( "REI+"); 
  REJID =  MD_newMeshDataId( "REJ+"); 
  IPERID =  MD_newMeshDataId( "IPER");

  void *tmps, *tmpm;
  int i, j, k, etype, gid, pid, mpid, mult = sizeof(P_int)/sizeof(int);

  int **ns = new int* [ numParts ];
  int **nr = new int* [ numParts ];
  void ****se = new void*** [ numParts ];
  void ****re = new void*** [ numParts ];
  map<int,int> *nsmaps = new map<int,int>[numParts];
  map<int,int> *nrmaps = new map<int,int>[numParts];
  map<int,int>::iterator nsiter, nriter;
  int *neigborS = new int[numParts];
  int *neigborR = new int[numParts];
  for(i=0; i<numParts; i++) {
      neigborS[i] = 0;
      neigborR[i] = 0;
      ns[i] = 0;
      nr[i] = 0;
      se[i] = 0;
      re[i] = 0;
  }

  int PeriodicFlag = 0, PeriodicFlagS = 0;;

  map<int,int>::iterator iter1, iter2;

  mesh = pmesh[0];

  for(ipart =0 ;ipart<numParts;ipart++){ //loop over all the parts on each proc
      
   info = ginfo[ipart];
   map<int,int> numS;
   
   pBdryPartIter BIter=PM_bdryPartIter(mesh, ipart, 0);

   while(vertex = (pVertex)BdryPartIter_next(BIter)){

//follwing starts a quest to determine if the entity owned on this "part" or not. This is more complex than just determining if it is owner by the process. This is due to how simmetrix and fmdb handle things differently

      int ownID = EN_ownerProc(vertex);
      int partid;
       if(ownID == PMU_rank() ) { //owned by this process 
          int ownpID = EN_isOwnerPart(vertex, ipart); //cant call this for other process
          if(ownpID) { // owned by this part (ipart)
            int globid = EN_gid(vertex);
            partid = PMU_part(globid);
          }
          else { //owned by this proc but different part
             pEntOrig eo = EN_original(vertex);
             pEntity vOrg = EntOrig_ent(eo);
             int orgid = EntOrig_gid(eo);
             partid = PMU_part(orgid);
             EntOrig_delete(eo);
          }
       }
       else { // owner by different proc
           pEntOrig eo = EN_original(vertex);
           pEntity vOrg = EntOrig_ent(eo);
           int orgid = EntOrig_gid(eo);
           partid = PMU_part(orgid);
           EntOrig_delete(eo);
       }
       
       int neigborID = ownID*numParts+partid;
       
//      int neigborID = EN_ownerProc(vertex);  //previous
       if(neigborID!=PMU_rank()*numParts+ipart){  //owned by others, need receive
           if(nrmaps[ipart].find(neigborID)==nrmaps[ipart].end()) //ID is not in the map yet
               nrmaps[ipart][neigborID] = neigborR[ipart]++;
       }

       else { //owned by current processor, need send          
           pEntCopies ec = EN_copies(pEntity(vertex)); //find remote copies
           int entSize = EntCopies_size(ec);
           for(int iec=0; iec<entSize;iec++){
               int ecgid = EntCopies_gid(ec, iec); //find gid for
               int proc = PMU_proc(ecgid); //procid 
               int partid = PMU_part(ecgid); //local partid

/* previous api commented below               
               if(nsmaps[ipart].find(ecgid)==nsmaps[ipart].end()) {
                   nsmaps[ipart][ecgid] = neigborS[ipart]++;
                   numS[ecgid]=1; 
               }
               else
                   numS[ecgid]++;
*/
//trying to deal with multiple parts here for simmetrix
               if(nsmaps[ipart].find(proc*numParts+partid)==nsmaps[ipart].end()) {
                   nsmaps[ipart][proc*numParts+partid] = neigborS[ipart]++;
                   numS[proc*numParts+partid]=1; 
               }
               else
                   numS[proc*numParts+partid]++;
               
           }
       }
   }

   BdryPartIter_delete(BIter);
  
   if(neigborS[ipart]){
       ns[ipart] = new int [neigborS[ipart]];
       se[ipart] = new void** [neigborS[ipart]];
       map<int, int>::iterator numSiter = numS.begin();
       for(i=0;i<neigborS[ipart];i++,numSiter++){
          ns[ipart][i] = 0;
//          se[ipart][i] = new void*[2*numS[numSiter->second]];
          se[ipart][i] = new void*[2*numSiter->second];
       }
   }
   
   if(neigborR[ipart]){
       nr[ipart] = new int[neigborR[ipart]];
       re[ipart] = new void**[neigborR[ipart]];
       for(i=0;i<neigborR[ipart];i++) {
           nr[ipart][i] = 0;
           re[ipart][i] = 0;
       }
   }


  /* first take care of periodic pairs in info that we made in localInfo */

  for (tmps=0, tmpm=0; (ent = (pEntity) PList_next(info->perSlv, &tmps))
                    && (ment = (pMatch) PList_next(info->perMst, &tmpm));)  {
    /* the master is offproc, we made sure of it in localInfo */

      // PMU_proc(), gives rank, given global partition id (gid)
      // Match_gid(), Returns global id of partition on which entity in ment
      // lies.
      if(PMU_rank==0 && ipart==0)
          std::cout<<"Periodic conditition is not supported yet \n";

      pid = PMU_proc(Match_gid(ment));

      addSendSegment(info, ent, Match_ent(ment), pid);

  }
  
  for (etype = 0; etype < 3; etype++)  {
      if ((etype==1 && !info->edgeson) || (etype==2 && !info->faceson)){

          continue;
      }
      // neiter (etype==1 && !info->edgeson) || (etype==2 && !info->faceson)
      ei = PM_bdryPartIter(mesh, ipart, etype);

      while (ent = BdryPartIter_next(ei))  {

          if (isPeriodic(ent, &ment))  {
              if(PMU_rank==0 && ipart==0)
                  std::cout<<"Periodic conditition is not supported yet \n";
              PeriodicFlag = 1;
              mpid = PMU_proc(Match_gid(ment));  
              ec = EN_copies(ent);
              for (i=0; i < EntCopies_size(ec); i++)  {
                  
                  // if the owner proc of the current copy is NOT PMU_rank(),
                  // i.e. it is off-proc
                  // fill up the send array with the copie's (off-proc) pointer
                  if ((gid = PMU_proc(EntCopies_gid(ec,i))) != PMU_rank())  {
                      pid = nsmaps[ipart][gid];
                    if (se[ipart][pid]) {  // pointer exists
                      se[ipart][pid][3*ns[ipart][pid]] = EntCopies_ent(ec, i);
                      se[ipart][pid][3*ns[ipart][pid]+1] = Match_ent(ment);
                      se[ipart][pid][3*ns[ipart][pid]+2] = (void*) mpid;
                      ns[ipart][pid]++;
                    } else {
                      printf("[%2d] setupGlobal: Writing to unallocated se[%d] \n",PMU_rank(),pid);
                      exit(1);
                    }
                  }
              }
#ifdef FMDB
              EntCopies_delete(ec);
#endif
          }
      }
      BdryPartIter_delete(ei);
  
  }// for (etype ...
  } // loop over the multiple parts for each proc

  // Communicates integer arrays across processes, only in parallel
  // collective call

  MPI_Allreduce(&PeriodicFlag, &PeriodicFlagS, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  

  if(PeriodicFlagS) {
      PMU_commuInt(ns, nr, nsmaps, nrmaps);
      
      for(ipart=0; ipart<numParts;ipart++){
          for (i=0; i < neigborR[ipart]; i++){
              if (nr[ipart][i] > 0) {
                  re[ipart][i] = new void* [ 3*nr[ipart][i] ];
              }
          }
      }
      
      
      // Communicates arrays across processes, only in parallel
      PMU_commuArr((void***)se, ns, (void***)re, nr, nsmaps, nrmaps, MPI_INT, 3*mult);
  
      for(ipart = 0; ipart < numParts; ipart++) {
          for(nriter = nrmaps[ipart].begin();nriter!=nrmaps[ipart].end();nriter++) {
              i = nriter->first;
              k = nriter->second;
              for (j=0; j < nr[ipart][k]; j++)  {
                  ent = (pEntity) re[ipart][k][3*j];
                  EN_attachDataInt(ent, REIID, i+1);  /* +1 is to avoid the 0 ambiguity */
                  EN_attachDataInt(ent, REJID, j);
              }
          }
      }
  } //if there is periodicity

  mesh = pmesh[0];
  for(ipart=0; ipart<numParts;ipart++){
    info = ginfo[ipart];

    for (etype = 0; etype < 3; etype++)  {
      if ((etype==1 && !info->edgeson) || (etype==2 && !info->faceson)){
          
          continue;
      }
      // neiter (etype==1 && !info->edgeson) || (etype==2 && !info->faceson)
      ei = PM_bdryPartIter(mesh, ipart, etype);

      while (ent = BdryPartIter_next(ei))  {
            
            int globid;
//again embark on the quest to determine if this ent is copy or not                    
        int ownID = EN_ownerProc(ent);
        int partid;
        if(ownID == PMU_rank() ) { //owned by a local part
           int ownpID = EN_isOwnerPart(ent, ipart);
           if(ownpID) { // owned by ipart too
             int globid = EN_gid(ent);
             partid = PMU_part(globid);
           }
           else { //some other part than ipart
              pEntOrig eo = EN_original(ent);
              pEntity vOrg = EntOrig_ent(eo);
              int orgid = EntOrig_gid(eo);
              partid = PMU_part(orgid);
              EntOrig_delete(eo);
           }
        }
        else {  //off proc
           pEntOrig eo = EN_original(ent);
           pEntity vOrg = EntOrig_ent(eo);
           int orgid = EntOrig_gid(eo);
           partid = PMU_part(orgid); 
           EntOrig_delete(eo);
        }
        pid = ownID*numParts+partid;

//       pid = EN_ownerProc(ent);  //previous api
       if (EN_isOnPartBdry(ent)&&(pid != PMU_rank()*numParts+ipart))  {  // found a slave
 
              if (EN_getDataInt(ent, REIID, &i))  {
                  int mapid = nrmaps[ipart][i-1];
                  EN_getDataInt(ent, REJID, &j);
                  em = (pEntity) re[ipart][mapid][3*j+1];
#ifdef ibm
                  mpid = (long) re[ipart][mapid][3*j+2];
#else                  
                  mpid = (P_int) re[ipart][mapid][3*j+2];
#endif               

// For MeshSim v5.4, off-proc periodic nodes were not captured in perSlv and perMst above so
// this section was needed to get them into the ilwork array.
// MeshSim V6.0 does classify the off-proc periodic nodes correctly so this addSendSegment will
// result in double accounting for the off-proc periodic nodes.  So will comment
// it out for now.
                  if (mpid != PMU_rank()*numParts+ipart) {
//                      addSendSegment(info, ent, em, mpid);//commented out !

		  }    
                  else
                      EN_attachDataPtr(ent, IPERID, (void*) em);
              }
              else  {
                  int myct, myct2, isown;
                  eo = EN_original(ent);
                  
                  ec = EN_copies(ent);
                  int size = EntCopies_size(ec);
                  pEntity CopEnt, OrgEnt;
                  OrgEnt = EntOrig_ent(eo);
                  globid = EntOrig_gid(eo); 
/*   The following code is not required but keeping for safekeeping 
                  for(int it=0; it<size; it++ ) {
                    CopEnt = EntCopies_ent(ec, it);
                    int gid = EntCopies_gid(ec, it);
                    int partid = PMU_part(gid);
                    int ownid = PMU_proc(gid);
                    if(globid == gid) {
                       OrgEnt = CopEnt;
                    }
                    else continue; 
                  }
*/                  
                  addSendSegment(info, ent, OrgEnt, pid);
                  EntOrig_delete(eo);

              }
          }
      }
      BdryPartIter_delete(ei); 
  }//for etype


  /* setup send stuff */
  for (i=0; i < neigborS[ipart]; i++)  {
    if (se[ipart][i]) {  // only delete memory with established pointers
      delete [] se[ipart][i];
      se[ipart][i] = 0;  // reset pointer to NULL
    } 
    ns[ipart][i] = 0;
  }

  for (i=0; i < neigborR[ipart]; i++)  {
    if (re[ipart][i]) {  // only delete memory with established pointers
      delete [] re[ipart][i];
      re[ipart][i] = 0;  // reset pointer to NULL
    }
    nr[ipart][i] = 0;
  }

   for(nriter=nrmaps[ipart].begin();nriter!=nrmaps[ipart].end();nriter++) {
    i = nriter->first;
    j = nriter->second;
    nr[ipart][j] = info->stask[i] ? info->stask[i]->numSeg : 0;
    if (nr[ipart][j] > 0) {
      re[ipart][j] = new void* [ nr[ipart][j] ];
      for (k=0; k < nr[ipart][j]; k++)
        re[ipart][j][k] = (pEntity) PList_item(info->stask[i]->ments, k);
    }
   }
  }//loop over multiple parts on each proc

  // Communicates integer arrays across processes, only in parallel
  // collective call
  PMU_commuInt(nr, ns, nrmaps, nsmaps);


  for(ipart=0;ipart<numParts;ipart++)
  for (i=0; i < neigborS[ipart]; i++) {
    if (ns[ipart][i] > 0) {
      se[ipart][i] = new void* [ ns[ipart][i] ];
    }
  }


  PMU_commuArr((void***)re, nr, (void***)se, ns, nrmaps, nsmaps, MPI_INT, mult);
  

  for(ipart=0;ipart<numParts;ipart++){
      info = ginfo[ipart];
  
  for(nsiter=nsmaps[ipart].begin();nsiter!=nsmaps[ipart].end();nsiter++) {
      i = nsiter->first;
      k = nsiter->second;
     for (j=0; j < ns[ipart][k]; j++)
         addRecvSegment(info, (pEntity) se[ipart][k][j], i);
  }


  for (i=0; i < neigborS[ipart]; i++)  {
    if (se[ipart][i]) {  // only delete memory with established pointers
      delete [] se[ipart][i];
      se[ipart][i] = 0;  // reset pointer to NULL
    }
  }

  for (i=0; i < neigborR[ipart]; i++)  {
    if (re[ipart][i]) {  // only delete memory with established pointers
      delete [] re[ipart][i];
      re[ipart][i] = 0;  // reset pointer to NULL
    }
  }

  if(neigborS[ipart]){     
      delete [] se[ipart];
      delete [] ns[ipart];
   }
  if(neigborR[ipart]){     
      delete [] re[ipart];
      delete [] nr[ipart];
  }

  }
  
  delete [] se;
  delete [] re;
  delete [] ns;
  delete [] nr;
  delete [] nsmaps;
  delete [] nrmaps;
  delete [] neigborR;
  delete [] neigborS;  
  
  }
  
#ifdef __cplusplus
}
#endif

#endif

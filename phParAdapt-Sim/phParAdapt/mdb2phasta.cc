////////////////////////////////////////////////////////////////////////
//
// This program generates the input files necessary to run ENSA in
// parallel mode. All boundary conditions and communication data
// structures are set up.
//
// C. Whiting
// Fall `97
//
// 
////////////////////////////////////////////////////////////////////////
#include "parallel.h"
#include <string.h>
#include <strings.h>
#include <fstream>
#include <string>
#include <iostream>
#include "func.h"
#include "attachData.h"
#include "ccfunc.h"
#include "BoundaryCondition.h"
#include "phParAdapt.h"

extern int isReorder;
int totNumNodes=0;
int totNumEdges=0;
int totNumFaces=0;
using namespace std;

extern "C" int procSize();
extern "C" int phParAdaptProcSize();

#ifdef SIM
#include "SimModel.h"

#ifdef MODELER_SHAPES
#include "modelerShapes.h"
#include "ShapesModel.h"
#endif


#ifdef DISCRETE
#include "SimDiscrete.h"
#endif

#ifdef DISCRETE
//extern "C" void GM_registerDiscrete();
//  extern "C" void SimDiscrete_start();
//  extern "C" void SimDiscrete_stop();
#endif

#ifdef SIM_PARASOLID
//#include "ParasolidModel.h"
//#include "modelerParasolid.h"
//#include "SimParasolidInt.h"
#include "SimParasolidKrnl.h"
//using namespace SCOREC_ParasolidModel;
#endif

#ifdef ACIS
#include "SimAcisInt.h"
#include "SimAcisKrnl.h"
#endif
#else // not SIM

// using namespace SCOREC_att;
// #include "MeshSim.h"
// #include "MeshSimInternal.h"

// #ifdef MODELER_DISCRETE
// #include "SimDiscrete.h"
// #endif

// #ifdef MODELER_SHAPES
// #include "MeshSimShapes.h"
// #endif

// #ifdef MODELER_PARASOLID
// #include "MeshSimParasolid.h"
// #endif

#endif



// old extern "C" void MA_init();
void printTime(double *eltime, int nt);

/* global variables !!!*/
// dangerous !!!
time_t tstart[MAXNT];
double eltime[MAXNT];

int globalP=1;
int fType=0; 
int rStart=0; 
int per=1; 
int intBC=0; 
int timing=0;
int wGraph=0; 
int prCd=0; 
int wPetsc=0; 
int zScale=0;
int justpart=0; 
int parted=1;
int ICtyp;
extern int adaptFlag;
extern int strategy;
extern int phastaIO;
extern int old_format;
extern int WRITEASC;
extern int lstep;
extern int NSFTAG;
extern int numParts;
extern int numTotParts;
extern pProgress prog;
extern pMeshDataId incorp;
extern pMeshDataId phasta_solution;

// argument is attribute filename, model, parallel mesh 
void mdb2phasta(char* fname, char* mname,pGModel model,std::vector<pParMesh> pmesh,int nshgTot )
{

  pPList  bdry;

  pMesh mesh;
  pParMesh pMesh;

  FILE    *fprocs;
  FILE    *fstart;
  FILE    *attFile;
  FILE    *chkfp;
  int ipart;
  int Nodecount = 0;
  int Edgecount = 0;
  int Facecount = 0;

  globalInfo **info = new globalInfo* [numParts];
  for(ipart=0;ipart<numParts;ipart++)
      info[ipart] = new globalInfo;

// PMU already initialized (?)
// old   MS_init();
// old   MA_init();


  tstart[0] = time(NULL);
  tstart[1] = time(NULL);

  // default attribute filename

//old  strcpy(fname, "geom.atdb");

  // process command line arguements
  // procArgs overwrites geom.sms by the mesh directory name specified
  // here: 2stage in the 2nd stage
  //procArgs(argc, argv, fname, mname);

  if (PMU_rank()==0){
    cout << "\n\nGenerating input files for PHASTA\n\n";
  }

  // get the on-proc mesh(es): here only 1 part per proc
//  mesh = PM_mesh(pmesh, 0);
  // cout<<"\n PM_mesh(pmesh, 0)  success\n"; 

  eltime[1] = difftime(time(NULL), tstart[1]);

  /* P-change */
  if (PMU_rank()==0){
      cout << "\nPolynomial order: " << globalP << "\n";
  }


  // read the internal boundary conditions (if present)
#ifdef PARALLEL
// This part need to be taken care later for multiple parts per proc. Min Zhou
//   if (intBC == 1){
//       readBC(model, pmesh);
//       cout<<"\n readBC(model, pmesh) success\n";
//   }
#else
  if (intBC == 1){
      readBC(model, mesh);
  }
#endif

  pAttachDataCommu adc;
  // only do this if mesh not prevoiusly adapted
  if(!adaptFlag){

      // tag used to identify data that is to be migrated 
      // the data attached via  incorp="smsNum" is later on
      // retrieved in:
      // 1) setup.cc     : 2-step parNSpre, writing intermediate ncrptmp.%d  
      // 2) writeNCorp.c : writing the final ncvec
      incorp = MD_newMeshDataId("smsNum");
      
#ifdef PARALLEL
// from /net/vistmp1/kjansen/P-NSpre-dist/phNSpre/src/
      adc = PM_newAttachDataCommu(1,0,1);
      if(adc){ //Note this function is currently only called for Simmetrix lib
               // supports only one part per processor
          MD_setMeshCallback(incorp, CBmigrateOut, pm_sendAnInt, adc);
          MD_setMeshCallback(incorp, CBmigrateIn, pm_recvAnInt, adc);
          pMesh = pmesh[0]; 
          PM_setMigrId(pMesh, incorp);
      }

      if(PMU_rank()==0){
          cout<<"\n PMU_newAttachDataCommu success\n";
          cout<<"\n MD_setMeshCallback success\n";
          cout<<"\n PM_setMigrId success\n";
      }
  }//(!adaptFlag)
  
#endif 

  //  Here we find nshgTot (total number of DOFs) and attach the global number to the mesh
  //  for migration (depending on the globalP)
  if(globalP >3){
      if(PMU_rank()==0){
          cout<<"\nerror in mdb2phasta: global poly order too high at the moment\n";
      }
      exit(1);
  }

#ifdef PARALLEL
  // deprecated 1st stage
  if(justpart == 1){
      if(PMU_rank()==0){
         cout<<"\nerror in mdb2phasta:just partitioning deprecated, exiting \n";
      }
      exit(1);
  }

#else//SERIAL
  vertexModes= M_numVertices(mesh);
  if(globalP >1){
      edgeModes=(globalP-1)* M_numEdges(mesh);
      if(globalP >2){
          faceModes=((globalP-2)*(globalP-1)/2)*M_numFaces(mesh);
      }
  }
  nshgTot=vertexModes+edgeModes+faceModes;
#endif

  // only when adaption is NOT used
  if(!adaptFlag){
      nshgTot = 0;
#ifdef SIM
/*      for(ipart=0;ipart<numParts;ipart++){
          mesh = PM_mesh(pmesh[ipart],0);
      // only deal with a partitioned mesh
      // using -P option == parted
      // read an already partitioned mesh and the tags. 
      // These were printed out
      // in setup.cc when the -J flag was on - now read it under the -P flag
      
      // in the -P case old way of calculating nshgTot will fail as
      // there was never a total mesh.  Find the maximum number in the
      // list that was read in and use mpi to deduce nshgTot
      char nctmpstr[200];
      FILE *fp;
      void* tmp;
      int ireadthis;
      pVertex v; 
      sprintf(nctmpstr, "%s/ncrptmp.%d", mname,PMU_rank()+1);
      
      fp = fopen(nctmpstr, "r"); 
      
      VIter vIter = M_vertexIter(mesh);
      while( v = VIter_next(vIter) ) {
          
          fscanf(fp,"%d",&ireadthis);
      
          EN_attachDataInt((pEntity)v,incorp,ireadthis);
          
          if(ireadthis>nshgTot) nshgTot=ireadthis;
      }
      VIter_delete(vIter);
      fclose(fp);
      nshgTot=nshgTot+1;  // C numbering
      }
      int nshgg;
      MPI_Allreduce(&nshgTot, &nshgg, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      nshgTot=nshgg;
*/
#ifdef DEBUG     
//      fprintf(stderr,"read in a partitioned mesh:my rank is %d and I think nshgTot equals %d \n",
//              PMU_rank(),nshgTot );
#endif
#endif
      
  }//end of attaching global id for migration
  /* periodic boundary conditions (could be moved after partitioning, if mesh
   * matching info were present in the model) */
  tstart[4] = time(NULL);
  setPeriodic(model);
  eltime[4] = difftime(time(NULL), tstart[4]);
  
  //  functions relating to global multiprocessor information
  for(ipart=0;ipart<numParts;ipart++)
      initGlobalInfo(info[ipart]);   // Has to be called before partitioning
  
  tstart[5] = time(NULL);
  
#ifdef PARALLEL
  // lets be sure that everybody can get the information it needs
  // regarding nshgTot and incorp
  // rank root (here 0) broadcasts a message  to all processes
  // of the group, here ALL processes  
  MPI_Bcast(&nshgTot, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  pMesh = pmesh[0];
  for(ipart=0;ipart<numParts;ipart++) { //multiple parts per proc
#ifdef FMDB     
      pMesh = pmesh[ipart];
      mesh = PM_mesh(pmesh[ipart],0);
#endif 
#ifdef SIM
      mesh = PM_mesh(pMesh, ipart);   
#endif      
      info[ipart]->nshgTot = nshgTot;
#ifdef DEBUG
      if(PMU_rank()==0 && ipart==0)
          fprintf(stderr, "WARNING: nshgTot is not maintained in this code! \n");
#endif
//       fprintf(stderr,"my rank is %d and I think nshgTot equals %d \n",
//               PMU_rank(),nshgTot );

  
      // pMeshDataId
      info[ipart]->incorp = incorp;

      
      /* The next command splits up the mesh so if you want to calculate */
      /* global stuff, you must do so before this point */
      
#ifdef PARALLEL
      setup(pMesh, info[ipart], ipart);  
  
      // printf("[%d] memory usage after setup: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
  
  
  if (PMU_rank()==0&&ipart==numParts-1){
      cout << "\nsetup(pmesh, info) success  \n";
  }
  
  // if -J option has been specified (initial stage)
  // set to 1 in procArgs 
  if (justpart){
      if (PMU_rank()==0){
          cout<<"partitioning deprecated in this module-use Mitosis\n";
      }
      exit(0);
  }

#else // not PARALLEL
  setup(mesh, info);    // Only sets some things in info
#endif
        
  VIter vIter= M_vertexIter( mesh );
  pVertex v;
  while ( v = VIter_next( vIter ) ) {
      if( EN_isOwnerPart  (( pEntity )v , ipart))
          Nodecount++;
  }
  VIter_delete(vIter);

  pEdge e;
  EIter eIter = M_edgeIter( mesh );
  while ( e = EIter_next( eIter ) ) {
      if( EN_isOwnerPart  (( pEntity )e , ipart))
          Edgecount++;
  }
  EIter_delete(eIter);
  pFace f;
  FIter fIter = M_faceIter( mesh );
  while ( f = FIter_next( fIter ) ) {
      if( EN_isOwnerPart  (( pEntity )f , ipart))
          Facecount++;
  }    
  FIter_delete(fIter);  

} //loop over the parts on each proc

#ifdef PARALLEL
  MPI_Allreduce(&Nodecount, &totNumNodes, 1, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&Edgecount, &totNumEdges, 1, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&Facecount, &totNumFaces, 1, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD);  

  eltime[5] = difftime(time(NULL), tstart[5]);

  tstart[10] = time(NULL);

  initLocalInfo(pmesh, info);
  if(PMU_rank()==0){
      printf("\n[%2d] initLocalInfo(pmesh, info) success\n",PMU_rank());
      cout<<"Time Usage of iniLocalInfo is :"<<time(NULL)-tstart[10] <<"\n";
   }
#ifndef ibm
  //printf("[%2d] memory usage after setup: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
#else
  totNumNodes = Nodecount;
  totNumEdges = Edgecount;
  totNumFaces = Facecount;
  initLocalInfo(pMesh,info);
  cout<<"\n  initLocalInfo(mesh, info) success\n";
#endif
  
  
  eltime[10] = difftime(time(NULL), tstart[10]);

  /**********************************************************************/
  /* setup the global communication tasks                               */
  /**********************************************************************/
#ifdef PARALLEL
  tstart[6] = time(NULL);
  setupGlobalTasks(pmesh, info);
  if(PMU_rank()==0){
      printf("\n[%2d] setupGlobalTasks(pmesh, info) success\n",PMU_rank());
      cout<<"Time Usage of setupGlobalTasks is :"<<time(NULL)-tstart[6]<<"\n";
   }
  //printf("[%2d] memory usage after setupGlobalTasks: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
  eltime[6] = difftime(time(NULL), tstart[6]);
#endif

  /* This is very expensive in terms of memory, try to get rid of this */
  /* allocate boundary integral lists */
  /**********************************************************************/
  /* creating the block structure                                       */
  /**********************************************************************/
  tstart[7] = time(NULL);

//  #ifdef PARALLEL
//    // second stage operation: mesh should be correctly assigned
//    // after being read in as parted 
//    // re-assign if in 1-stage (!parted)  although partition already took place
//    if(parted == 0 && justpart == 0){
//        mesh=PM_mesh(pmesh,0);
//    }
//  #endif
  double timeofwriteEnsa = 0.0;
  
  pMesh = pmesh[0];
  for(ipart = 0; ipart<numParts; ipart++){ //loop over multiple parts on each proc

  mesh = PM_mesh(pMesh, ipart);
  bdry = PList_new();
  genblock(mesh, info[ipart], bdry);
  if(PMU_rank()==0 && ipart==numParts-1)
      cout<<"\n genblock(mesh, info, bdry)  success\n";

  eltime[7] = difftime(time(NULL), tstart[7]);

  if(isReorder)
      R_reordering(mesh, ipart);

/**********************************************************************/
  /* write the files for ENSA                                           */
  /**********************************************************************/
  tstart[8] = time(NULL);
#ifdef PARALLEL
  writeEnsaFiles(pMesh, info[ipart], bdry, ipart);
  
#else
  writeEnsaFiles(mesh, info, bdry);
#endif  

  timeofwriteEnsa += time(NULL)-tstart[8];
  if(PMU_rank()==0&&ipart==numParts-1){
      cout<<"\n["<<PMU_rank()<<"] writeEnsaFiles success\n" << endl;
      cout<<"Time Usage of writeEnsaFiles is: "<<timeofwriteEnsa<<"\n";
   }
  eltime[8] = difftime(time(NULL), tstart[8]);
  PList_delete(bdry);


   double tmptime = time(NULL);
#ifdef PARALLEL
  printInfo(pMesh, info[ipart], ipart);
#else
  printInfo(mesh, info);
#endif

  if(PMU_rank()==0 && ipart==numParts-1){
      cout<<"\n printInfo success\n";
  //    cout<<"Time Usage for printInfo is: "<<time(NULL)-tmptime<<"\n";
  }
  if (WRITEASC)
    printPeriodicBC(mesh);

  // need to rewrite the mesh file since some edge directions could
  // have changed due to switching for periodicity etc...
  // for facemodes (higher order)
  if (info[ipart]->edgeson) {
    printf("Rewriting the SMS file, some edge/face orientations "
           "could have changed while taking care of periodicity \n");
#ifdef SIM
    M_write(mesh, "geom2.sms", 0, prog);
#else
    M_write(mesh, "geom2.sms", 0);
#endif
  }

  freeGlobalInfo(info[ipart]);
  delete info[ipart];

  } //loop over multi-part on each proc

  delete [] info;

  // 3rd arg. (en_type) : 0 for vertex, 1 for edge etc.
  if(adaptFlag && strategy!=8){
      cleanAttachedData(mesh,phasta_solution,0);
  }

  MD_deleteMeshDataId(phasta_solution);
  if(!adaptFlag) {
    MD_deleteMeshDataId(incorp);
    AttachDataCommu_delete(adc);
  }

//#ifdef PARALLEL
//   PM_delete(pmesh);
// #else
//   M_delete(mesh);
// #endif

  if (PMU_rank() == 0)  {
    eltime[0] = difftime(time(NULL), tstart[0]);

    cout << "\nElapsed time: " << (double)eltime[0]/60.0 << " min\n";
    if (timing)
      printTime(eltime, MAXNT);

#ifdef PARALLEL
    cout << "\nSee 'info.out.<rank+1>' for mesh information\n\n";
#else
    cout << "\nSee 'info.out' for mesh information\n\n";
#endif

    /* write some files needed for ensa */

    fprocs=fopen("numpe.in", "w");
    fstart=fopen("numstart.dat", "w");

    fprintf(fprocs, "%d\n", PMU_size()*numParts);
    fprintf(fstart, "%d\n", lstep);

    fclose(fprocs);
    fclose(fstart);
  }



}


void printTime(double *eltime, int nt)
{
  int i;
  FILE *fp=stdout;

  /* convert from seconds to minutes */
  for (i=0; i < nt; i++)
    eltime[i] /= 60.0;

  fprintf(fp, "\n\nTiming statistics (min):\n");
  fprintf(fp, "------------------\n");
  fprintf(fp, "\nTotal Time:\t%f \n", eltime[0]);

  fprintf(fp, "\nActivity:       \tTime (min)\t%% total\n");
  fprintf(fp, "---------       \t----------\t-------\n");
  fprintf(fp, "load mesh       \t%f\t%f\n", eltime[1], eltime[1]/eltime[0]);
  fprintf(fp, "readBC          \t%f\t%f\n", eltime[2], eltime[2]/eltime[0]);
  fprintf(fp, "setPeriodic     \t%f\t%f\n", eltime[4], eltime[4]/eltime[0]);
  fprintf(fp, "setup           \t%f\t%f\n", eltime[5], eltime[5]/eltime[0]);
  fprintf(fp, "  initLocalInfo \t%f\t%f\n", eltime[10], eltime[10]/eltime[0]);
  fprintf(fp, "  partitionMesh \t%f\t%f\n", eltime[11], eltime[11]/eltime[0]);
  fprintf(fp, "  set dof info  \t%f\t%f\n\n", eltime[12], eltime[12]/eltime[0]);

  fprintf(fp, "setupGlobalTasks\t%f\t%f\n", eltime[6], eltime[6]/eltime[0]);
  fprintf(fp, "createBoundary  \t%f\t%f\n", eltime[7], eltime[7]/eltime[0]);
  fprintf(fp, "writeEnsaFiles  \t%f\t%f\n", eltime[8], eltime[8]/eltime[0]);
  fprintf(fp, "  getConnectivit\t%f\t%f\n", eltime[14], eltime[14]/eltime[0]);
  fprintf(fp, "  attach BC's   \t%f\t%f\n", eltime[15], eltime[15]/eltime[0]);
  fprintf(fp, "  writeFiles    \t%f\t%f\n\n", eltime[16], eltime[16]/eltime[0]);
  fprintf(fp, "\n\n");

  fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////////
// adapt.cc
//
// adption main module
//
// J.Mueller/O.Sahni 2004-2007
///////////////////////////////////////////////////////////////////////////////////
#include "mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <iostream>
#include <fstream>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "SimAdvMeshing.h"
#include "MeshSimAdapt.h"
#ifdef FMDB
#include "MeshAdapt.h"
#include "phPartitionCallbacks.h"
#endif

#if ( defined SIM_PARASOLID )
  #include "SimParasolidKrnl.h"
#elif ( defined DISCRETE )
  #include "SimDiscrete.h"
#endif

#include "phParAdapt.h"
#include "phReadWrite.h"
#include "attachData.h"
#include "MeshSimInternal.h"

extern "C" int readLicenseFile(char*);
extern "C" int procSize();
// avoid elements at (bdry.) with no interior nodes/dofs
extern void fix4NodesOnSurface(pMesh mesh);

#ifdef __cplusplus
extern "C" {
#endif

int setSizeFieldFromAttachedData(pMesh mesh, 
                                 pMSAdapt simAdapter, 
                                 pMeshDataId attachedDataID);



int CBcounter = 0;
int delDblArraycounter = 0;
int delDblcounter = 0;

pMeshDataId phasta_solution;
pMeshDataId errorIndicatorID;
pMeshDataId ybarID;
pMeshDataId nodalGradientID;
pMeshDataId nodalHessianID;
pMeshDataId nodalSizeID;
pMeshDataId wallStressID;
pMeshDataId wallDistID;
pMeshDataId SeparatedID;
pMeshDataId ShockID;
pMeshDataId isOrgNodeID;
pMeshDataId nodalVorticityID;
pMeshDataId OrgSizeID;
pMeshDataId nodalDirectionID;
pMeshDataId QualitySizeID;
pMeshDataId modes;
pMeshDataId incorp;
pMeshDataId localGradientID; 
pMeshDataId localPatchVolID;
pMeshDataId localHessianID;
pMeshDataId numSurroundNodesID;
pMeshDataId locMaxInterpolErrorID;
pMeshDataId globMaxInterpolErrorID;
pMeshDataId interfaceMetricID; //CWS

pMeshDataId meshSizeID;

// to communicate solution data for nodes fixed on partition bdry.
pMeshDataId partBdryFix4SolutionID;

// currently used in isotropic sizefield case
pMeshDataId oldMeshSizeID;

extern char outputFormat[];
extern int timeStepNumber;
extern int isSizeLimit;
extern int numSplit;
extern int nYbarVars;

extern int preLBforAdaptivity;
extern double masterProcWgt;
extern int DisplacementMigration;
extern int isBLAdapt;
extern int isThickAdapt;
extern int dwalMigration;
int BCInflowFace = 0;
int BCInflowFaceTag = 0;

// for the switch of errors 
// either in restart.n.m or error.n.m
int errorName;
int numVars;
extern int lstep;
extern pProgress prog;
time_t wtimePoints[32];
clock_t start_t, diff;

rusage cputimePoints[2];

int 
adapt(  // parallel mesh
            pParMesh pmesh,
            //serial mesh
            pMesh mesh,
            // model
            pGModel model,
            //time step
	    int timeStepNumber,
	    // strategy is to specify 
	    // how to do adaptation (i.e., size-field or tag driven)
	    // 1-2 : size-field driven (for anisotropic) 
	    // 3-4 : tag driven (for isotropic)
	    // 5-6 : size-field driven (for isotropic)
	    // < 0 sets a manual mesh size-field
	    int strategy,
	    // factor is the constant appearing in the error expression
	    // for tag driven it is used to define threshold for refinement
	    // for size-field driven it is used to define the error tolerance
	    double factor,
	    // number of solution (field) variables (5 for incompressible)
	    int ndof,
	    // number of variables for error indicators (EI)
	    // (e.g., 5 for ybar & 10 for residual-based)
	    int nvar,
	    // the maximal mesh edge length allowed in mesh size 
	    double hmax,
	    // the minimal mesh edge length allowed in mesh size
	    double hmin,
	    // option is used to decide how to compute the error value
	    // provides different choices like analytic hessian, manual size-field etc.
	    // for isotropic (tag or size-field driven :
	    // use 3 EI for flow problem or use 1 EI for scalar problem
	    int option)
{  
  // number of solution variables (globally visible)
  numVars = ndof; 

  // polynomial order
  int poly = 1;   
  // number of nodes (as considering only linear basis)
  int nshg;       
  lstep = timeStepNumber;
  char model_file[256];
  char mesh_file[256];
//  char moutfile[256];

  char solution_file[256];
  char error_indicator_file[256];
  char ybar_file[256];
  char wall_stress_tag[256], wall_stress_file[256];


//  sprintf(moutfile, "refined_p.sms");

//  pMesh mesh;
/*
  if(PM_verify  (  pmesh ,0, prog) == 0){
      if (PMU_rank() == 0){
          printf("\nerror in adapt.cc: invalid parallel mesh read in\n");
      }
      SimPartitionedMesh_stop();
      exit(1);
  }
*/
  wtimePoints[2] = time(0);

  // assume each proc has EXACTLY ONE partition on it !!!
//  mesh = PM_mesh(pmesh, 0);
 

  wtimePoints[3] = time(0);

  // make sure there is only one part per proc
  if(PM_numParts(pmesh) > 1){
      if(PMU_rank() == 0){
          printf("\nerror in adapt:only one part per proc at the moment allowed \n");
      }
      SimPartitionedMesh_stop();
      exit(1);
  }

  // get the total number of verts/rgns overal procs/parts
  int nshgTot = getNSHGTOT( mesh);
  int nRgnsTot, nRgnsLoc = M_numRegions(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nRgnsLoc, &nRgnsTot, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
#ifdef DEBUG
  if(PMU_rank()==0) {
    printf("\nTotal number of nodes   (before adaptation): %d\n",nshgTot);
    printf("\nTotal number of regions (before adaptation): %d\n",nRgnsTot);
  }

  printf("\n-- Before Adaptation...\n");
  printf(" [%2d]: (BA) local # of vertices: %d\n", PMU_rank(),M_numVertices(mesh));
  printf(" [%2d]: (BA) local # of edges   : %d\n", PMU_rank(),M_numEdges(mesh));
  printf(" [%2d]: (BA) local # of faces   : %d\n", PMU_rank(),M_numFaces(mesh));
  printf(" [%2d]: (BA) local # of regions : %d\n", PMU_rank(),M_numRegions(mesh));
#endif

  wtimePoints[4] = time(0);

  nshg=M_numVertices(mesh);
  phasta_solution = MD_newMeshDataId("restart solution");
  errorIndicatorID = MD_newMeshDataId("error indicator");
  ybarID = MD_newMeshDataId("ybar");
  localPatchVolID = MD_newMeshDataId("local Patch Vol"); 
  localGradientID  = MD_newMeshDataId("local Gradient");
  localHessianID  = MD_newMeshDataId("local Hessian"); 
  nodalGradientID = MD_newMeshDataId("nodal Gradient");
  nodalHessianID  = MD_newMeshDataId("nodal Hessian");
  nodalSizeID = MD_newMeshDataId("nodal Size");
  wallStressID = MD_newMeshDataId("wall stress");
  wallDistID = MD_newMeshDataId("wall Dist");
  SeparatedID = MD_newMeshDataId("separated BL");
  ShockID = MD_newMeshDataId("Shock Detect");
  isOrgNodeID = MD_newMeshDataId("is Originating node");
  nodalVorticityID = MD_newMeshDataId("nodal vorticity");
  OrgSizeID = MD_newMeshDataId("Original Size");
  QualitySizeID = MD_newMeshDataId("Quality Size"); 
  nodalDirectionID = MD_newMeshDataId("nodal Direction"); 
  numSurroundNodesID = MD_newMeshDataId("number of surrounding nodes"); 
  locMaxInterpolErrorID = MD_newMeshDataId("local interpol error");
  globMaxInterpolErrorID = MD_newMeshDataId("global interpol error"); 
  interfaceMetricID = MD_newMeshDataId("interface metric");// CWS

  modes = MD_newMeshDataId("number of modes");// required for higher order

  pAttachDataCommu adc,adcSN;
///*
  //if(strategy!=7 && strategy!=8) { //CWS
  if(strategy!=7 ) {

    // only when required to be migrated over procs -- required for solution callback!
    int ndwal = 0;
    if(dwalMigration) {
      ndwal = 1;
    }
    adc = AttachDataCommu_new(sizeof(double)/sizeof(int),0,numVars+ndwal);
/*    
    if(adc){
        MD_setMeshCallback(phasta_solution , CBmigrateOut, pm_sendDblArray, adc );
        MD_setMeshCallback(phasta_solution , CBmigrateIn,  pm_recvDblArray, adc );
        MD_setMeshCallback(phasta_solution , CBdelete, delDblArray, NULL );
        PM_setMigrId( pmesh,phasta_solution  );
    }
*/    
    // cout<<"\n["<<PMU_rank()<<"]: PM_setMigrId(solution) success ...\n";
//*/
    incorp = MD_newMeshDataId("smsNum");
    // only when required to be migrated over procs
    adcSN = AttachDataCommu_new(1,0,1);
//     MD_setMeshCallback(incorp, CBmigrateOut, pm_sendAnInt, adcSN);
//     MD_setMeshCallback(incorp, CBmigrateIn, pm_recvAnInt, adcSN);
//     PM_setMigrId(pmesh, incorp);
    // cout<<"\n["<<PMU_rank()<<"]: PM_setMigrId(incorp) success ...\n";
  }

#ifdef SIM
  MSA_CallbackFunc*  phastaTransferFn = phastaTransferTopSIM;
#endif
#ifdef FMDB
  CBFunction phastaTransferFn = phastaTransferTopSCOREC;
#endif

  wtimePoints[5] = time(0);

  // parallel adapt object
  pMSAdapt simAdapter;

  SMeshChanges* mco;
  int ndofTot;
  // if strategy<0 (command line argument as of now) 
  // then sets a manual size field  
  switch(strategy) {
  case 1 :
  case 2 : { //anisotropic adaptation
    if(PMU_rank()==0)
      if(option!=9)
	printf("\nStrategy chosen for ANISOTROPIC adaptation : size-field driven\n");
      else
	printf("\nStrategy chosen for ISOTROPIC adaptation  (based on local hmin of mesh metric from hessians) : size-field driven\n");
    
    sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);
    if(isThickAdapt){
       sprintf(wall_stress_file,"restart.%i.%i",lstep,PMU_rank()+1);
       sprintf(wall_stress_tag,"wall shear stresses");
    }
    
    char error_tag[256], ybar_tag[256];
    if(option==1 || option==9) {
      if(strategy == 1) {
      	if(PMU_rank()==0)
	         printf("\nUsing ybar to compute hessians...\n\n");
	      sprintf(ybar_tag,"ybar");
         if(PMU_size()==1) {         
      	   sprintf(ybar_file,"ybar.%i.%i",lstep,PMU_rank()+1);
         } else {
      	   sprintf(ybar_file,"restart.%i.%i",lstep,PMU_rank()+1);
         }
      }
      else if (strategy == 2) {
	      if(PMU_rank()==0)
      	  printf("\nUsing numerical/computed hessians (i.e, from phasta)...\n\n");
	        sprintf(error_tag,"hessains");
           sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
      }
      
      printf("\n Reading files:\n");
      printf(" ...%s (for \"solution\")\n",solution_file);
      printf(" ...%s (for \"%s\")\n\n",ybar_file,ybar_tag);
//      if(isThickAdapt)
//        printf(" ...%s (for \"wall shear stress\")\n",wall_stress_file);
    }
    // strategy 1, option 2:  stationary  temperature hessians computed
    // from solution
    if(option==2) {
       sprintf(error_tag,"solution");// for stationary  temperature 
       sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
       printf("\nUsing TEMPERATURE field to compute hessians...\n\n");
    }

    if(option == 3) { //Isotropic Adaptation but passed as Anisotropic size field
       if(strategy == 1) {
         if(PMU_rank()==0)
          printf("\nUsing Error fields for isotropic adaption passed as anisotropic...\n\n"); 
       sprintf(error_tag,"errors");
	    sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
       }
      printf("\n Reading files:\n");
      printf(" ...%s (for \"solution\")\n",solution_file);
      printf(" ...%s (for \"%s\")\n\n",error_indicator_file,error_tag);
    }

    if(option == 4) { //Hessians for direction and errors for size field
       if(strategy == 1) {
         if(PMU_rank()==0)
           printf("\nUsing Error fields for size field and hessians for directions of anisotropic...\n\n"); 
       sprintf(error_tag,"errors");
       sprintf(ybar_tag,"ybar");

       if(PMU_size()==1) {         
   	    sprintf(error_indicator_file,"errors.%i.%i",lstep,PMU_rank()+1);
          sprintf(ybar_file,"ybar.%i.%i",lstep,PMU_rank()+1);
       } else {
   	    sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
          sprintf(ybar_file,"restart.%i.%i",lstep,PMU_rank()+1);
       }
       }
      printf("\n Reading files:\n");
      printf(" ...%s (for \"solution\")\n",solution_file);
      printf(" ...%s (for \"%s\")\n\n",error_indicator_file,error_tag);
      printf(" ...%s (for \"%s\")\n\n",ybar_file,ybar_tag);
    }

    wtimePoints[6] = time(0);

    // 2nd arg. (1/0 : size-field/tag  driven; 3rd:1/0 predict load balanc/not )

#ifdef DEBUG    
/*    
    if(PMU_size()==1) {
      pMesh meshMerge;
      meshMerge = M_createFromParMesh(pmesh, 3, prog);
      M_write(meshMerge, "mesh_new.sms", 0, prog);
      M_release(meshMerge);
    }   
*/    
#endif

//    if(PMU_size()==1) {
//      simAdapter = MSA_new(mesh,1);
//    }
//    else {
      simAdapter = MSA_new(pmesh,1);
//   }

#ifdef FMDB
    if(nvar){
        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(nvar);
        if(isBLAdapt)
           simAdapter->setPreLBflag(0);
        else 
	   simAdapter->setPreLBflag(1);
    }
#endif        

#ifdef SIM
  MSA_setAdaptBL(simAdapter, isBLAdapt);
  int localAdapt=1; // 0 if you want the whole mesh to adapt, 1 if you want only where you have set the size field
  MSA_setLocal(simAdapter, localAdapt);	
  if(isBLAdapt==1) {
     MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);
     MSA_setBLMinLayerAspectRatio(simAdapter, 1.0);
  }
// MSA_setBLSnapping(simAdapter, 0);
  
  if(PMU_size()>1 && isBLAdapt==1) {
//     MSA_setCoarsenMode(simAdapter, 0);
//     MSA_setBLSnapping(simAdapter, 0);
     MSA_setBLMinLayerAspectRatio(simAdapter, 0.0);
     MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
  }

  MSA_setCoarsenMode(simAdapter, 0);
  MSA_setMaxIterations(simAdapter, numSplit);
  MSA_setBoundaryMeshModification(simAdapter, isBLAdapt);

#endif

    // invoke solution transfer
//   MSA_setCallback(simAdapter, phastaTransferFn, 9, NULL);// 8 for regions, 1 for vertices
    wtimePoints[7] = time(0);

    if(option>0 && strategy==1) {

      wtimePoints[8] = time(0);

      // attaching the solution to the original mesh
      double *sol;
      readArrayFromFile(solution_file,"solution",sol);
      attachArray(sol,mesh,phasta_solution,ndof,poly);
      delete [] sol;

      // read ybar (and compute/use hessians of ybar) 
      double *error_indicator;
      double *ybar;
      readArrayFromFile(ybar_file,ybar_tag,ybar);
      attachArray(ybar,mesh,ybarID,nYbarVars,poly);
      delete [] ybar;
      if(option == 4) {
        readArrayFromFile(error_indicator_file,error_tag,error_indicator);
        attachArray(error_indicator,mesh,errorIndicatorID,nvar,poly);
        delete [] error_indicator;
      }

      M_writeVTKFile(mesh,"nodalEI", ybarID, 6);
      if(isThickAdapt){
         // read wall shear stress  
           WallStress(mesh);

         for(int i=0; i<10;i++) 
            SmoothWallField(mesh, wallStressID, 1);

//         M_writeVTKFile(mesh, "SmoothwallStress", wallStressID, 1);

          Vorticity(pmesh, mesh);
          TotalBLHeight(mesh);
  
//         DetectShock(pmesh, mesh);
         for(int i=0;i<10;i++)
            SmoothWallField(mesh, wallDistID, 3);
         
         DistOffWall(mesh, simAdapter);
         CalcYPlus(mesh);

//         delete [] wall_stress;
      }
      
      
      wtimePoints[9] = time(0);

      wtimePoints[10] = time(0);

      // calculating hessians for ybar field
      // first reconstruct gradients and then the hessians 
      // also deals with boundary issues &
      // applies smoothing procedure for hessians
      // (simple average : arithmetic mean)
      if (option != 3)
         hessiansFromSolution(pmesh,mesh,lstep);     

#ifdef DEBUG      
//      M_writeVTKFile(mesh, "nodalGrad", nodalGradientID, 3);
//      M_writeVTKFile(mesh,"nodalHess", nodalHessianID, 6); 
      M_writeVTKFile(mesh,"nodalEI2", ybarID, 6);
#endif      
      wtimePoints[11] = time(0);
    }
    else if (strategy == 2) { // cannont use analytic hessian in this case
//        // option : choice of variable (whose hessians are taken)
//        // use the hessians computed from phasta
//        double *hessiansFromRestart;
//        readArrayFromFile(error_indicator_file,error_tag,hessiansFromRestart);
//        double *hessians = new double[nshg*6];
//        V_getHessians(hessiansFromRestart,mesh,nvar,option,hessians);
//        delete [] hessiansFromRestart;
//        attachArray(hessians,mesh,nodalHessianID,6,poly);
//        delete [] hessians;

//        // once option is used in V_getHessians 
//        // it is changed so that we don't use analytic hessian (in sizefield.cc)
//        option = 1;

//        // simple average over a patch surrounding the vertex 
//        SmoothHessians(mesh);      
    }

    wtimePoints[12] = time(0);
    
    // compute mesh size-field using hessian strategy (anisotropic adaptation)
    // and set it at each vertex
    if(option == 4) {
      transformToScalarErrorVal(mesh,nvar,option);
      setSizeFieldUsingCombined(pmesh,mesh,simAdapter,factor,hmax,hmin,option);
    }
    if(option != 3 && option != 4) {
      setSizeFieldUsingHessians(pmesh,mesh,simAdapter,factor,hmax,hmin,option);
    }  
    if(option == 3) {  
      transformToScalarErrorVal(mesh,nvar,option);
      setIsotropicSizeField(pmesh,mesh,simAdapter,factor,hmax,hmin,option);
    }
    wtimePoints[13] = time(0);
  }// case 2:  anisotropic adaptation
  break;
  case 3 :
  case 4 : { // isotropic adaptation (tag driven)
    printf("\nStrategy chosen for ISOTROPIC adaptation : tag driven\n");
    sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);

    char error_tag[256];
    if(strategy == 3) {
      sprintf(error_tag,"error");
      if(PMU_size()==1) {
         sprintf(error_indicator_file,"error.%i.%i",lstep,PMU_rank()+1);
      } else { 
         sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
      }
    }
    else if (strategy == 4) {
//      sprintf(error_tag,"time derivative of solution");
      sprintf(error_tag,"error");
      sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);

    }

    printf("\n Reading files:\n");
    printf(" ...%s (for \"solution\")\n",solution_file);
    printf(" ...%s (for \"%s\")\n\n",error_indicator_file,error_tag);
    
    // (2nd 0tag/1 sizefield;   3rd 1 : predictive load balancing, 0 not)

    if(PMU_size()==1) {
        simAdapter = MSA_new(pmesh,0);
    }
    else {
      simAdapter = MSA_new(pmesh,0);
   }

#ifdef FMDB
    if(nvar){
        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(nvar);
    }
#endif   
//   MSA_setCallback(simAdapter, phastaTransferFn, 9, NULL);// 8 for regions, 1 for vertices
#ifdef SIM
    MSA_setAdaptBL(simAdapter, isBLAdapt);
    MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);
// MSA_setBLSnapping(simAdapter, 0);
#endif

    // attaching the solution to the original coarse mesh
    double *sol;
    readArrayFromFile(solution_file,"solution",sol);
    attachArray(sol,mesh,phasta_solution,ndof,poly);
    delete [] sol;
    
    // read residual based error indicators (10 quantities)
    double *error_indicator;
    readArrayFromFile(error_indicator_file,error_tag,error_indicator);
    attachArray(error_indicator,mesh,errorIndicatorID,nvar,poly);
    delete [] error_indicator;

    // this has to be done for all parts, so:
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i=0; i<20; i++) {
       SmoothErrorIndicators(mesh,option);
    }
    // tag the entities to be refinement (for isotropic refinement)
    // factor is used to evaluate the threshold for refinement
    // as of now do not use hmin and hmax
    // tagEntitiesForRefinement uses applyMarkingStrategy
#ifdef DEBUG    
//    M_writeVTKFile(mesh,"nodalEI", errorIndicatorID, nvar);
#endif    
    tagEntitiesForRefinement(mesh,simAdapter,factor,hmax,hmin,option);
  }
  break;
  case 5 :
  case 6 : { //isotropic adaptation (size-field driven)
    printf("\nStrategy chosen for ISOTROPIC adaptation : size-field driven\n");

    sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);

    char error_tag[256];
    // read residual based error indicators (10 quantities)
    // only for option=1 which sets size field based on error indicators
    // option=2 sets based on interface proximity
    // (for 2, load in solution for error, although it is not used)
    if(strategy == 5) {
      if(option==1 || option==10 || option==11) { 
        sprintf(error_tag,"errors");
        if(PMU_size()==1) {
           sprintf(error_indicator_file,"errors.%i.%i",lstep,PMU_rank()+1);
        } else {
           sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
        }
      } else if (option==2) {
        sprintf(error_tag,"solution");
        sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
      } else if (option==5) {
        sprintf(error_tag,"error"); 
        sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
      }
        
    }
    else if (strategy == 6) {
      sprintf(error_tag,"errors");
      sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);
    }

    printf("\n Reading files:\n");
    printf(" ...%s (for \"solution\")\n",solution_file);
    printf(" ...%s (for \"%s\")\n\n",error_indicator_file,error_tag);
    
    wtimePoints[6] = time(0);

    // args: pParMesh, 0/1 (tag/sizefield), 0/1 non-predictive/predictive load balancing
    if(PMU_size()==1) {
        simAdapter = MSA_new(pmesh,1);
    }
    else {
      simAdapter = MSA_new(pmesh,1);
   }

#ifdef FMDB
    if(nvar){
        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(nvar);
    }
#endif   
//    MSA_setCallback(simAdapter, phastaTransferFn, 9, NULL);// 8 for regions, 1 for vertices
    MSA_setAdaptBL(simAdapter, isBLAdapt);

#ifdef SIM
//     MSA_setCoarsenMode(simAdapter, 0);
     MSA_setMaxIterations(simAdapter, numSplit);
     if(isBLAdapt) {
        MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
        MSA_setBLMinLayerAspectRatio(simAdapter, 0.0);
     }
     MSA_setBoundaryMeshModification(simAdapter, isBLAdapt);

#endif
    wtimePoints[7] = time(0);

    wtimePoints[8] = time(0);

    // attaching the solution to the original coarse mesh
    double *sol;
    readArrayFromFile(solution_file,"solution",sol);
//    attachArray(sol,mesh,phasta_solution,ndof,poly);
//    delete [] sol;

    int ndwal = 0;
    double* dwal;
    if(dwalMigration) {
      readArrayFromFile(solution_file,"dwal",dwal);
      ndwal = 1;
    }

    ndofTot = ndof+ndwal;
    double *SolTot = new double[ndofTot*nshg];

    for(int inode=0;inode<nshg;inode++) {
        for(int idof=0;idof<ndof;idof++)
            SolTot[inode*ndofTot+idof] = sol[inode*ndof+idof];
        for(int idof=0;idof<ndwal;idof++) {
            SolTot[inode*ndofTot+ndof+idof] = dwal[inode*ndwal+idof];
        }
            
    }

    attachArray(SolTot,mesh,phasta_solution,ndofTot,poly);
 
    if(dwalMigration) delete[] dwal;
    delete [] sol;
    delete [] SolTot;

    
    // read residual based error indicators (10 quantities)
    double *error_indicator;
    readArrayFromFile(error_indicator_file,error_tag,error_indicator);

    if (option==10) {
      double *ybar_indicator;
      readArrayFromFile(error_indicator_file,"ybar",ybar_indicator);
      int nshg = M_numVertices(mesh);
      for(int inode=0;inode<nshg;inode++) {
         //3 because it is the 4th variable to be overwritten with the last variable 12 of the ybar field. WARNING: HARD CODED time avg eddy viscosity overwriting the first diffusive flux:w
         error_indicator[inode*10+3] = ybar_indicator[inode*13+12]; 
       }
       delete [] ybar_indicator;
    }
    else if (option==11) {
      double *dwal_indicator;
      readArrayFromFile(error_indicator_file,"dwal",dwal_indicator);
      double *ybar_indicator;
      readArrayFromFile(error_indicator_file,"ybar",ybar_indicator);
      int nshg = M_numVertices(mesh);
      for(int inode=0;inode<nshg;inode++) {
         //WARNING: HARD CODED dwal overwriting the first diffusive flux:w
         error_indicator[inode*10+4] = ybar_indicator[inode*13+12]; 
         error_indicator[inode*10+3] = dwal_indicator[inode]; 
       }
       delete [] dwal_indicator;
    }
       

    attachArray(error_indicator,mesh,errorIndicatorID,nvar,poly);
    delete [] error_indicator;

    wtimePoints[9] = time(0);

    wtimePoints[10] = time(0);

    // currently we are smoothing the residual based part 
    // of the set of all EIs given through phasta (which might be already smoothed)
    // the quantity that is smoothed (arithmetic mean) and attached:
    // res_u^2 + res_v^2 + res_w^2
    // (the first three components from the set of ten)
    // the ID used to attach the smoothed EIs is "errorIndicatorID"
    // data previously attached under that ID is deleted
//    SmoothErrorIndicators(mesh,option);
    
//    writeSmoothEIs(mesh);
    
    // transform To Scalar Error Val --> required by setIsotropicSizeField
    // transforms the 10 attached EI into ONE single value
    // according to individual needs
    // error will be attached via errorIndicatorID
//KEJ    for (int i=0; i<20; i++) {
//KEJ       SmoothErrorIndicators(mesh,option);
//KEJ    }
//KEJ tried the above but it seems to be broken as it needs an option set which is not set in the location elsewhere in this routine I copied from....suspect routines have changed since last used
    
    transformToScalarErrorVal(mesh,nvar,option);

    wtimePoints[11] = time(0);

    wtimePoints[12] = time(0);

    // compute mesh size-field using residues (isotropic adaptation)
    // and set it at each vertex 
    // uses a single scalar nodal value attached via errorIndicatorID
    setIsotropicSizeField(pmesh,mesh,simAdapter,factor,hmax,hmin,option); 

    wtimePoints[13] = time(0);

    //M_writeVTKFile(mesh, "nodalSolution", phasta_solution, 6);
#ifdef DEBUG
//      M_writeVTKFile(mesh, "nodalSolution", phasta_solution, 6);
#endif    
  }   
  break;
  case 7: { //octree splitting only for generating BIG meshes

    // args: pParMesh, 0/1 (tag/sizefield), 0/1 non-predictive/predictive load balancing
/*
    int ndisp = 0;
    if(DisplacementMigration)
        ndisp = 3;
    sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);
      // attaching the solution to the original mesh
    double *sol, *disp;
    readArrayFromFile(solution_file,"solution",sol);
    if(DisplacementMigration)
        readArrayFromFile(solution_file, "displacement", disp);
    double *SolTot = new double[(ndof+ndisp)*nshg];

    for(int inode=0;inode<nshg;inode++) {
        for(int idof=0;idof<ndof;idof++)
            SolTot[inode*(ndisp+ndof)+idof] = sol[inode*ndof+idof];
        for(int idof=0;idof<ndisp;idof++)
            SolTot[inode*(ndisp+ndof)+ndof+idof] = disp[inode*ndisp+idof];
    }
    ndofTot = ndof+ndisp;

    attachArray(SolTot,mesh,phasta_solution,ndofTot,poly);

    delete [] sol;
    if(DisplacementMigration)
       delete [] disp;
    delete [] SolTot;
*/
    if(PMU_size()==1) {
        simAdapter = MSA_new(mesh,0);
    }
    else {
      simAdapter = MSA_new(pmesh,0);
   }

#ifdef SIM
   MSA_setAdaptBL(simAdapter, isBLAdapt);
//  MSA_setPrebalance(simAdapter, 0);

 if(isBLAdapt == 1) {  
   if(PMU_size()>1) {

//     MSA_setCoarsenMode(simAdapter, 0);
     MSA_setPrebalance(simAdapter, 0);
     MSA_setBLSnapping(simAdapter, 0);
     MSA_setBLMinLayerAspectRatio(simAdapter, 0.0);
     MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);
  }
 }  
#endif

#ifdef FMDB
    if(ndofTot){
        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(ndof);
    }
#endif
//      MSA_setCallback(simAdapter, phastaTransferFn, 9, NULL);
//    MSA_setCallback(simAdapter, NULL,9,NULL);

    pEdge edge;
    EIter eIter = M_edgeIter(mesh);
/*
    pVertex vert;
    VIter vIter = M_vertexIter(mesh);
    while(vert = VIter_next(vIter)) {
       if(EN_isBLEntity(vert)) cout << "is BL Entity" << endl;
    }
*/
    while(edge = EIter_next(eIter)) { 
      MSA_setRefineLevel(simAdapter,(pEntity)edge,1);
    }
    EIter_delete(eIter);
//    PM_write(pmesh, "mesh_size.sms", prog);
  }
  break;
  case 8: {
    if ( PMU_rank() == 0 ) {
       printf("\n CASE 8\n"); 
       printf("\n Strategy chosen for reading size-field from file to drive adaptation\n",PMU_rank());
    }

    sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);

    char error_tag[256];
    sprintf(error_tag,"error");
    sprintf(error_indicator_file,"restart.%i.%i",lstep,PMU_rank()+1);

    char interfaceMetric_file[256];
    sprintf(interfaceMetric_file,"gradphi.%i.%i",lstep,PMU_rank()+1); 

    if ( PMU_rank() == 0 ) {
      printf("\n Reading files:\n");
      printf(" ...%s (for \"solution\")\n",solution_file);
      printf(" ...%s (for \"%s\")\n\n",error_indicator_file,error_tag);
      printf(" ...%s (for \"metric\")\n\n",interfaceMetric_file); 
    }
    
    // args: pParMesh, 0/1 (tag/sizefield), 0/1 non-predictive/predictive load balancing
    if(PMU_size()==1) {
        simAdapter = MSA_new(pmesh,1);
    }
    else {
      simAdapter = MSA_new(pmesh,1);
   }

#ifdef FMDB
    if(nvar){

        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(nvar);
    }
#endif   
//    MSA_setCallback(simAdapter, NULL,9,NULL);// 8 for regions, 1 for vertices

    // attaching the solution to the original coarse mesh
    double *sol;
    readArrayFromFile(solution_file,"solution",sol);
    attachArray(sol,mesh,phasta_solution,ndof,poly);
    delete [] sol;

    // attaching the interface metric to the original coarse mesh
    double *interfaceMetric;
    //HACK we are using the write_gradphi function from phSolver 
    readArrayFromFile(interfaceMetric_file,"gradphi",interfaceMetric); 
    attachArray(interfaceMetric,mesh,interfaceMetricID,9,poly); 
    delete [] interfaceMetric;

    //setIsotropicSizeField(pmesh,mesh,simAdapter,factor,hmax,hmin,option); 

    setSizeFieldFromAttachedData(mesh,simAdapter,interfaceMetricID); 

  }
  break;
  default : {
    if(strategy<0) {
      printf("\nStrategy chosen for adaptation : manual size-field\n");

      // 2nd arg. (1 : size-field driven adaptation)
           // attaching the solution to the original mesh
      double *sol;
      char error_tag[256];
      sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);
//      readArrayFromFile(solution_file,"solution",sol);
//      attachArray(sol,mesh,phasta_solution,ndof,poly);
//      delete [] sol;
	   sprintf(error_tag,"ybar");
	   sprintf(error_indicator_file,"ybar.%i.%i",lstep,PMU_rank()+1);
      double *error_indicator;
      if(option==-1) {
         readArrayFromFile(error_indicator_file,error_tag,error_indicator);
         attachArray(error_indicator,mesh,errorIndicatorID,nvar,poly);
         delete [] error_indicator;
      }

    if(PMU_size()==1) {
        simAdapter = MSA_new(pmesh,1);
    }
    else {
      simAdapter = MSA_new(pmesh,1);
   }

#ifdef FMDB
    if(nvar){
        simAdapter->setSolutionID(phasta_solution);
        simAdapter->setnVar(nvar);
    }
#endif 
      // set mesh size-field manually (anisotropic adaptation) at each vertex
      // strategy provides different choices
      setManualSizeField(pmesh, mesh,simAdapter,strategy, option);
//      MSA_setCallback(simAdapter, NULL,9,NULL);
//      NormalizedEdgeLength(mesh);

#ifdef SIM
    printf("isBLAdapt=%i\n\n",isBLAdapt);
    MSA_setAdaptBL(simAdapter, isBLAdapt); 
   if(isBLAdapt == 1) { 
     MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);

      if(PMU_size()>1) {

//     MSA_setCoarsenMode(simAdapter, 0);
//     MSA_setBLSnapping(simAdapter, 0);
     MSA_setBLMinLayerAspectRatio(simAdapter, 0.0);
     MSA_setExposedBLBehavior(simAdapter, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);

    }
   }
#endif
    }
    else {
      printf("\nSpecify a correct (adaptation) strategy (adapt.cc)\n");
      exit(-1);
    }
    } 
    break;
  }

  if(preLBforAdaptivity)
    mesh = PM_mesh(pmesh,0);

  if(strategy==1) { 
    cleanAttachedData(mesh,ybarID,0);
    if(option==4) 
      cleanAttachedData(mesh,errorIndicatorID,0);
  }
  if(strategy==3 || strategy==4)
    cleanAttachedData(mesh,errorIndicatorID,0);
  else if (strategy==5 || strategy==6)
    cleanAttachedData(mesh,errorIndicatorID,0,0);

//   if(masterProcWgt>0.) {
//     int tnp = PM_totNParts(pmesh);
//     double *partwts = new double[tnp];

//     if(tnp>1 && PMU_rank()==0)
//       printf("\n Master processor weight : %f\n\n",masterProcWgt);
    
//     partwts[0] = masterProcWgt;
//     for(int iPart=1; iPart<tnp; iPart++)
//       partwts[iPart] = 1.;
    
//     if(tnp>1)
//       MSA_setPartWts(simAdapter,partwts); 
    
//     delete [] partwts;
//   }

  MPI_Barrier(MPI_COMM_WORLD);

  wtimePoints[14] = time(0);
  
  getrusage(RUSAGE_SELF,&cputimePoints[0]);

  start_t = clock();

  // adaptation
#ifdef SIM
  MSA_adapt(simAdapter, prog);
#else
  MSA_adapt(simAdapter);
#endif
  MSA_delete(simAdapter);

  getrusage(RUSAGE_SELF,&cputimePoints[1]);

  wtimePoints[15] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);
  diff = clock()-start_t;
  cout << "Time for adaptation: " <<diff/CLOCKS_PER_SEC <<"\n";

  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[16] = time(0);

#ifdef SIM
 // takes case of bad brdy. elements (elements with no interior nodes)
//  fix4NodesOnSurface(mesh);
#endif
  wtimePoints[17] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[18] = time(0);

#ifdef SIM
  // fix4NodesOnSurface introduces new vertices for which 
  // a solution is needed
  //if(strategy>0 && (strategy!=7 && strategy!=8)) { // CWS
  if(strategy>0 && (strategy!=7 )) {
    partBdryFix4SolutionID = MD_newMeshDataId("part. bdry. fix 4 solution");
//    fix4SolutionTransfer(mesh);
    MPI_Barrier(MPI_COMM_WORLD);

    // do not know whether this would happen or not
    // (i.e., new nodes created by fix4NodesOnSurface at partition bdry.)
    // in any case does not hurt
//    commuFix4SolutionTransfer(pmesh,mesh);
    MPI_Barrier(MPI_COMM_WORLD);

    MD_deleteMeshDataId(partBdryFix4SolutionID);
  }
#endif

  wtimePoints[19] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
  printf("-- Adaptation Done...\n");
  printf(" [%2d]: (AA) local # of vertices: %d\n", PMU_rank(),M_numVertices(mesh));
  printf(" [%2d]: (AA) local # of edges   : %d\n", PMU_rank(),M_numEdges(mesh));
  printf(" [%2d]: (AA) local # of faces   : %d\n", PMU_rank(),M_numFaces(mesh));
  printf(" [%2d]: (AA) local # of regions : %d\n", PMU_rank(),M_numRegions(mesh));
  printf("\n [%2d]: Number of sol. trans. callbacks: %d\n\n",PMU_rank(),CBcounter);
  printf("\n [%2d]: Number of del. array  callbacks: %d",PMU_rank(),delDblArraycounter);
  printf("\n [%2d]: Number of del. single callbacks: %d\n\n",PMU_rank(),delDblcounter);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[20] = time(0);
  

#ifdef FMDB
/*  phPartitionCallbacks phCB(phasta_solution,ndofTot);
  M_loadbalance(mesh,phCB);
  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[31] = time(0);
  if(!PMU_rank())
      std::cout<<"Time for M_loadbalance is t = "<<wtimePoints[31]-wtimePoints[20]<<"(s)\n";
*/
#endif

  if(masterProcWgt>0.) {
    int tnp = PM_totalNumParts(pmesh);
    if(tnp>1) {

      printf("\nPartitioning mesh equally after adaptation \n\n");

      pPartitionOpts popts = PartitionOpts_new();
      PartitionOpts_setTotalNumParts(popts,tnp);
      PartitionOpts_setPartWtEqual(popts);
      
      PM_partition(pmesh,popts,prog);
      
      PartitionOpts_delete(popts);

      printf("\nMesh partitioned equally after adaptation \n\n");
      
      mesh = PM_mesh(pmesh,0);

      printf("-- Equally partitioned mesh after adaptation...\n");
      printf(" [%2d]: (PAA) local # of vertices: %d\n", PMU_rank(),M_numVertices(mesh));
      printf(" [%2d]: (PAA) local # of edges   : %d\n", PMU_rank(),M_numEdges(mesh));
      printf(" [%2d]: (PAA) local # of faces   : %d\n", PMU_rank(),M_numFaces(mesh));
      printf(" [%2d]: (PAA) local # of regions : %d\n", PMU_rank(),M_numRegions(mesh));
    }
  }

  wtimePoints[21] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[22] = time(0);


// Must write out dependicies(smd) prior to resurce (sms)
#ifdef SIM
//  GM_write(model, "thisModel.smd",0, prog);
  if(PMU_size()==1) {
     M_write(mesh, "mesh_out.sms", 0, prog);
     M_writeVTKFile(mesh, "mesh_out", phasta_solution, ndof);
//  } else {
//     PM_write(pmesh, "mesh_out.sms", prog);
  }
#endif
#ifdef DEBUG
//   M_writeVTKFile(mesh, "mesh_out", 0, 0);
//   M_writeVTKFile(mesh, "mesh_out", phasta_solution, ndof);
#endif   
/*   
  pMesh meshMerge;
  meshMerge = M_createFromParMesh(pmesh, 3, prog);
  M_write(meshMerge, "mesh_outMerge.sms", 0, prog);
  M_writeVTKFile(meshMerge,"mesh_outMerge", phasta_solution, ndof);
  M_release(meshMerge);
*/
  // write parallel mesh
  if(PMU_rank()==0)
    printf("\n Writing out the parallel mesh ...\n\n");
#ifdef SIM
//  PM_write(pmesh, moutfile,0);
#endif
#ifdef FMDB
//  PM_write(pmesh, moutfile);
#endif
  wtimePoints[23] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);

  // create a new directory
  // put the adapted mesh and solution in this new directory
  char sss[100];
  sprintf(sss,"mkdir %d",lstep);
  system(sss);
  sprintf(sss,"%d",lstep);
  chdir(sss);

  MPI_Barrier(MPI_COMM_WORLD);
  wtimePoints[24] = time(0);

#ifdef DEBUG
  // write each part into a separate file
  // latest mesh file for each partition seems to have some trouble during visualization
  int version = 2; 
  int gproc;
  char stringName[255];
  if(PMU_rank()==0)
    printf("\n Writing out the partitions ...\n\n");
  for(  int i =0 ; i < PM_numParts(pmesh); i++ ) {
    gproc = PMU_totalPartFromGid( PM_totalNumParts(pmesh), PMU_gid( PMU_rank(), i )) ;
    mesh = PM_mesh(pmesh, i);
    sprintf(stringName, "refined_p.%d.sms",gproc+1);
//    M_write( mesh, stringName, version);
  }

#endif//DEBUG

  wtimePoints[25] = time(0);
  MPI_Barrier(MPI_COMM_WORLD);

  wtimePoints[26] = time(0);

  if(strategy>0 && (strategy!=7 || strategy!=8)) {

//     int nshg_fine = M_numVertices(mesh);

//     double* solution;
//     getAttachedArray(solution,mesh,phasta_solution,ndof,poly);

//     char output_file[256];
//     sprintf(output_file,"restart.%i.%i",lstep,gproc);
//     writeArrayToFile(output_file,"solution",outputFormat,"write",
// 		     nshg_fine,ndof,lstep,solution);
  
//     delete [] solution;
  }
  
  wtimePoints[27] = time(0);

  wtimePoints[28] = time(0);

  if(BCInflowFace) {
    if(PMU_rank()==0)
      printf(" Info: printing inflow face info. for model face with tag : %d\n",BCInflowFaceTag);
    BCInflowFaceInfo(model,pmesh,mesh);
  }

  wtimePoints[29] = time(0);
#ifdef DEBUG
  if( (option>0 && strategy==1) ||
      (strategy==5 || strategy==6) ) {
    printf("\n[%d] Info: printing timing statistics in log file \n",PMU_rank());
    printTimeStatsToFile(strategy);
  }
#endif

  //if(!(strategy==7 || strategy==8)) { //CWS
  if(!(strategy==7) ) {
    AttachDataCommu_delete(adc);
    AttachDataCommu_delete(adcSN);
  }

  MD_deleteMeshDataId(incorp);
  MD_deleteMeshDataId(errorIndicatorID);
  MD_deleteMeshDataId(modes);

  MD_deleteMeshDataId(nodalGradientID);
  MD_deleteMeshDataId(nodalHessianID);
  MD_deleteMeshDataId(localPatchVolID);
  MD_deleteMeshDataId(localGradientID);
  MD_deleteMeshDataId(localHessianID);
  MD_deleteMeshDataId(numSurroundNodesID);
  MD_deleteMeshDataId(locMaxInterpolErrorID);
  MD_deleteMeshDataId(globMaxInterpolErrorID);  

  MD_deleteMeshDataId(interfaceMetricID);   //CWS
  nshgTot = getNSHGTOT(mesh);
  nRgnsLoc = M_numRegions(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nRgnsLoc, &nRgnsTot, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

  if(PMU_rank()==0) {
    printf("\nTotal number of nodes   (after adaptation): %d\n",nshgTot);
    printf("\nTotal number of regions (after adaptation): %d\n",nRgnsTot);
  }



  return nshgTot;
} 

void printTimeStatsToFile(int strategy) 
{
  int gproc = PMU_gid(PMU_rank(),0);
  
  char tStatsFile[256];
  sprintf(tStatsFile,"timeStats.log.%i",gproc+1);
  FILE *tFile=fopen(tStatsFile,"w");

  fprintf(tFile,"\n");
  fprintf(tFile,"\nWall Time Stats for Partition [%d] (in secs):\n\n",gproc+1);
  fprintf(tFile,"PM_read()                        : %6.2f\n",
	  difftime(wtimePoints[1],wtimePoints[0]));
  fprintf(tFile,"PM_mesh()                        : %6.2f\n",
	  difftime(wtimePoints[3],wtimePoints[2]));
  fprintf(tFile,"Data migr. set-up                : %6.2f\n",
	  difftime(wtimePoints[5],wtimePoints[4]));
  fprintf(tFile,"MSA_new() & MSA_setCallBack() : %6.2f\n",
	  difftime(wtimePoints[7],wtimePoints[6]));
  fprintf(tFile,"readArrayFromFile()              : %6.2f\n",
	  difftime(wtimePoints[9],wtimePoints[8]));

  if(strategy==1) {
    fprintf(tFile,"hessiansFromSolution()           : %6.2f\n",
	    difftime(wtimePoints[11],wtimePoints[10]));
    fprintf(tFile,"setSizeFieldUsingHessians()      : %6.2f\n",
	    difftime(wtimePoints[13],wtimePoints[12]));
  }
  else if(strategy==5 || strategy==6) {
    fprintf(tFile,"transformToScalarErrorVal()      : %6.2f\n",
	    difftime(wtimePoints[11],wtimePoints[10]));
    fprintf(tFile,"setIsotropicSizeField()          : %6.2f\n",
	    difftime(wtimePoints[13],wtimePoints[12]));
  }
  
  fprintf(tFile,"*********************************************\n");
  fprintf(tFile,"MSA_adapt()                      : %6.2f\n",
	  difftime(wtimePoints[15],wtimePoints[14]));
  fprintf(tFile,"MSA_adapt() [CPU Time]           : %6.2f\n",
	  difftime((cputimePoints[1].ru_utime).tv_sec+
		   (cputimePoints[1].ru_stime).tv_sec,
		   (cputimePoints[0].ru_utime).tv_sec+
		   (cputimePoints[0].ru_stime).tv_sec));
  fprintf(tFile,"*********************************************\n");
  fprintf(tFile,"fix4NodesOnSurface()             : %6.2f\n",
	  difftime(wtimePoints[17],wtimePoints[16]));
  fprintf(tFile,"fix4SolutionTransfer()           : %6.2f\n",
	  difftime(wtimePoints[19],wtimePoints[18]));
  if(masterProcWgt>0.) {
    fprintf(tFile,"PM_partition() (AA)              : %6.2f\n",
	    difftime(wtimePoints[21],wtimePoints[20]));
  }
  fprintf(tFile,"PM_write()                       : %6.2f\n",
	  difftime(wtimePoints[23],wtimePoints[22]));
  fprintf(tFile,"M_writeSMS()                     : %6.2f\n",
	  difftime(wtimePoints[25],wtimePoints[24]));
  fprintf(tFile,"writeArrayToFile()               : %6.2f\n",
	  difftime(wtimePoints[27],wtimePoints[26]));
  fprintf(tFile,"BCInflowFaceInfo()               : %6.2f\n",
	  difftime(wtimePoints[29],wtimePoints[28]));

  fclose(tFile);
}

void check(pMesh mesh)
{
  int nv=0, ne=0, nf=0, nr=0, unknown=0;
  
  VIter vit=M_vertexIter(mesh);
  pVertex vertex;
  while( vertex=VIter_next(vit) ) {
    switch( V_whatInType(vertex) ) {
    case Tvertex: nv++; break;
    case Tedge: ne++; break;
    case Tface: nf++; break;
    case Tregion: nr++; break;
    default:
      unknown++;
    }
  }
  VIter_delete(vit);

  printf("# of vertices on gv: %d\n",nv);
  printf("# of vertices on ge: %d\n",ne);
  printf("# of vertices on gf: %d\n",nf);
  printf("# of vertices on gr: %d\n",nr);
  printf("# of vertices of unknown classification: %d\n",unknown);
}

/**
 * @brief     Sets the size field using a size field attached to the mesh.
 * @remark    The 'mesh' vertices are iterated over and the the attached data with identity 'attachedDataID' is read.  It is assumed that this data is a 3x3 matrix conforming to the Simmetrix format for the metric.  Please refer to the Simmetrix SimModSuite MeshAdapt documentation for this format specification.
 *    <-- INSERT BLANK LINE -->
 * @param     mesh (In) part mesh
 * @param     simAdapter (InOut) adaptation object
 * @param     attachedDataID (In) identifier of attached data 
 * @return
 */
int setSizeFieldFromAttachedData(pMesh mesh, 
                                 pMSAdapt simAdapter, 
                                 pMeshDataId attachedDataID) 
{     
  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
      //retrieve the interface metric
      double *vtxInterfaceMetric;
      if(!EN_getDataPtr((pEntity)vertex, attachedDataID,(void**)&vtxInterfaceMetric)){
         cout<<"\nERROR in setSizeFieldFromAttachedData(...) : requested data is not attached to vertex" << endl;
         V_info(vertex);
         exit(0);
      }

      double simMetric[3][3];

      for(int iRow=0; iRow<3; iRow++) 
      {
         for(int iDir=0; iDir<3; iDir++) 
         {
            simMetric[iRow][iDir] = vtxInterfaceMetric[iRow*3+iDir];
         }
      }
      MSA_setAnisoVertexSize(simAdapter, vertex, simMetric); //CWS
  }
  VIter_delete(vIter);
}

#ifdef __cplusplus
}
#endif

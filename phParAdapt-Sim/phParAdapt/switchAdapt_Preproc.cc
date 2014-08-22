///////////////////////////////////////////////////////////////////////////
// switchAdapt_Preproc.cc
// switch adaptor or preprocessor 
//
// in case the adaptor is called ahead, mesh construction
// is to be handled differently
//
// also a datastructure has to be provided for the solution transfer
//
// J.Mueller/O.Sahni 2004-2007
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include "phParAdapt.h"
#include "func.h"
#include "ccfunc.h"

#include "MeshSimInternal.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#include "SimAdvMeshing.h"
#include "SimMeshTools.h"
#endif
#ifdef FMDB
#include "myAttribute.h"
#include "pmZoltanCallbacks.h"
#include "phPartitionCallbacks.h"
#include "AOMD_cint.h"
#include "mAOMD.h"
#include "pmZoltanCallbacks.h"
#ifdef ENTITY_GROUP
#include "AOMD_EntGrp.h"
#include "mEntityGroup.h"
#include "BLUtil.h"
#endif
#endif
//#include "M_writeVTKFile.h"
#include "MeshSimAdapt.h"
#include "phReadWrite.h"
#include "attachData.h"
#if ( defined SIM_PARASOLID )
  #include "SimParasolidKrnl.h"
#elif ( defined PARASOLID )
    #include "modelerParasolid.h"
#elif ( defined DISCRETE )
  #include "SimDiscrete.h"
#endif

#ifdef FMDB
  using namespace AOMD;
#endif

extern "C" int procSize();

#ifdef __cplusplus
extern "C" {
#endif

extern int adaptFlag;
extern int timeStepNumber;
extern int rStart;
extern int multipleRestarts;
extern int strategy;
extern int ensa_dof;// number of field(=solution) variables
extern int nSolVars;
extern int nErrorVars;
extern double* wght;
extern double factor;
extern double hmax;
extern double hmin;
extern int adaptOption;
extern int numTotParts;
extern int numParts;
extern int preLBforAdaptivity;
extern double masterProcWgt;
extern int SolutionMigration;
extern int DisplacementMigration;
extern pMeshDataId phasta_solution;
extern char gname[256];
pProgress prog;
int lstep;
extern time_t wtimePoints[32];
extern int isBLAdapt;
extern int isThickAdapt;
extern int dwalMigration;
extern int buildMapping;

//pAManager SModel_attManager(pModel model);

int deleteDblArray(void*, pAttachDataId, int, void **attdata, void*)
{
     // Note: attdata is the pointer to the pointer that was attached
     delete static_cast<double*>(*attdata);
       return 1;                                      // 1 for success
}

int
switchAdapt_Preproc(int argc, char *argv[]){

    // process command line arguments
    // procArgs overwrites geom.sms by the mesh directory name specified
    // fname:attribute file name ,mname mesh directory name 
    procArgs(argc, argv);

    assignGlobalVars();

    lstep = timeStepNumber;

    // model and mesh declaration
    pGModel model; 
    pParMesh pmesh;
    pMesh mesh;
    int ndof = nSolVars, poly = 1;
    int ndwal = 0;
    int nmapping = 0;
    int ndisp = 0;

    if(!adaptFlag)
       phasta_solution = MD_newMeshDataId("restart solution");

//    pNativeModel nModel; //This variable is not used Min Zhou
    pACase acase;
    char    fname[100];
    char    mname[100];
  // default attribute/mesh filename
    strcpy(fname, "geom.spj");

    if((adaptFlag && (strategy!=7 || strategy!=8))|| multipleRestarts ){
        strcpy(fname, "geomNOIC.spj");
    }

    strcpy(mname, "geom.sms");

/*
//spj file stuff   

    pAManager attmngr = AMAN_load(fname);  // - guess that this replaces AMAN_retrieve()
    
    if (attmngr == 0){
      if (PMU_rank() == 0){
	fprintf(stderr, "could not open attribute file %s\n", fname);
	exit(-1);
      }           
    }
    else{
      if (PMU_rank() == 0){
	printf("\n AMAN_load(fname) success\n");
      }
    }
    acase = AMAN_findCase(attmngr, "geom");
    if (acase == NULL){
      char* casename="geom";
      if ( ( acase = AMAN_findCase(attmngr, casename)) == NULL) {
	if(PMU_rank() == 0) {
          printf("[%d] WARNING: could not find attribute case %s\n",PMU_rank(),casename);
        }
        casename="problem definition";
	if(PMU_rank() == 0) {
          printf("[%d] Attempting to find case %s\n",PMU_rank(),casename);
        }
        if ( ( acase = AMAN_findCase(attmngr, casename)) == NULL) {
	   if(PMU_rank() == 0) {
             printf("[%d] Error: could not find attribute case %s\n",PMU_rank(),casename);
           }
           exit (-1);
        }
      }
      else{
	if(PMU_rank() == 0)
	  printf("\n AMAN_findCase success\n");
      }
    }
//end spj stuff
*/

    // associate the attribute case with the model
#if  ( defined SIM_PARASOLID )
   if(PMU_rank()==0) printf("\n model file = %s\n", gname);
   MPI_Barrier(MPI_COMM_WORLD);
   pParasolidNativeModel pnModel;
   pnModel = ParasolidNM_createFromFile(gname,0);  	// create Parasolid native model from part file

   if(NM_isAssemblyModel(pnModel)) {
     pGAModel amodel = GAM_createFromNativeModel(pnModel, prog);
     NM_release(pnModel);
     model = GM_createFromAssemblyModel(amodel, NULL, prog);
     GM_release(amodel);
   }
   else {
//    model = GM_createFromNativeModel(pnModel, prog);
    model = GM_load("geom.smd", pnModel, prog);
    NM_release(pnModel);
   }
#endif  
#if (defined SIM_GEOMSIM)
           model = GM_load("geom.smd", NULL, prog);
#endif

//New way with smd
    pAManager attmngr = SModel_attManager(model);
    acase = AMAN_findCase(attmngr, "geom");
    if (acase){
       if(PMU_rank()==0) printf("Found case, setting the model\n");
       MPI_Barrier(MPI_COMM_WORLD);
       AttCase_setModel(acase,model);
    } else {
        printf("Case not found, rename case to geom\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(PMU_rank()==0) printf("AttCase_associate(acase,prog)\n");
    AttCase_associate(acase,prog);
    MPI_Barrier(MPI_COMM_WORLD);
    if(PMU_rank()==0) printf("AttCase_associate(acase,prog); ... done\n");

//New smd end
//*/

/*
//Old way with spj files
    AttCase_setModel  ( acase, model );   

    // associate the case with the model (case: e.g. set of BCs)
    AttCase_associate(acase, prog);
//spj file end    
*/

#ifdef PARALLEL

    if(PMU_rank()==0) {
      printf("\n");
      printf("Meshing    Library build ID : %s\n",SimMeshing_buildID());
      printf("PMesh      Library build ID : %s\n",SimPartitionedMesh_buildID());
      printf("MeshTools  Library build ID : %s\n",SimMeshTools_buildID());
      printf("AdvMeshing Library build ID: %s\n",SimAdvMeshing_buildID());
    }
    
    double tmptime;
    wtimePoints[0] = time(0);
#ifdef SIM
    pPartitionOpts pOpts;
    pAttachDataCommu adc;
    int id; 

    if(PMU_size()==1) {
      pmesh = PM_load(mname,sthreadDefault,model,prog);
      if(adaptFlag && isThickAdapt)
         mesh = M_load(mname, model, prog);
//         mesh = PM_mesh(pmesh, 0);
      else 
         mesh = PM_mesh(pmesh, 0);
    } else {
      pmesh = PM_load("parts.sms",sthreadDefault,model,prog);
      mesh = PM_mesh(pmesh, 0);
    }
#endif

//   VolumeDebug(mesh);
     EdgeLength(mesh);
#ifdef FMDB
    if(numTotParts < PMU_size())
        cout<<"ERROR: numTotParts < nproc, only numTotParts>=nproc is supported \n";

    pmesh = MS_newMesh(model);
    vector<pParMesh> meshes;

    PM_load(pmesh,mname);
    mesh = pmesh;

    tmptime = time(0);
    if(PMU_rank()==0)
      cout<<"Time Usage of PM_load is: "<<tmptime-wtimePoints[0]<<"\n";	
 
#ifdef ENTITY_GROUP
    if(isBLAdapt)
      Mesh_InitBLs(pmesh, model);
#endif
#endif

    if(PM_verify(pmesh, 1, sthreadDefault, prog) == 0)
       cout << "Mesh is not valid!!!! " << endl;

    if(SolutionMigration && !adaptFlag) //migrate solution from restart files
    {
        ifstream numstart("numstart.dat");
        int stepnum;
        numstart>>stepnum;
        char solutionfile[64];
        sprintf(solutionfile,"restart.%d.%d", stepnum,PMU_rank()+1);
        double *sol, *disp;
        readArrayFromFile(solutionfile,"solution", sol);
        ndisp = 0;
        if(DisplacementMigration)
        {
            readArrayFromFile(solutionfile,"displacement", disp); 
            ndisp = 3;
        } 
        double* dwal;
        ndwal = 0;
        if(dwalMigration) {
          readArrayFromFile(solutionfile,"dwal",dwal);
          ndwal = 1;
        }

        int nshg = M_numVertices(mesh);

        double* mapping;
        nmapping = 0;
        if(buildMapping && !adaptFlag) {
           nmapping = 2;
           mapping = new double[nshg*nmapping];
           for(int inode=0;inode<nshg;inode++){
             mapping[inode*nmapping]   = (double) PMU_rank(); //input part (starting from 0)
             mapping[inode*nmapping+1] = (double) inode; //vtx number (starting from 0)
           } 
        }

        int totsize = (ndof+ndisp+ndwal+nmapping);
        double *SolTot = new double[totsize*nshg];

        for(int inode=0;inode<nshg;inode++) {
           for(int idof=0;idof<ndof;idof++)
              SolTot[inode*totsize+idof] = sol[inode*ndof+idof];
           for(int idof=0;idof<ndisp;idof++)
              SolTot[inode*totsize+ndof+idof] = disp[inode*ndisp+idof];
           for(int idof=0;idof<ndwal;idof++)
              SolTot[inode*totsize+ndof+ndisp+idof] = dwal[inode*ndwal+idof];
           for(int idof=0;idof<nmapping;idof++)
              SolTot[inode*totsize+ndof+ndisp+ndwal+idof] = mapping[inode*nmapping+idof];
        }       

        attachArray(SolTot, mesh, phasta_solution,totsize, poly);
        delete[] sol;
        delete[] SolTot;
        if(DisplacementMigration)
            delete[] disp;
        if(dwalMigration)
          delete[] dwal;
        if(buildMapping && !adaptFlag)
          delete[] mapping;

        numstart.close();        
#ifdef DEBUG        
//        M_writeVTKFile(mesh, "InitialErrors", phasta_solution, ndof);
//        M_writeVTKFile(mesh, "InitialSolution", phasta_solution, ndof);
#endif        
    }
#ifdef SIM
    if(SolutionMigration!=0 && !adaptFlag) {
      
       int mult = sizeof(double)/sizeof(int);
       //adc = PM_newAttachDataCommu(mult, 0, ndof);
       adc = PM_newAttachDataCommu(mult, 0, ndof+ndisp+ndwal+nmapping);

       MD_setMeshCallback(phasta_solution, CBmigrateOut, pm_sendDblArray, adc);
       MD_setMeshCallback(phasta_solution, CBmigrateIn, pm_recvDblArray, adc);
       MD_setMeshCallback(phasta_solution, CBdelete, deleteDblArray, 0);

       PM_setMigrId(pmesh, phasta_solution);

    }
    pOpts = PM_newPartitionOpts();
    PartitionOpts_setTotalNumParts(pOpts, numTotParts);

    if(numTotParts!=PMU_size() && !adaptFlag) { //multiparts per proc
      PM_partition(pmesh, pOpts, sthreadNone, prog);
      PartitionOpts_delete(pOpts);
    }
    numParts = PM_numParts(pmesh);


    if(PMU_rank()==0) {
      printf("checking the mesh after partitioning");
    }
    if(PM_verify(pmesh, 1, sthreadDefault, prog) == 0)
        cout << "Mesh is not valid!!!! " << endl;
#endif

#ifdef FMDB
    meshes.push_back(pmesh);
    
    if(numTotParts>PMU_size()&&!adaptFlag){ //multiparts per proc
        numParts = numTotParts/PMU_size();
        int remainPart = numTotParts - numParts*PMU_size();
        if(PMU_rank()<remainPart)
            numParts++;
        PM_setMaxNumParts(numParts);        
        tmptime = time(0);
        if(SolutionMigration){
            ndisp = 0;
            if(DisplacementMigration){
                ndisp = 3;
            }
            ndwal = 0;
            if(dwalMigration) {
              ndwal = 1;
            }
            nmapping = 0;
            if(buildMapping && !adaptFlag) {
              nmapping = 2
            }
            phPartitionCallbacks zlb(phasta_solution, ndof+ndisp+ndwal+nmapping);
    	    zlb.setLocalNumParts(numParts);
            M_loadbalance2(meshes,zlb);
        }
        else{
            AOMD::zoltanCB zlb;
            zlb.setLocalNumParts(numParts);
            M_loadbalance2(meshes, zlb);
        }
       if(PMU_rank()==0)
          cout<<"Time Usage of M_loadbalance2 is :"<<time(0)-tmptime<<"\n";
    }
    else 
        PM_setMaxNumParts(1);
    wtimePoints[1] = time(0);
    PM_write2(meshes, "geom_.sms");

    if(PM_verify(pmesh, 1, sthreadDefault, prog) == 0)
          cout << "Mesh is not valid!!!! " << endl;
          
 /* 
    if(PM_verify  (  pmesh ,0,sthreadDefault,prog ) == 0){
      if (PMU_rank() == 0){
	printf("\nerror in adapt.cc: invalid parallel mesh read in\n");
      }
      SimPartitionedMesh_stop();
      exit(1);
    }
 */   
#endif
#endif
//to calculate minimum and maximum edge lengths

/*
//checking if pyramids are part of BL
     RIter rit=M_regionIter(mesh);
     pRegion region;
     int Count = 0;
     while(region=RIter_next(rit) ) {
        int rType=R_topoType(region);
        if (rType==3) {
           if(!EN_isBLEntity(region)) {
              printf("Pyramid not part of the BL!!!\n");
              Count++;
           }
        }
     }
     printf("Count of pyramids: %d\n",Count);
     RIter_delete(rit);
     exit(0);
            //
*/
    // use the adaptor first
    // created new mesh, new restart files
    // total number of DOFS defined here
    int nshgTot=0;

    if(adaptFlag){        
#ifndef ibm
//        printf("\n[%2d] memory usage before mesh adaptation: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
        // flags needed by the preprocessor
        // to directlty take the ICs from the restart files
        // reads files: ./restart.%d.%d
        rStart=1;
        nshgTot = adapt (pmesh,
                         mesh,
                         model,
                         timeStepNumber,
                         strategy,
                         factor,
                         // number of field(=solution) variables
                         ensa_dof,
                         // number of variables for error indicators (EI)
                         // (e.g., 5 for ybar & 10 for residual-based)
                         nErrorVars,
                         hmax,
                         hmin,
                         adaptOption);


    if(PM_verify(pmesh, 1, sthreadDefault, prog) == 0)
        cout << "Mesh is not valid!!!! " << endl;
#ifndef ibm
  //      printf("\n[%2d] memory usage after mesh adaptation: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
#ifdef FMDB
        meshes.clear();
        meshes.push_back(pmesh);

        if(numTotParts>PMU_size()){ //multiparts per proc
            numParts = numTotParts/PMU_size();
            PM_setMaxNumParts(numParts);        
            ndisp = 0;
            if(DisplacementMigration)
                ndisp = 3;
            ndwal = 0;
            if(dwalMigration) {
              ndwal = 1;
            }
            nmapping = 0;
            if(buildMapping && !adaptFlag) {
              nmapping = 2
            }
            phPartitionCallbacks zlb(phasta_solution, ndof+ndisp+ndwal+nmapping;
    	    zlb.setLocalNumParts(numParts);
            M_loadbalance2(meshes,zlb);
        }
#endif
    }

#ifdef SIM

    PM_write(pmesh, "mesh_parts.sms", sthreadNone, prog);
#endif

    if(nErrorVars)
      delete [] wght; 
    
    vector<pParMesh> meshes;
    meshes.push_back(pmesh);

    // continue with the usual preprocessing
    mdb2phasta (fname,mname,model, meshes,nshgTot );
#ifndef ibm
    //printf("\n[%2d] memory usage after running preprocessor: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif

    AMAN_release(attmngr);
    /* Delete mesh */
    for(int ipart=0;ipart<numParts;ipart++)
        //Mesh need to be release here Min Zhou
//        M_release((pUnstructuredMesh) meshes[ipart]);

    /* Delete the model. */
        //Model need to be recleased here Min Zhou
//    GM_release(model);
    
    return 1;
    
}


#ifdef __cplusplus
}
#endif   

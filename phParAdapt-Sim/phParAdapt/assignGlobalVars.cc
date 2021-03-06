#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <iostream> 
#include <fstream>
#include <string>

#include "Input.h"
#include "phParAdapt.h"
using namespace std;
/* global variables */
extern double sizeRatio=0.5;
extern double ratioThresh=0.75;
extern int localAdapt=0;
// 0: all mesh is adapted
// 1: only where sizes set are modified
extern int coarsenMode=2;
// 0: No coarsening will be done on the mesh.
// 1: By default an initial coarsening pass is done. Mode 1 disables this pass.
// 2: All coarsening is done as needed (default).
extern int numSmooth=10;
extern int AnisoSimmetrix=1;
extern int globalP;
extern int fType;
extern int rStart;
extern int per;
extern int intBC;
extern int timing;
extern int wGraph;
extern int prCd;
extern int  zScale;
int WRITEASC=0;
int ensa_dof=5;
int CUBES=0;
int rRead=5;
int NSFTAG=-1; 
int lStart=0; 
int FortFormFlag=0; 
int phastaIO=1;
int numParts=1;
int numTotParts = 0;
int SolutionMigration = 0;
int DisplacementMigration = 0;
int isReorder = 0;
int numSplit = 0;
int isBLAdapt = 0;
int isSizeLimit = 0;
int MaxLimitFact = 16;
int MinLimitFact = 16;
int isThickAdapt = 0;
double rhoinp = 1.225;
double muinp = 1.78e-5;
int isMCOPI = 0;
int old_format=0;

int dwalMigration = 0;
int buildMapping = 0;

char mname[256];
char fname[256];
char gname[256];

// for adaptation 
int strategy;
int nSolVars;
int nErrorVars;
int nYbarVars;
double factor;
double hmax;
double hmin;
int adaptOption;
double refthreshold;
char outputFormat[100];
int adaptFlag=0;//use adaptor first
double* wght;
double epsilon = 0.001;
double max_size = 0.2;
double min_size = 0.1;
int size_flag = 0;

int preLBforAdaptivity;
double masterProcWgt;

int timeStepNumber=0;

char* iformat = "binary";
char* oformat;

extern int justpart;
extern int parted;
extern int errorName;
// should be removed in later versions
// now still used in attachIBC.cc
// the issue is whether the option should be kept
// to be able to read in a SINGLE restart file
// or if that should be moved to the partitioner
int multipleRestarts = 0;
int numRestartsToRead = 1;

char phversion[256];
char meshoutname[256];
char phastaVersion[256];

int numelX,Idirection;
int BYPASS=0;
int SONFATH = 0;

void 
assignGlobalVars()
{
    ifstream inf("adapt.inp");
    if (!inf) {
        cerr << "Error: adapt.inp not found" << endl;
        exit(0);
    }
    numTotParts = PMU_size();
   
    Input inputObject(inf);

//******************my changes, I am reading files adapt.inp Azat 09.16.04 12.58PM
//******************Help Part***************************************************

    refthreshold  = inputObject.threshold(); // reference threshold value(fraction of GlobalMaxError) above which the edges would be marked
    nSolVars = inputObject.inp_nSolVars();
    nErrorVars   = inputObject.inp_nErrorVars();//Number of error indicators to be read
    nYbarVars   = inputObject.inp_nYbarVars();//Number of error indicators to be read
    double *wght1 = inputObject.weights();   // Weights of the error indicators
    wght = new double[nErrorVars];
    for(int iEvar=0; iEvar<nErrorVars; iEvar++)
      wght[iEvar] = wght1[iEvar];
    //sets the global polynomial order
    globalP = inputObject.inp_globalP();
    numSmooth = inputObject.inp_numSmooth();
    AnisoSimmetrix = inputObject.inp_AnisoSimmetrix();
    localAdapt = inputObject.inp_localAdapt();
    coarsenMode = inputObject.inp_coarsenMode();
    ratioThresh = inputObject.inp_ratioThresh();
    sizeRatio = inputObject.inp_sizeRatio();

    // time step number
    timeStepNumber = inputObject.inp_timeStepNumber();
    numelX = inputObject.inp_numelX();
    NSFTAG = inputObject.inp_NSFTAG();
    ensa_dof =inputObject.inp_ensa_dof();
    strcpy(fname,inputObject.inp_fname());
    strcpy(mname,inputObject.inp_mname());
    strcpy(gname,inputObject.inp_gname() ) ;
    Idirection = inputObject.inp_Idirection();
    BYPASS  =inputObject.inp_BYPASS();
    fType  = inputObject.inp_fType();
    zScale =  inputObject.inp_zScale();
    adaptFlag  = inputObject.inp_adaptFlag();
    errorName  = inputObject.inp_errorName();
    SONFATH  = inputObject.inp_SONFATH();
    lStart  = inputObject.inp_lStart();
    rRead   = inputObject.inp_rRead();
    rStart   = inputObject.inp_rStart();
    strategy   = inputObject.inp_strategy();
    factor   = inputObject.inp_factor();
    
    nSolVars   = inputObject.inp_nSolVars();
    nErrorVars   = inputObject.inp_nErrorVars();
    hmax   = inputObject.inp_hmax();
    hmin   = inputObject.inp_hmin();
    adaptOption = inputObject.inp_adaptOption ();
    multipleRestarts   = inputObject.inp_multipleRestarts();
    per   = inputObject.inp_per();
    prCd   = inputObject.inp_prCD();
    timing   = inputObject.inp_timing();
    wGraph   = inputObject.inp_wGraph();
    strcpy(phastaVersion , inputObject.inp_phastaVersion());
    old_format   = inputObject.inp_old_format();
    FortFormFlag   = inputObject.inp_FortFormFlag();
    strcpy(outputFormat ,inputObject.inp_outputFormat());

    oformat = outputFormat; 

    CUBES   = inputObject.inp_CUBES();
    intBC   = inputObject.inp_intBC();
    strcpy(phversion , inputObject.inp_version());
    WRITEASC   = inputObject.inp_WRITEASC();
    phastaIO   = inputObject.inp_phastaIO();
    parted     =  inputObject. inp_parted();
    // changes from Joe
    epsilon  = inputObject.inp_epsilon();
    max_size = inputObject.inp_max_size();
    min_size = inputObject.inp_min_size();
    size_flag = inputObject.inp_size_flag();

    preLBforAdaptivity =  inputObject. inp_preLBforAdaptivity();
    masterProcWgt         =  inputObject. inp_masterProcWgt();
    numTotParts = inputObject.inp_numTotParts();
    SolutionMigration = inputObject.inp_SolutionMigration();
    DisplacementMigration = inputObject.inp_DisplacementMigration();
    isReorder = inputObject.inp_Reorder();
    numSplit = inputObject.inp_numSplit();
    isBLAdapt = inputObject.inp_BLAdapt();
    isSizeLimit = inputObject.inp_isSizeLimit();
    MaxLimitFact = inputObject.inp_MaxLimitFact();
    MinLimitFact = inputObject.inp_MinLimitFact();
    isThickAdapt = inputObject.inp_ThickAdapt();
    rhoinp = inputObject.inp_rhoinp();
    muinp = inputObject.inp_muinp();
    isMCOPI = inputObject.inp_MCOPI();

    dwalMigration = inputObject.inp_dwalMigration();
    buildMapping = inputObject.inp_buildMapping();
}

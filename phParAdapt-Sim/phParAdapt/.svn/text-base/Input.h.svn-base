#ifndef H_Input
#define H_Input

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;

class Input {
public:
  Input (ifstream &inf);
  ~Input () {delete [] d_wght;}
  int inp_nSolVars () { return dnSolVars;}
  int inp_nErrorVars () { return dnErrorVars;}
  int inp_nYbarVars () { return dnYbarVars;}
  double *weights () { return d_wght; }
  double threshold() { return d_threshold; }
//this input is requiring for new release of phParAdapt (with size+tag field adaptivity)
int inp_globalP () {return dglobalP;}
int inp_timeStepNumber () {return dtimeStepNumber;}
int inp_numelX (){return dnumelX;}
int inp_NSFTAG (){return dNSFTAG;}
int inp_ensa_dof (){return densa_dof;}
char *inp_fname () {return dfname;}
char *inp_mname () {return dmname;}
char *inp_gname () {return dgname;}
int inp_Idirection () {return dIdirection;}
int inp_BYPASS  () {return dBYPASS;}
int inp_fType () {return dfType;}
double inp_zScale () {return dzScale;}
int inp_adaptFlag () { return dadaptFlag;}
int inp_errorName () { return derrorName;}
int inp_SONFATH   () { return dSONFATH;}
int inp_lStart    () { return dlStart;}
int inp_rRead     () { return drRead;}
int inp_rStart    () { return drStart;}
int inp_strategy  () { return dstrategy;}
double inp_factor () { return dfactor;}
double inp_hmax   () { return dhmax;}
double inp_hmin   () { return dhmin;}
int inp_adaptOption ()  { return dadaptOption ;}
int inp_multipleRestarts () { return dmultipleRestarts;}
int inp_per () { return dper;}
int inp_prCD () {return dprCD;}
int inp_timing () { return dtiming;}
int inp_wGraph () { return dwGraph;}
char *inp_phastaVersion () { return dphastaVersion;}
int inp_old_format () { return dold_format;}
int inp_FortFormFlag () { return dFortFormFlag;}
char*  inp_outputFormat () { return doutputFormat;}
int inp_CUBES () { return dCUBES;}
int inp_intBC () { return dintBC;}  
char* inp_version () { return dversion;}
int inp_WRITEASC () { return dWRITEASC;}
int inp_phastaIO () { return dphastaIO;}  
int inp_parted() {return dparted;}
int inp_preLBforAdaptivity() {return dpreLBforAdaptivity;}
double inp_masterProcWgt() {return dmasterProcWgt;}
//additions from Joe
double inp_epsilon () {return depsilon;}
double inp_max_size () {return dmax_size;}
double inp_min_size () {return dmin_size;}
int inp_size_flag () {return dsize_flag;}
 int inp_numTotParts(){return dnumTotParts;}
 int inp_SolutionMigration(){return dSolutionMigration;}
 int inp_DisplacementMigration(){return dDisplacementMigration;}
 int inp_Reorder(){return disReorder;}
 int inp_numSplit(){return dnumSplit;}
 int inp_BLAdapt(){return disBLAdapt;}
 int inp_isSizeLimit(){return disSizeLimit;}
 int inp_MaxLimitFact(){return dMaxLimitFact;}
 int inp_MinLimitFact(){return dMinLimitFact;}
 int inp_ThickAdapt(){return disThickAdapt;}
 double inp_rhoinp   () { return drhoinp;}
 double inp_muinp   () { return dmuinp;}
 int inp_dwalMigration() {return ddwalMigration;}
 int inp_buildMapping() {return dbuildMapping;}
 int inp_MCOPI() {return disMCOPI;}
private:
 //************this part for the local variables*****************
double depsilon;
double dmax_size;
double dmin_size;
int dsize_flag;
 int dglobalP;
 int dtimeStepNumber;
 int dnumelX;
 int dNSFTAG;
 int densa_dof;
 char dfname[256];
 char dmname[256];
 char dgname[256];
 int dIdirection;
 int dBYPASS;
 int dfType;
 double dzScale;
 int dadaptFlag;
 int derrorName;
 int dSONFATH;
 int dlStart;
 int drRead;
 int drStart;
 int dstrategy;
 double dfactor;
 int dnSolVars;
 int dnErrorVars;
 int dnYbarVars;
 double dhmax;
 double dhmin;
 int dmultipleRestarts;
 int dper;
 int dprCD; 
 int dtiming;
 int dwGraph;
 char dphastaVersion[256];
 int dold_format;
 int dFortFormFlag;
 char doutputFormat[256];
 int dCUBES;
 int dintBC;
 char dversion[256];
 int dWRITEASC;
 int dphastaIO;
 double *d_wght;
 double d_threshold;
 int dadaptOption;
 int dparted;
 int dpreLBforAdaptivity;
 double dmasterProcWgt;
 int dnumTotParts;
 int dSolutionMigration;
 int dDisplacementMigration;
 int disReorder;
 int dnumSplit;
 int disBLAdapt;
 int disSizeLimit;
 int dMaxLimitFact;
 int dMinLimitFact;
 int disThickAdapt;
 double drhoinp;
 double dmuinp;
 int disMCOPI;
 int ddwalMigration;
 int dbuildMapping;
};


#endif


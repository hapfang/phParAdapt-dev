//////////////////////////////////////////////////////////////////////////
// This function collects the information from the mesh database into 
// the arrays which will be used by phasta.                          
//
// called in mdb2phasta
//
// all operations are performed on the calling proc's local mesh (=partition)
///////////////////////////////////////////////////////////////////////////
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"
#include <iostream>
#include "EnsaArrays.h"
#include "EnsaParameters.h"
using namespace std;

/* externally defined global variables */
extern time_t tstart[MAXNT];
extern double eltime[MAXNT];
extern int per, ICtyp, intBC, fType, zScale, ensa_dof, rStart, lStart, rRead;
extern int numRestartsToRead;// in procArgs
extern int multipleRestarts;// in procArgs
extern int lstep;
extern int FortFormFlag;
extern int SONFATH;
extern int adaptFlag;
extern int strategy;
extern int numTotParts;
extern int numParts;
extern int DisplacementMigration;
extern int dwalMigration;
/* global variables */

#ifdef PARALLEL
extern void writeEnsaFiles(pParMesh pmesh, globalInfo *info, pPList bdry, int ipart)
#else
extern void writeEnsaFiles(pMesh mesh, globalInfo *info, pPList bdry)
#endif
{
#ifdef PARALLEL
  pMesh mesh = PM_mesh(pmesh, ipart);
#endif
  int i;
  int nshgReqst;
  char fname[30];

  //lstep = FortFormFlag ? 1 : 0;

  // info->nshg calc now moved to localInfo for efficiency

  EnsaParameters ePar(3, info->nen, info->nfaces, 0, ensa_dof, 7+ensa_dof,
                      ensa_dof+1, fType);
  /*
     These parameters changed with the addition of scalar equations.
     The formulas should key off of numVars (I see that ndof is used elsewhere
     in the code so...).
     numVars gets set to ensa_dof by an input variable from the command line
     numEBC = 12 changed to ensa_dof + 7 (12 is for theta)
     numNBC = 6 changed to ensa_dof + 1 (1 extra for turbulence wall)
  */


  EnsaArrays earr(info, ePar);

  int *ifath=0, *nsons=0;

  // write processor mapping arrays
#ifdef PARALLEL
  writeNCorp(pmesh, info, &ifath, &nsons, ipart);
#else
  writeNCorp(mesh, info, &ifath, &nsons);
#endif

  if(PMU_rank()==0&&ipart==numParts-1)
      cout<<"\n["<<PMU_rank()<<"] writeNCorp success\n";
  
#ifdef PARALLEL
  /* generate the local work array for parallel communication */
  genILWork(info, earr.getILWORK(),ipart);

  if(PMU_rank()==0&&ipart==numParts-1)
    cout<<"\n["<<PMU_rank()<<"] genILWork success \n";

#endif

  double zscale[3]={0.0, 0.0, 0.0};
  if (zScale) {
      MPI_Barrier(MPI_COMM_WORLD);// just to be on the safe side
      if (PMU_rank() == 0){
          printf("\nEnter the scale factor for the x, y, z coordinate: ");
          cin >> zscale[0];
          cin >> zscale[1];
          cin >> zscale[2];
      }
  }
  MPI_Barrier(MPI_COMM_WORLD);// just to be on the safe side
  MPI_Bcast(&zscale, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  double* qTot;

  int nv=ensa_dof;
  int ndisp = 0;
  if(DisplacementMigration)
      ndisp = 3;
  int ndwal = 0;
  if(dwalMigration) 
    ndwal = 1;
  int numnp=M_numVertices(mesh);
  int nrquest = ndisp+ndwal+rRead;

  if(lStart){ //lin-2-quad
    // for now we use this to get a quadratic solution
    // load in the auxilliary linear mesh
    /*    pMesh lmesh=M_new(1, model);
    M_load(lmesh, "geoml.sms");
    int numnpl=M_numVertices(lmesh);
    double* qtmpl = new double [ nv*numnpl ];
    restart(qtmpl, numnpl, rRead, &lstep);
    int nshg = numnp + M_numEdges(mesh);
    qtmp = new double [ nv*nshg ];
    lin2quad (mesh, lmesh, qtmpl, qtmp, nshg); */
  }
  else if (rStart && adaptFlag == 0){
    // use linear portion
    if (rStart == 3) {
      nshgReqst = numnp;
    } else if (rStart == 4) {
      // use quadratic portion
      nshgReqst = numnp + M_numEdges(mesh);
    }
    // file contains all modes (rStart==1)
    else{
        nshgReqst = info->nshgTot;
    }

    qTot = new double [ nv*(info->nshgTot) ];
    for (i=0; i< nv*(info->nshgTot); i++){
        qTot[i]=0;
    }
    //
    //  We have zeroed the qTot array because restart is going to fill
    //  only the portion of this array that we request it to based on
    //  our choice of nshgReqst and Read (below called nvReqst for
    //  consistency. Note that the array is shaped nv*nshgTot because
    //  later, higher order modes will need to index into this TOTAL
    //  solution array for all variables.  If they were not read it
    //  will find zero's there which we assume will be ok or replaced
    //  by attribute expressions. 
    //
    restart(qTot, info->nshgTot, nshgReqst, rRead, &lstep, mesh, ipart);
  }
//  else if (rStart>0 && adaptFlag==1 && strategy!=7){
   else if (rStart>0 && adaptFlag==1){
      
      // request LOCAL number of Shapefuns
      nshgReqst = info->nshg;
      qTot = new double [ (nv+ndisp+ndwal)*(info->nshg) ];

      for (i=0; i< (nv+ndisp+ndwal)*(info->nshg); i++){
        qTot[i]=0;
      }
      // greps/reads solution from adaptor/file: implemented in restart.cc
      restart(qTot, info->nshg, nshgReqst, nrquest, &lstep, mesh, ipart);
      if(PMU_rank() == 0){
            cout<<"\n  solution retriever from mesh: ADAPTED success ...\n";
      }
  }
  // adaptFlag!=0 while rStart hasn't been specified
  else if(rStart==0 && adaptFlag==1) {
      if(PMU_rank() ==0 ){
          cerr<<"\nerror in writeEnsaFiles(): wanted to use adapted restart data \n"
              <<" while restart option hasn't been specified - \n";
      }
      exit(1);
  }
  // rStart==0 AND adaptFlag==0
  else{
      qTot = new double [ nv*(info->nshg) ];
  }

  // gets the coordinates passed into  earr's **x
  // implemented in getData.cc
  getX(mesh, earr.getX());


  // per - flag is set on the commandline, declared in mdb2phasta
  // per=1 is default 
  if (per == 1){
      
      // only if master AND slave are on same proc/part
      // iper will include them
      // off-proc master/slave pairs are treated as
      // "a   communication" (see SST comment in attachPeriodicBC)
      attachPeriodicBC(mesh, info, earr.getIPER());
  }

  // assign the number of ess. BCs
  // earr is EnsaArrays=EnsaArrays(globalInfo,EnsaParameters)
  // earr's nBC, iBC, BC are instantiated to
  // arrays of size info->nshg, value 0
  info->numpbc = attachEssentialBC(mesh, info, qTot,
                                   earr.getNBC(),
                                   earr.getIBC(),
                                   earr.getBC());//nBC,iBC,BC

#ifdef DEBUG
  printf("\n[%d] in writeEnsaFiles:number of essential BCs: %d \n", PMU_rank()*numParts+ipart,info->numpbc );       
#endif

  attachNaturalBC(mesh, earr.getIBCB(), earr.getBCB(), ePar.getNUMNBC());


  // earr.q (getQ()) was set in 
  attachInitialCondition(mesh, info, qTot, earr.getQ());

  if (qTot) delete [] qTot;

  getConnectivity(mesh, bdry, earr.getIEN(), earr.getIENB(), info);

  // writing the EnsaArrays into files (also uses phastaIO)
  earr.write(mesh, info, ePar, zscale, ifath, nsons, ipart);

//    if (!FortFormFlag)  {
//      /* we need to include the ncorp into the geombc database */
//      /* Now saving this array (info->ncvec) (it's made in writeNCorp) */
//      char geomstr[20];
//      sprintf(geomstr, "geombc.dat.%d", PMU_rank()+1);
//      if (!FortFormFlag)  {
//      /* we need to include the ncorp into the geombc database */
//      /* Now saving this array (info->ncvec) (it's made in writeNCorp) */
//      char geomstr[20];
//      sprintf(geomstr, "geombc.dat.%d", PMU_rank()+1);
//      FILE *fileptr = fopen(geomstr, "a");
//      fprintf(fileptr, "mode number map from partition to global : < %d > %d \n",
//              (info->nshg)*sizeof(int)+sizeof(char), info->nshg);
//      fwrite(info->ncvec, sizeof(int), info->nshg, fileptr);
//      free(info->ncvec);
//      fprintf(fileptr, "\n");
//      fclose(fileptr);
//    }
  
  if (SONFATH) {
// see writeNCorp -- function is diabled 
//     free(ifath);
//     free(nsons);
  }

  if (PMU_rank() == 0) { // clean up old files if they exist 
    int done = 0;
    for (i = PMU_size()+1; !done; i++) {
      sprintf(fname, "geombc.dat.%d", i);
      remove(fname);
      sprintf(fname, "restart.0.%d", i);
      done = remove(fname);
    }
    for (i = PMU_size()+1; !done; i++) {
      sprintf(fname, "geom.dat.%d", i);
      remove(fname);
      sprintf(fname, "bc.dat.%d", i);
      remove(fname);
      sprintf(fname, "restart.1.%d", i);
      done = remove(fname);
    }
  }
}

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "parallel.h"
#include "func.h"
#include "EnsaArrays.h"
#include "EnsaParameters.h"
#include <strings.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

#include "phastaIO.h"
/* Add other architectures as appropriate */
// fortran stuff taken out
//  #if defined(SGI) || defined(SUN4)
//    #define writeFortran wrunfmt_
//  #elif defined(IBM)
//    #define writeFortran wrunfmt
//  #endif



// fortran stuff taken out
//  extern "C"
//  void writeFortran(int*, int*, int*, int*, int*, int*, int*,
//                    int*, int*, int*, int*, int*, int*, int*,
//                    int*, int*, int*, int*, int*, int*, double*, int*, int*,
//                    double*, double*, double*, int*);
using namespace std;
// #ifndef SIM
// using namespace SCOREC_mesh;
// #endif
extern int getnumBshape(blockKey bKey);

// declared in procArgs (default "binary")
extern char* oformat;

extern int lstep;
extern int WRITEASC;
extern int phastaIO;
extern int old_format;
extern int FortFormFlag;
extern int zScale;
extern forwardblock Iblock;
extern forwardblock Bblock;
extern blockKey *RIblock;
extern blockKey *RBblock;
extern int* NBblock;
extern int* NIblock;
extern int* Nshape;
extern int* NshapeB;
extern int SONFATH;
extern char phversion[];
extern int totNumNodes;
extern int totNumEdges;
extern int totNumFaces;
extern int numTotParts;
extern int numParts;
extern int DisplacementMigration;
extern int dwalMigration;
extern int buildMapping;
extern int adaptFlag;
// extern char options[];


// globalInfo *info initialized  in initGlobalInfo
// EnsaParameters constructed in writeEnsaFiles()
// a single instance (on each proc) EnsaArrays constructed in writeEnsaFiles()
EnsaArrays::EnsaArrays(globalInfo *info, const EnsaParameters &par)
  : _info(info)
{
  int i, j, k, l, num, tmpblk;

  int nsd = par.getNSD();
  int numVars = par.getNUMVARS();
  int numEBC = par.getNUMEBC();
  int numNBC = par.getNUMNBC();

  int ndisp = 0;
  if(DisplacementMigration)
      ndisp = 3;
  int ndwal = 0;
  if(dwalMigration)
      ndwal = 1;
  int nmapping = 0;
  if(buildMapping && !adaptFlag)
      nmapping = 2;

  num = info->numnp;// number nodal points
  x = new double*[num];
  for (j=0; j < num; j++) x[j] = new double[nsd];//nsd=dim=3

  num = info->nshg;
  q = new double*[num];
  for (j=0; j < num; j++) q[j] = new double[numVars+ndisp+ndwal+nmapping];

  tmpblk = Iblock.size();
  ien = new int** [tmpblk];
  for(j=0; j < tmpblk; j++){
    num = NIblock[j];
    ien[j] = new int*[num];
    for (k=0; k < num; k++)
      ien[j][k] = new int[Nshape[j]];
  }

  tmpblk = Bblock.size();
  ienb = new int**[tmpblk];
  BCB = new double**[tmpblk]; /* saving calls to size and unnecessary */
  iBCB = new int**[tmpblk]; /* loops by allocating some stuff here */
  for(j=0; j<tmpblk; j++){
    num = NBblock[j];
    ienb[j] = new int*[num];
    BCB[j] = new double* [num];
    iBCB[j] = new int* [num];
    for (k=0; k < num; k++) {
      ienb[j][k] = new int[NshapeB[j]];
      BCB[j][k] = new double [numNBC];
      for(l=0;l<numNBC; l++) BCB[j][k][l] = 0.0;
      iBCB[j][k] = new int [2];
      iBCB[j][k][0] = iBCB[j][k][1] = 0;
    }
  }

  num = info->nshg;
  nBC = new int[num];

  num = info->nshg;
  iBC = new int[num];

#ifdef DEBUG
  printf("\n[%d] in EnsaArrays::EnsaArrays iBC size: %d \n", PMU_rank(),num);
#endif

  num = info->nshg;
  BC = new double*[num];
  for (j=0; j < num; j++) BC[j] = new double[numEBC];

  num = info->nshg;
  iper = new int[num];

  nlworkCalc(info);
  ilwork = new int[info->nlwork];

  for (i=0; i < info->nshg; i++){
    iBC[i] = 0;
    iper[i] = -1;
    nBC[i] = -1;
    for (j=0; j < numEBC; j++)
      BC[i][j] = 0.0;
  }
  for (i=0; i < info->nshg; i++)
    for (j=0; j < numVars+ndisp+ndwal+nmapping; j++)
      q[i][j] = 0.0;
}

void EnsaArrays::fixInd(globalInfo *info)
{
  int i, j, k, num, num2;

  int tmpblk = Iblock.size();
  for(j=0; j < tmpblk; j++){
    num = NIblock[j];
    num2 = Nshape[j];
    for (k=0; k < num; k++)
      for(i =0 ; i < num2; i++){
        if(ien[j][k][i] >=0) ien[j][k][i]++;
        else ien[j][k][i]--;
      }
  }

  tmpblk = Bblock.size();
  for(j=0; j < tmpblk; j++){
    num = NBblock[j];
    num2 = NshapeB[j];
    for (k=0; k < num; k++)
      for(i =0; i < num2; i++){
        if(ienb[j][k][i] >=0) ienb[j][k][i]++;
        else ienb[j][k][i]--;
      }
  }
  for (i=0; i < info->nshg; i++){
    iper[i]++;
    nBC[i]++;
  }
}

void EnsaArrays::write(pMesh mesh, globalInfo *info, const EnsaParameters &par,
                       double zscale[3], int *ifath, int *nsons, int ipart)
{
  int i, j, k, ind, *ienf, *iBCBf, *ienbf;
  char fname[30], proc[5];
  double *xf, *BCf, *BCBf, *qf, *qdisp, *qdwal;
  double *mapping_partid, *mapping_vtxid;
  FILE *fgmasc;
  FILE *fbcasc;
  FILE *frsasc;
  int tmpblk, tmpblkb, nelblk, bpoly, bnen, bnsh;
  int nElts, nShps;
  // convert from C to FORTRAN indexing
  fixInd(info);

  int numnp = info->numnp;
  int nshg = info->nshg;
  int nshgTot = info->nshgTot;
  int numel = info->numel;
  int numelb = info->numelb;
  int numpbc = info->numpbc;
  int nlwork = info->nlwork;
  int nshgOwn = info->nshgOwn;
  int nsd = par.getNSD();
  int nen = par.getNEN();
  int nfaces = par.getNFACES();
  int numflx = par.getNUMFLX();
  int numVars = par.getNUMVARS();
  int numEBC = par.getNUMEBC();
  int numNBC = par.getNUMNBC();
  int fType = par.getFTYPE();

  xf = new double[nsd*numnp];
  BCf = new double[numEBC*numpbc];
  qf = new double[numVars*nshg];
  if(DisplacementMigration)
      qdisp = new double[3*nshg];

  if(dwalMigration) 
      qdwal = new double[nshg];

  if(buildMapping && !adaptFlag) {
      mapping_partid = new double[nshg];
      mapping_vtxid = new double[nshg];
  }
  //////////////////////////////////////////////////////////////////////////////////////
  // FortFormFlag == 0
  //////////////////////////////////////////////////////////////////////////////////////
  if (!FortFormFlag)  {
    ind = 0;
    for (j=0; j < nsd; j++)
      for (i=0; i < numnp; i++)
        xf[ind++] = x[i][j];
    ind = 0;
    for (j=0; j < numEBC; j++)
      for (i=0; i < numpbc; i++)
        BCf[ind++] = BC[i][j];
    ind = 0;
    for(j=0;j<numVars;j++)
      for (i=0; i < nshg; i++)
        qf[ind++] = q[i][j];
    ind = 0;
    if(DisplacementMigration)
        for(j=0;j<3;j++)
            for(i=0;i<nshg;i++)
                qdisp[ind++] = q[i][j+numVars];
    ind = 0;
    if(dwalMigration) {
        for(j=0;j<1;j++) {
            for(i=0;i<nshg;i++) {
                int ndisp=0;
                if (DisplacementMigration) ndisp=3;
                qdwal[ind++] = q[i][j+ndisp+numVars];
            }
        }
    }
    ind = 0;
    if(buildMapping && !adaptFlag) {
//        for(j=0;j<2;j++) {
            for(i=0;i<nshg;i++) {
                int ndisp=0;
                if (DisplacementMigration) ndisp=3;
                int ndwal=0;
                if (dwalMigration) ndwal=1;
                mapping_partid[ind] = q[i][0+ndisp+ndwal+numVars];
                mapping_vtxid[ind]  = q[i][1+ndisp+ndwal+numVars];
                //printf("Writing mapping - rank: %d - vtx: %d - pID: %f - vID: %f\n",PMU_rank(),i+1,mapping_partid[ind]+1.0,mapping_vtxid[ind]+1.0);
                ind++;
            }
//       }
    } 
    ind = 0;
    if (zScale) {
      for (i=0; i < numnp; i++) {
        xf[i*3+0] *= zscale[0];
        xf[i*3+1] *= zscale[1];
        xf[i*3+2] *= zscale[2];
      }
    }
  }// if (!FortFormFlag)
  //////////////////////////////////////////////////////////////////////////////////////
  // FortFormFlag == 1 
  //////////////////////////////////////////////////////////////////////////////////////
  else  {
      if(PMU_rank() == 1){
          cout<<"\n FortFormFlag feature deprecated: exiting\n";
      }
      exit(0);

//      ind = 0;
//      for (i=0; i < numnp; i++)
//        for (j=0; j < nsd; j++)
//          xf[ind++] = x[i][j];
//      ind = 0;
//      for (i=0; i < numpbc; i++)
//        for (j=0; j < numEBC; j++)
//          BCf[ind++] = BC[i][j];
//      ind = 0;
//      for (i=0; i < nshg; i++)
//        for (j=0; j < numVars; j++)
//          qf[ind++] = q[i][j];
//      ind = 0;
//      if (zScale) {
//        for (i=0; i < numnp; i++) {
//          xf[i] *= zscale[0];
//          xf[i+numnp] *= zscale[1];
//          xf[i+2*numnp] *= zscale[2];
//        }
//      }//if (zScale)

//      /* Calc of nshape and nshapeb may be wrong.. 
//       * cf ~yaworm..3_4_0..writeEnsaFiles.cc  -- SST */
//      int fdof, nshape = -1, nshapeb = -1;
//      tmpblk = Iblock.size();
//      for (k=0; k < tmpblk; k++)  {
//        if (nshape < Nshape[k])
//          nshape = Nshape[k];
//      }
//      FIter mfiter = M_faceIter(mesh);
//      pFace face;
//      while ( (face = FIter_next(mfiter)) )  {
//        if (F_whatInType(face) != Gregion)  {   /* bdry face */
//          fdof = F_numEdges(face);
//          if (fdof  > nshapeb)  /* only works for linear */
//            nshapeb =  fdof;
//        }
//      } 

//      ienf  = new int[nshape*numel];
//      ienbf = new int[nshape*numelb];
//      iBCBf = new int[numelb*2];
//      BCBf  = new double[numelb*nshapeb*numNBC];

//      ind=0;
//      tmpblk = Iblock.size(); /* Number of topologies */
//      for (k=0; k < tmpblk; k++) {
//        nElts = NIblock[k]; /* Number of elements for topology k */
//        nShps = Nshape[k]; /* Number of modes per element */
//        for (i=0; i < nElts; i++) {
//          for (j=0; j < nShps; j++)
//            ienf[ind++] = ien[k][i][j];
//          for (j=nShps; j < nshape; j++)
//            ienf[ind++] = ien[k][i][nShps-1];
//        }
//      }
//      ind=0;
//      tmpblk = Bblock.size();
//      for (k=0; k < tmpblk; k++) {
//        nElts = NBblock[k];
//        nShps = NshapeB[k];
//        /* NshapeB is the number of modes in the entire boundary element (not
//         * just on the boundary face) */
//        for (i=0; i < nElts; i++) {
//          for (j=0; j < nShps; j++)
//            ienbf[ind++] = ienb[k][i][j];
//          for (j=nShps; j < nshape; j++)
//            ienbf[ind++] = ienb[k][i][nShps-1];
//        }
//      }
//      ind = 0;
//      for (k=0; k < tmpblk; k++) {
//        nElts = NBblock[k];
//        for (i=0; i < nElts; i++)
//          for (j=0; j < 2; j++)
//            iBCBf[ind++] = iBCB[k][i][j];
//      }
//      ind = 0;
//      for (k=0; k < tmpblk; k++) {
//        nElts = NBblock[k];
//        nShps = RBblock[k].nenbl;
//        for (i=0; i < nElts; i++) {
//          for (j=0; j < numNBC*nShps; j++)
//            BCBf[ind++] = BCB[k][i][j];
//          for (j=numNBC*nShps; j < numNBC*nshapeb; j++)
//            BCBf[ind++] = BCB[k][i][numNBC*nShps-1];
//        }
//      }
//      writeFortran(&numVars, &numnp, &nshg, &nsd, &numel, &numelb, &nen,
//              &nfaces, &numflx, &numpbc, &nlwork, &nshape, &nshapeb, &PMU_rank(),
//              &nshgOwn, ilwork, nBC, iBC, iper, iBCBf, xf, ienf, ienbf,
//              BCf, qf, BCBf, &fType);
//      delete[] ienf;
//      delete[] ienbf;
//      delete[] iBCBf;
//      delete[] BCBf;
  }// else FortFormFlag == 1

  ////////////////////////////////////////////////////////////////////////////////////////////
  /* write ascii files for debugging */
  ////////////////////////////////////////////////////////////////////////////////////////////
  if(WRITEASC) {
    sprintf(proc, "%d", PMU_rank()+1);
    strcpy(fname, "geom.asc.");
    strcat(fname, proc);
    fgmasc = fopen(fname, "w");

    strcpy(fname, "bc.asc.");
    strcat(fname, proc);
    fbcasc = fopen(fname, "w");

    strcpy(fname, "restart.asc.");
    strcat(fname, proc);
    frsasc = fopen(fname, "w");
    tmpblk = Iblock.size(); /* nbr of blocks */

    fprintf(fgmasc, "numnp, nshg, nsd, numel, numelb, nen, nfaces, numflx, nlwork\n");
    fprintf(fgmasc, "%d %d %d %d %d %d %d %d %d %u %d\n",
            numnp, nshg, nsd, numel, numelb, nen,
            nfaces, numflx, tmpblk, Bblock.size(), nlwork);
    for (i=0; i < nlwork; i++)
      fprintf(fgmasc, "%d\n", ilwork[i]);

    fprintf(fgmasc, "\nNodal coordinates:\n");
    fprintf(fgmasc, "------------------\n");
    for(i=0; i < numnp; i++)
      fprintf(fgmasc, "%f %f %f\n", x[i][0],
              x[i][1], x[i][2]);

    int nblock = 0;
    int nbblock = 0;
    fprintf(fgmasc, "\nInterior element-nodal connectivity:\n");
    fprintf(fgmasc, "------------------------------------\n");
    for(k=0; k < tmpblk; k++){
      fprintf(fgmasc, "%d %d %d %d %d %d %d\n", NIblock[k],
              RIblock[k].nen, RIblock[k].maxpoly, Nshape[k],
              getnumBshape(RIblock[k]), RIblock[k].nenbl,
              RIblock[k].lcsyst);
      fprintf(fgmasc, "********************\n");
      for(i=0; i < NIblock[k]; i++){
        for(j=0; j < Nshape[k]; j++)
          fprintf(fgmasc, "%d ", ien[k][i][j]);
        fprintf(fgmasc, "\n");
      }
      fprintf(fgmasc, "\n ---------------------------------------------- \n");
      nblock++;
    }

    fprintf(fgmasc, "\nBoundary element-nodal connectivity:\n");
    fprintf(fgmasc, "------------------------------------\n");
    tmpblk = Bblock.size();
    for(k=0; k < tmpblk; k++){
      fprintf(fgmasc, "%d %d %d %d %d %d %d\n", NBblock[k],
              RBblock[k].nen, RBblock[k].maxpoly, NshapeB[k],
              getnumBshape(RBblock[k]), RBblock[k].nenbl,
              RBblock[k].lcsyst);
      fprintf(fgmasc, "********************\n");
      for(i=0; i < NBblock[k]; i++){
        for(j=0; j < NshapeB[k]; j++)
          fprintf(fgmasc, "%d ", ienb[k][i][j]);
        fprintf(fgmasc, "\n");
      }
      fprintf(fgmasc, "\n ---------------------------------------------- \n");
      nbblock++;
    }

    fprintf(fbcasc, "\nBoundary condition mapping array (nBC):\n");
    fprintf(fbcasc, "---------------------------------------\n");
    fprintf(fbcasc, "%d\n", numpbc);
    for(i=0; i < nshg; i++)
      fprintf(fbcasc, "%d\n", nBC[i]);

    fprintf(fbcasc, "\nBoundary condition codes array (iBC):\n");
    fprintf(fbcasc, "-------------------------------------\n");
    for(i=0; i < numpbc; i++)
      fprintf(fbcasc, "%d\n", iBC[i]);

    fprintf(fbcasc, "\nBoundary condition array (BC):\n");
    fprintf(fbcasc, "------------------------------\n");
    for(i=0; i < numpbc; i++){
      for(j=0; j < numEBC; j++)
        fprintf(fbcasc, "%f ", BC[i][j]);
      fprintf(fbcasc, "\n");
    }

    fprintf(fbcasc, "\nPeriodic masters array (iper):\n");
    fprintf(fbcasc, "------------------------------\n");
    for(i=0; i < nshg; i++)
      fprintf(fbcasc, "%d\n", iper[i]);

    fprintf(fbcasc, "\niBCB and BCB are written Staggered for each block in unf:\n");
    fprintf(fbcasc, "\nNatural boundary condition codes array (iBCB):\n");
    fprintf(fbcasc, "----------------------------------------------\n");
    for(k=0; k < tmpblk; k++){
      fprintf(fbcasc, "%d %d %d %d\n", NBblock[k], RBblock[k].nen,
              RBblock[k].maxpoly, NshapeB[k]);
      fprintf(fbcasc, "********************\n");
      for(i=0; i < NBblock[k]; i++)
        for(j=0; j< 2 ; j++)
          fprintf(fbcasc, "%d\n", iBCB[k][i][j]);
    }
    fprintf(fbcasc, "\nNatural boundary condition array (BCB):\n");
    fprintf(fbcasc, "---------------------------------------\n");
    for(k=0; k < tmpblk; k++){
      fprintf(fbcasc, "%d %d %d %d\n", NBblock[k], RBblock[k].nen,
              RBblock[k].maxpoly, NshapeB[k]);
      fprintf(fbcasc, "********************\n");
      for(i=0; i < NBblock[k]; i++)
        for(j=0;j<numNBC; j++)
          fprintf(fbcasc, "%f\n", BCB[k][i][j]);
    }

    for (i=0; i < nshg; i++){
      for (j=0; j < numVars; j++)
        fprintf(frsasc, "%f ", q[i][j]);
      fprintf(frsasc, "\n");
    }

    fclose(fbcasc);
    fclose(fgmasc);
    fclose(frsasc);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Now the Unformatted stuff OR phastaIO
  // if above FortFormFlag==1 performs necessary steps (besides fixInd)  towards
  // phastaIO
  ////////////////////////////////////////////////////////////////////////////////////////////
  if(phastaIO==1){
                  
      char fname[255];
      
      int magic_number = 362436;
      int* mptr = &magic_number;
      int fgeom, frest;
      int iarray[10];
      int size, nitems;
      char dirname[25];
      char filename[255];
      int isize;
      int maxshg;
      
      int numdir;
      int myrank;
      int mysize;

      numdir = 2048;
      myrank = PMU_rank();
      mysize = PMU_size();
      // dir setup for parallel case
      //sprintf(dirname,"%d-procs_case/",numTotParts);
      //sprintf(dirname,"%d/",PMU_rank());
     
      if (ipart == 0) { // We create the directory(ies) that are required
        if (myrank == 0 ) {
          sprintf(dirname,"%d-procs_case/",numTotParts);
          //printf("Creating a directory -  rank, ipart, dirname: %d %d %s\n",myrank,ipart,dirname);       
          if ( mkdir( dirname, 00755 ) ) {
            cerr << "overwriting directory " << dirname << " - Rank " << myrank << "- ipart " << ipart <<  endl;
          }
 
          sprintf(filename,"%snumpe.in",dirname);
          ofstream npe(filename);
          npe << numTotParts << endl ;
          npe.close();
      
          sprintf(filename,"%snumstart.dat",dirname);
          ofstream numstart(filename);
          numstart<< lstep << endl ;
          numstart.close();
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (mysize > numdir) { // if mysize > numdir, then we will write the phasta files in sub dir. Create them now
          sprintf(dirname,"%d-procs_case/%d/",numTotParts,myrank%numdir);
          // what for :unlink( dirname );
          if(myrank < numdir){
            //printf("Creating a directory -  rank, ipart, dirname: %d %d %s\n",myrank,ipart,dirname);       
            if ( mkdir( dirname, 00755 ) ) {
              //cerr << "overwriting directory " << dirname << endl;
              cerr << "overwriting directory " << dirname << " - Rank " << myrank << "- ipart " << ipart <<  endl;
            }
          }
          MPI_Barrier(MPI_COMM_WORLD);
        } // myrank < numdir

      } // ipart == 0
 
      if(mysize > numdir) { 
          sprintf(dirname,"%d-procs_case/%d/",numTotParts,myrank%numdir);
      }
      else { // size is < numdir so we write directly in %d-proc_case
          sprintf(dirname,"%d-procs_case/",numTotParts);
      }
     
      //printf("EnsaArrays: myrank %d - ipart %d - dirname %s\n",myrank,ipart,dirname);
 
      bzero( (void*)filename,255);
      sprintf(filename,"%sgeombc.dat.%d",dirname,PMU_rank()*numParts+ipart+1 );
      openfile_( filename, "write", &fgeom );
      
      bzero( (void*)filename, 255 );
      sprintf(filename,"%srestart.%d.%d",dirname, lstep, PMU_rank()*numParts+ipart+1 );
      openfile_( filename, "write", &frest );
      
      
      
      tmpblk = Iblock.size(); // number of blocks //
      tmpblkb = Bblock.size(); // num boundary blocks //
      
      // before anything we put in the standard headers //
      writestring_( &fgeom,"# PHASTA Input File Version 2.0\n");
      writestring_( &frest,"# PHASTA Input File Version 2.0\n");
      
      writestring_( &fgeom, "# Byte Order Magic Number : 362436 \n");
      writestring_( &frest, "# Byte Order Magic Number : 362436 \n");
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"# Output generated by NSpre version: %s \n", phversion);
      writestring_( &fgeom, fname );
      writestring_( &frest, fname );
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"# CmdLineArgs : \n");
      writestring_( &fgeom, fname );
      writestring_( &frest, fname );
      
      time_t timenow = time ( &timenow);
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"# %s\n", ctime( &timenow ));
      writestring_( &fgeom, fname );
      writestring_( &frest, fname );
      
      int one=1;
      
      size = 1;
      nitems = 1;// length of iarray
      iarray[ 0 ] = 1;
      writeheader_( &fgeom, "byteorder magic number ",
                    (void*)iarray, &nitems, &size, "integer", oformat );
      
      writedatablock_( &fgeom, "byteorder magic number ",
                       (void*)mptr, &nitems, "integer", oformat );
      
      writeheader_( &frest, "byteorder magic number ",
                    (void*)iarray, &nitems, &size, "integer", oformat );
      
      writedatablock_( &frest, "byteorder magic number ",
                       (void*)mptr, &nitems, "integer", oformat );
      
      
      ///////////////////////////////////////////////////////////////////////////////      
      // writing the restart //
      ///////////////////////////////////////////////////////////////////////////////
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of modes : < 0 > %d\n", nshg);
      writestring_( &frest, fname );
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of variables : < 0 > %d\n", numVars);
      writestring_( &frest, fname );
      
      size =  numVars*nshg;
      nitems = 3;// length of iarray
      iarray[ 0 ] = nshg;
      iarray[ 1 ] = numVars;
      iarray[ 2 ] = lstep;
      
      writeheader_( &frest, "solution ",
                    ( void* )iarray, &nitems, &size,"double", oformat );
      
      nitems = numVars*nshg;
      writedatablock_( &frest, "solution ",
                       ( void* )(qf), &nitems, "double", oformat );
     if(DisplacementMigration){
          size =  3*nshg;
          nitems = 3;// length of iarray
          iarray[ 0 ] = nshg;
          iarray[ 1 ] = 3;
          iarray[ 2 ] = lstep;
          
          writeheader_( &frest, "displacement ",
                        ( void* )iarray, &nitems, &size,"double", oformat );
          
          nitems = 3*nshg;
          writedatablock_( &frest, "displacement ",
                       ( void* )(qdisp), &nitems, "double", oformat );     
      }
      if(dwalMigration){
          size = nshg;
          nitems = 3;// length of iarray
          iarray[ 0 ] = nshg;
          iarray[ 1 ] = 1;
          iarray[ 2 ] = lstep;
          
          writeheader_( &frest, "dwal",
                        ( void* )iarray, &nitems, &size,"double", oformat );
          
          nitems = nshg;
          writedatablock_( &frest, "dwal",
                       ( void* )(qdwal), &nitems, "double", oformat );     
      }
      if(buildMapping && !adaptFlag){
          size = nshg;
          nitems = 3;// length of iarray
          iarray[ 0 ] = nshg;
          iarray[ 1 ] = 1;
          iarray[ 2 ] = lstep;
          
          writeheader_( &frest, "mapping_partid",
                        ( void* )iarray, &nitems, &size,"double", oformat );
          //writeheader_( &frest, "mapping_partid",
          //              ( void* )iarray, &nitems, &size,"integer", oformat );
          
          nitems = nshg;
          writedatablock_( &frest, "mapping_partid",
                       ( void* )(mapping_partid), &nitems, "double", oformat );    
          //writedatablock_( &frest, "mapping_partid",
          //             ( void* )(mapping_partid), &nitems, "integer", oformat );    

          size = nshg;
          nitems = 3;// length of iarray
          iarray[ 0 ] = nshg;
          iarray[ 1 ] = 1;
          iarray[ 2 ] = lstep;
          
          writeheader_( &frest, "mapping_vtxid",
                        ( void* )iarray, &nitems, &size,"double", oformat );
          //writeheader_( &frest, "mapping_vtxid",
          //              ( void* )iarray, &nitems, &size,"integer", oformat );
          
          nitems = nshg;
          writedatablock_( &frest, "mapping_vtxid",
                       ( void* )(mapping_vtxid), &nitems, "double", oformat );     
          //writedatablock_( &frest, "mapping_vtxid",
          //             ( void* )(mapping_vtxid), &nitems, "integer", oformat );     
 
 
      }
      
      closefile_( &frest, "write" );
      // finished writing the restart //
      
      ///////////////////////////////////////////////////////////////////////////////
      // first write the headers with global data
      ///////////////////////////////////////////////////////////////////////////////
      
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of nodes : < 0 > %d \n", numnp);
      writestring_( &fgeom, fname );//has old equiv
      
      
      
      // required?
      int localNumNodes;
      int numNodes = totNumNodes;
      pVertex v;
      int count = 0;
//       VIter vIter = M_vertexIter( mesh );
//       while ( v = VIter_next( vIter ) ) {
//           if( EN_isOwnerPart  (( pEntity )v , 0))//only one part per proc   
//               count++;
//       }
//       VIter_delete(vIter);

      
// #ifdef PARALLEL
//       MPI_Allreduce(&count, &numNodes, 1, MPI_INT, 
//                     MPI_SUM, MPI_COMM_WORLD);
// #endif
    
      bzero( (void*)fname, 255 );
      sprintf( fname,"number of nodes in the mesh : < 0 > %d \n",
               numNodes);// has old equiv
      writestring_( &fgeom, fname );
      
      
      // required?
      count=0;
      int localNumEdges;
      int numEdges = totNumEdges;
//       pEdge e;
//       EIter eIter = M_edgeIter( mesh );
//       while ( e = EIter_next( eIter ) ) {
//           if( EN_isOwnerPart  (( pEntity )e , 0))//only one part per proc   
//               count++;
//       }
//       EIter_delete(eIter);

// #ifdef PARALLEL
//       MPI_Allreduce(&count, &numEdges, 1, MPI_INT, 
//                     MPI_SUM, MPI_COMM_WORLD);  
// #endif
      
      bzero( (void*)fname, 255 );
      sprintf( fname,"number of edges in the mesh : < 0 > %d \n",
               numEdges );// NO old equiv !
      writestring_( &fgeom, fname );

      // required? YES, at least some of the stuff
      count=0;
      int localNumFaces;
      int numFaces=totNumFaces;
//       pFace f;
//       FIter fIter = M_faceIter( mesh );
//       while ( f = FIter_next( fIter ) ) {
//           if( EN_isOwnerPart  (( pEntity )f , 0))only one part per proc   
//               count++;
//       }    
//       FIter_delete(fIter);

// #ifdef PARALLEL
//       MPI_Allreduce(&count, &numFaces, 1, MPI_INT, 
//                     MPI_SUM, MPI_COMM_WORLD);    
// #endif
      
      
      bzero( (void*)fname, 255 );
      sprintf( fname,"number of faces in the mesh : < 0 > %d \n",
               numFaces);// NO old equiv !
      writestring_( &fgeom, fname );
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of modes : < 0 > %d \n", nshg);
      writestring_( &fgeom, fname );// has old equiv
      
      bzero( (void*)fname, 255 );// mind the typo [sic!]
      sprintf(fname,"number of shapefunctions soved on processor : < 0 > %d \n", info->nshgOwn);
      writestring_( &fgeom, fname );// NOT in that property --> see below
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of global modes : < 0 > %d \n", nshgTot);
      writestring_( &fgeom, fname );// has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of interior elements : < 0 > %d \n", numel);
      writestring_( &fgeom, fname );// has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of boundary elements : < 0 > %d \n", numelb);
      writestring_( &fgeom, fname );// has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"maximum number of element nodes  : < 0 > %d \n", nen);
      writestring_( &fgeom, fname );// has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of interior tpblocks : < 0 > %d \n", tmpblk); 
      writestring_( &fgeom, fname ); // has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of boundary tpblocks : < 0 > %d \n", tmpblkb);
      writestring_( &fgeom, fname ); // has old equiv
      
      bzero( (void*)fname, 255 );
      sprintf(fname,"number of nodes with Dirichlet BCs : < 0 > %d \n", numpbc);
      writestring_( &fgeom, fname );// has old equiv

      printf("MeshInfo: part %d - nodes %d\n",PMU_rank()*numParts+ipart+1, numnp);
      printf("MeshInfo: part %d - elms %d\n",PMU_rank()*numParts+ipart+1, numel);
      
      ////////////////////////////////////////////////////////////////////////////////
      // coordinates
      ////////////////////////////////////////////////////////////////////////////////
      size = 3*numnp;
      nitems = 2 ;// length of iarray
      iarray[ 0 ] = numnp;
      iarray[ 1 ] = 3;
      writeheader_( &fgeom, "co-ordinates ", (void*)iarray, &nitems, &size,
                    "double", oformat  );
      
      nitems = 3*numnp;// xf <==> e->x ???
      writedatablock_( &fgeom, "co-ordinates ", (void*)(xf), &nitems,
                       "double", oformat );
      

      // parallel stuff: below code taken from 3.6.6 PartitionDriver
      // used data to be modified by the local data structures
      ////////////////////////////////////////////////////////////////////////
      // no new equivalent for these  ==> check parallel --> yes
      
      bzero((void*)filename, 255 );
      sprintf(filename,"number of processors : < 0 > %d \n", numTotParts);
      writestring_( &fgeom, filename );
      
      ////////////////////////////////////////////////////////////////////////  
      
       
      ////////////////////////////////////////////////////////////////////////
      // no new equivalent for these   ==> check parallel   --> yes
      bzero( (void*)filename, 255 );// nlwork <==> ilwork.size()
      sprintf(filename,"size of ilwork array : < 0 > %d\n", nlwork);
      writestring_( &fgeom, filename );
      
      ////////////////////////////////////////////////////////////////////////         
      
      
      ////////////////////////////////////////////////////////////////////////
      // no new equivalent for these   ==> check parallel --> yes
      // 
      
      // in contrast to NSpre2 this is used under a conditional
      if(nlwork > 0){
          isize = nlwork;
          nitems = 1;
          iarray[ 0 ] =  nlwork;
          writeheader_( &fgeom, "ilwork ", (void*)iarray,
                        &nitems, &isize, "integer", oformat );
          
          
          nitems = nlwork;
          writedatablock_( &fgeom, "ilwork ", (void*)(ilwork), &nitems,
                           "integer", oformat );
          
      }
      ///////////////////////////////////////////////////////////////////////   
      
      
      
      
      ///////////////////////////////////////////////////////////////////////     
      // implemented in 3.6.6 PartitionDriver --  but here hidden in
      // writeEnsaFiles 
      // ncvec allocated and filled in writeNCorp(..)
      // ncvec for master is size maxshg*nproc with nproc=PMU_size()
      // ncvec for others are size maxshg 
      // this is required because previously this header was written elsewhere
#ifdef PARALLEL
      
      /* Determine the largest (MPI_MAX) nshg among all processors */
      MPI_Allreduce(&info->nshg, &maxshg, 1, MPI_INT, 
                    MPI_MAX, MPI_COMM_WORLD);
#else
      maxshg = info->nshg;
#endif
      // (a loops from 0 PMU_size()-1)  ncorp[a].size() ;
      //isize = PMU_rank()?maxshg:PMU_size()*maxshg;
      
      isize=nshg;
      nitems = 1;// length of iarray
      //iarray[ 0 ] = PMU_rank()?maxshg:PMU_size()*maxshg;
      iarray[ 0 ] = nshg;
      writeheader_( &fgeom, " mode number map from partition to global",
                    (void*)iarray, &nitems, &isize, "integer", oformat );
      
      
//          int* fncorp = new int [  ncorp[a].size() ];
//          int count = 0;
//          for( map< int, int>::iterator iter = ncorp[a].begin();
//               iter != ncorp[a].end();
//               iter++ )
//              fncorp[count++] = (*iter).second ;
      
      // ncvec initialized in writeNCorp.c
      // nitems = PMU_rank() ?maxshg:PMU_size()*maxshg; // lengh of each proc's
      //ncvec
      nitems=nshg;
      writedatablock_( &fgeom, " mode number map from partition to global",
                       (void*)(info->ncvec), &nitems, "integer", oformat );
      ///////////////////////////////////////////////////////////////////////////////
      
      
      
      ///////////////////////////////////////////////////////////////////////////////
      // interior connectivity
      // now we enter the loop over interior blocks where we write ien //
      ///////////////////////////////////////////////////////////////////////////////
      char keyphrase[100]; // to store the dynamic string for each tpblock //
      int bnenbl, bnshlb, blcsyst;
      for(i=0; i< tmpblk; i++) {     // for interior each block //
          bnen = RIblock[i].nen;      // nen of the block -- topology //
          bpoly = RIblock[i].maxpoly; // polynomial order of the block //
          nelblk = NIblock[i] ;       // numel of this block //
          bnsh  = Nshape[i] ;         // nshape of this block // 
          bnshlb = getnumBshape(RIblock[i]);
          bnenbl = RIblock[i].nenbl;
          blcsyst = RIblock[i].lcsyst;
          
          // generate the key phrase specific to this block //
          
          generate_keyphraseNew( keyphrase,"connectivity interior ", &RIblock[i]);
          
          size = nelblk*bnsh;
          nitems = 7;
          iarray[ 0 ] = nelblk;
          iarray[ 1 ] = bnen;
          iarray[ 2 ] = bpoly;
          iarray[ 3 ] = bnsh;
          iarray[ 4 ] = bnshlb;
          iarray[ 5 ] = bnenbl;
          iarray[ 6 ] = blcsyst;
          
          
          ienf = new int [nelblk*bnsh];
          int ind = 0;
          for(j =0; j< bnsh; j++)
              for(k=0; k< nelblk; k++)
                  ienf[ind++] = ien[i][k][j];
          
          
          writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                        "integer", oformat );
          
          nitems = nelblk*bnsh;//ienf <==>  e->ien_sms[i] ???
          writedatablock_( &fgeom, keyphrase, (void*)(ienf), &nitems,
                           "integer", oformat );
          
          
          delete [] ienf;
// not used in old pNSpre
//          generate_keyphraseNew( keyphrase,"ien to sms ", &RIblock[i]);
          
//          size = nelblk;
//          nitems = 1;
//          writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
//                        "integer", oformat );

//          nitems = nelblk;// ienf <==>  e->ien_sms[i] ???
//          writedatablock_( &fgeom, keyphrase, (void*)(ienf), &nitems,
//                           "integer", oformat );
      }
      
      ///////////////////////////////////////////////////////////////////////////////
      // boundary connectivity
      // similar to the above we also write the boundary blocks //
      // along with these we also write the Natural boundary conditions //
      ///////////////////////////////////////////////////////////////////////////////
      for(i=0; i< tmpblkb; i++){
          
          bnen = RBblock[i].nen;
          bpoly = RBblock[i].maxpoly;
          bnenbl = RBblock[i].nenbl;
          blcsyst = RBblock[i].lcsyst;
          nelblk = NBblock[i]; // num of elements in this block //
          bnsh = NshapeB[i]; // num of shapefuncs for each element //
          bnshlb = getnumBshape(RBblock[i]);
          
          
          ienf = new int [nelblk*bnsh];
          iBCBf = new int [nelblk*2];
          BCBf = new double [nelblk*numNBC];
          
          for(k=0; k< nelblk; k++){
              for(j=0; j< bnsh; j++) ienf[k+j*nelblk]=ienb[i][k][j];
              for(j=0;j<2;j++) iBCBf[k+j*nelblk] = iBCB[i][k][j];
              for(j=0;j<numNBC;j++) BCBf[k+j*nelblk] =BCB[i][k][j];
          }
          
          iarray[ 0 ] = nelblk;
          iarray[ 1 ] = bnen;
          iarray[ 2 ] = bpoly;
          iarray[ 3 ] = bnsh;
          iarray[ 4 ] = bnshlb;
          iarray[ 5 ] = bnenbl;
          iarray[ 6 ] = blcsyst;
          iarray[ 7 ] = numNBC;
          
          // generate the key phrase specific to this block //
          generate_keyphraseNew( keyphrase,"connectivity boundary ",
                                 &RBblock[i]);
          size = nelblk*bnsh;

          if (old_format) size = nelblk*bnsh + nelblk*2 + nelblk*numNBC;
          nitems = 8;
          writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                        "integer", oformat );
          
          nitems = nelblk*bnsh;//ienf <==> e->ienb[i]
          writedatablock_( &fgeom, keyphrase, (void*)(ienf), &nitems,
                           "integer", oformat );
          
          if (!old_format){
              
              generate_keyphraseNew(keyphrase,"nbc codes ", &RBblock[i]);
              size = nelblk*2;
              nitems = 8;
              
              writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                            "integer", oformat );
          }
          
          nitems = nelblk*2;
          writedatablock_( &fgeom, keyphrase, (void*)(iBCBf), &nitems,
                           "integer", oformat );
          
          if (!old_format){
              
              generate_keyphraseNew(keyphrase,"nbc values ", &RBblock[i]);
              size = nelblk*numNBC;
              nitems = 8;
              
              writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                            "double", oformat );
              
          }
          nitems = nelblk*numNBC;
          writedatablock_( &fgeom, keyphrase, (void*)(BCBf), &nitems,
                           "double", oformat );
          
          delete[] ienf;
          delete [] iBCBf;
          delete [] BCBf;
          
      }//for(i=0; i< tmpblkb; 
      
      ///////////////////////////////////////////////////////////////////////////////
      //    the unblocked Essential Boundary conditions //
      ///////////////////////////////////////////////////////////////////////////////
      size = nshg;
      nitems = 1;
      iarray[ 0 ] = nshg;
      writeheader_( &fgeom , "bc mapping array ", (void *)iarray, &nitems,
                    &size, "integer", oformat );
      
      nitems = nshg; //nBC <==> e->nBC
      writedatablock_( &fgeom, "bc mapping array ", (void*)nBC, &nitems ,
                       "integer", oformat );
      
      size = numpbc;
      nitems = 1;
      iarray[ 0 ] = numpbc;
      writeheader_( &fgeom , "bc codes array ", (void *)iarray, &nitems,
                    &size, "integer", oformat );

      nitems = numpbc;
      writedatablock_( &fgeom, "bc codes array ", (void*)iBC, &nitems ,
                       "integer", oformat );
      
      size = numpbc*(numVars+7);
      nitems = 1;
      iarray[ 0 ] = numpbc*(numVars+7);
      writeheader_( &fgeom , "boundary condition array ", (void *)iarray, &nitems,
                    &size, "double", oformat );

    nitems = numpbc*(numVars+7);
    writedatablock_( &fgeom, "boundary condition array ", (void*)(BCf),
                     &nitems , "double", oformat );
    
    // this is taken up again in PartitionDriver
    size = nshg;
    nitems = 1;
    iarray[ 0 ] = nshg;
    writeheader_( &fgeom , "periodic masters array ", (void *)iarray, &nitems,
                  &size, "integer", oformat );
    
    nitems = nshg;
    writedatablock_( &fgeom, "periodic masters array ", (void*)iper,
                     &nitems , "integer", oformat );
    
    // parallel stuff: taken from PartitionDriver.cc
    if ( SONFATH > 0 ) {
        cout<<"\nERROR: SONFATH not implemented at the moment\n";
        exit(1);
    }
    
    closefile_( &fgeom, "write" );
    
    
    
  }//phastaIO==1 
  
  
  delete[] Nshape;
  delete[] NshapeB;
  delete[] RIblock;
  delete[] RBblock;

  delete[] xf;
  delete[] BCf;
  delete[] qf;
  if(DisplacementMigration)
      delete [] qdisp;
  if(dwalMigration)
      delete [] qdwal;
  if(buildMapping && !adaptFlag) {
      delete [] mapping_partid;
      delete [] mapping_vtxid;
  }
}//void EnsaArrays::write


EnsaArrays::~EnsaArrays()
{
  int j, k, num, tmpblk;

  num = _info->numnp;
  for (j=0; j < num; j++) delete[] x[j];
  delete[] x;

  num = _info->nshg;
  for (j=0; j < num; j++) delete[] q[j], delete[] BC[j];
  delete[] q, delete[] BC;

  tmpblk = Iblock.size();
  for (j=0; j < tmpblk; j++) {
    num = NIblock[j];
    for (k=0; k < num; k++) delete[] ien[j][k];
    delete[] ien[j];
  }
  delete[] ien;

  tmpblk = Bblock.size();
  for (j=0; j < tmpblk; j++) {
    num = NBblock[j];
    for (k=0; k < num; k++)
      delete[] ienb[j][k], delete[] BCB[j][k], delete[] iBCB[j][k];
    delete[] ienb[j], delete[] BCB[j], delete[] iBCB[j];
  }
  delete[] ienb, delete[] BCB, delete[] iBCB;
  delete[] nBC, delete[] iBC, delete[] iper;
  delete[] NIblock, delete[] NBblock;
  delete[] ilwork;

  //clean up Iblock, this part is critical for multiple parts on each proc
  forwardblock::iterator iblockiter, bblockiter;
  Iblock.clear();
  Bblock.clear();
}

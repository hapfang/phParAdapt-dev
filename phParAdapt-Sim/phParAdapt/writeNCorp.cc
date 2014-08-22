/**********************************************************************/
/* write the array that gives the global vertex number for each local */
/* number                                                             */
/**********************************************************************/
#include "parallel.h"
#include "func.h"
#include <assert.h>
#ifdef SIM
#include "MeshSim.h"
#endif

// fortran stuff taken out Feb. 2004-JM
/*  #ifdef SUN4 */
/*  #define  SonFath sonfath_ */
/*  #elif sgi */
/*  #define  SonFath sonfath_ */
/*  #elif ibm6000 */
/*  #define  SonFath sonfath */
/*  #endif */

#ifdef __cplusplus
extern "C" {
#endif

/* global variables */
//extern void SonFath(int*,int*,int*,int*,int*,int*, int*, int*);

extern int per,ICtyp;
extern int SONFATH, FortFormFlag;
extern pMeshDataId NDOFID;


//////////////////////////////////////////////////////////////////////////////////////////////
// this function basically fills  ncvec:
// the mapping of 
// ncvec(on-proc vertex number) = GlobalDOFStartingAt
// called in writeEnsaFiles.cc
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PARALLEL
void writeNCorp(pParMesh pmesh, globalInfo *info, int **nsons, int **ifath, int ipart)
#else
void writeNCorp(pMesh mesh, globalInfo *info, int **nsons, int **ifath)
#endif
{
  int sanity;
  int maxshg;
/*    int *colmaj; */
  int  *vnumnp;
  pVertex v; pEdge e; pFace f; pRegion r;
  VIter vIter;
  void *tmp;
  int i, j, nd, st, count = 0;
  FILE *fp;
#ifdef PARALLEL
  pMesh mesh = PM_mesh(pmesh, ipart);
  /* Determine the largest (MPI_MAX) nshg among all processors */
  // This Allreduce may need to be taken care of for multi-parts per proc. Min Zhou
  MPI_Allreduce(&info->nshg, &maxshg, 1, MPI_INT, 
                MPI_MAX, MPI_COMM_WORLD);
#else
  maxshg = info->nshg;
#endif

  /* ncvec for master is size maxshg*nproc */
  /* ncvec for others are size maxshg */
  info->ncvec =
    new int [maxshg];
  /* ncvec(on-proc vertex number) = GlobalDOFStartingAt1 */
  /* (DOF=Degree of Freedom) */

  vIter = M_vertexIter(mesh);
  while (v = VIter_next(vIter)) { 
/*      //  was set in mdb2phasta (for adaptFlag=0) */
    if(EN_getDataInt(v,info->incorp,&sanity))
      info->ncvec[count]=sanity+1;
    else
      info->ncvec[count]=EN_id((pEntity)v);
    count++;
/*      // alternative :was set in localInfo */
/*      info->ncvec[count++] = EN_dataI(v, "EQST") + 1; */   // obsolete function

#ifdef DEBUG
//    printf("\n[%d]in writeNCorp: retrieving data attached to a vertex: %d: \n", PMU_rank(), sanity+1);
#endif

  }
  
  VIter_delete(vIter);
  //    info->ncvec[count++] = EN_dataI(v, "smsNum") + 1;   // obsolete function

  /* stack edge DOFs onto ncvec next*/
  /* for each edge, loop over DOFs and insert GlobalDOF */ 

  if (info->edgeson)  {
    EIter eIter = M_edgeIter(mesh);
    while(e = EIter_next(eIter)) {
      if (EN_getDataInt(e, NDOFID, &nd) )  {   
        EN_getDataInt(e,info->incorp,&st);
        st++;

        for (i=0; i < nd; i++)
          info->ncvec[count++] = st + i;
      }
    }
    EIter_delete(eIter);
  }

  /* stack face DOFs onto ncvec next*/
  if (info->faceson)  {
    FIter fIter = M_faceIter(mesh);
    while (f = FIter_next(fIter)) {
      if (EN_getDataInt(f, NDOFID, &nd) ) {
        EN_getDataInt(f,info->incorp,&st);
        st++;

        for (i=0; i < nd; i++)
          info->ncvec[count++] = st + i;
      }
    }
    FIter_delete(fIter);
  }

  /* stack region DOFs onto ncvec next*/
  if (info->regnson)  {
    RIter rIter = M_regionIter(mesh);
    while(r = RIter_next(rIter)) {
      if (EN_getDataInt(r, NDOFID, &nd) ) {
        EN_getDataInt(r,info->incorp,&st);
        st++;

        for (i=0; i < nd; i++)
          info->ncvec[count++] = st + i;
      }
    }
    RIter_delete(rIter);
  }//regnson
  
  /* Make sure our count is consistent */
  assert(count == info->nshg);
  /* Fill any remaining entries up to maxshg with zeros */
  while (count < maxshg)
    info->ncvec[count++] = 0;

#ifdef PARALLEL
  /*         sendbuffer ,  count,  type  , recvbuffer , count(perrecv) */
//  MPI_Gather(info->ncvec, maxshg, MPI_INT, info->ncvec, maxshg,
  /*         type   ,  rank of receiver,  communicator    */
//             MPI_INT,       0,           MPI_COMM_WORLD);
#endif

#ifdef DEBUG
  
//  for(i=0; i<maxshg; i++){
//      printf("\n[%d]in writeNCorp:  ncvec[%d] = %d \n", PMU_rank(),i, info->ncvec[i]);
//  }
#endif

  if (PMU_rank()==0)  { /* master processor */
    /* Array should be column majored for fortran */
/*      colmaj = (int*) malloc(sizeof(int)*PMU_size()*maxshg); */
/*      for (i=0; i < PMU_size(); i++)  { */
/*        for (j=0; j < maxshg; j++) */
/*          colmaj[j*PMU_size() + i] = info->ncvec[i*maxshg + j]; */
/*      } */
  /* write the global ncorp.out file */
    if(FortFormFlag) {  /* fortran style */


        printf("\nFortFormFlag feature deprecated: exiting\n");
        
        exit(0);
   

/*  #ifdef sun4_5 */
/*        wncorp_((info->ncvec),&PMU_size(),&maxshg,&(info->nshgTot)); */
/*  #elif sgi */
/*        wncorp_((info->ncvec),&PMU_size(),&maxshg,&(info->nshgTot)); */
/*  #elif ibm6000 */
/*        wncorp((info->ncvec),&PMU_size(),&maxshg,&(info->nshgTot)); */
/*  #endif */
    }//FortFormFlag
/*      else { */
/*        fp = fopen("ncorp.out", "w"); */
/*        fprintf(fp, "%d %d %d\n", PMU_size(), maxshg, info->nshgTot); */
/*        fwrite(colmaj, sizeof(int), PMU_size()*maxshg, fp); */
/*        free(colmaj); */
/*        fputs("\n", fp); */
/*        fclose(fp); */
/*      } */
  }//PMU_rank()==0

  /* following is WRONG in all probability since info->ncvec is row major now,
   * plus only the zeroth proc has the "global" version of it */
  if (SONFATH > 0) { /* SONFATH has the total # of fathers, tfath */

      if(PMU_rank() ==0){
          printf("\n SONFATH feature deprecated: exiting\n");
      }
      exit(0);

/*      vnumnp = (int *)malloc(PMU_size()*sizeof(int)); */

      /* Try to alloc and free stuff in the same function..*/
      /* Use int* instead of int***/

/*      *nsons = (int *)malloc(SONFATH*sizeof(int)); */
/*      *ifath = (int *)malloc(maxshg*PMU_size()*sizeof(int)); */

/*    #if 0  /* This wont work any more coz numnp is not a ptr any more */ 
/*      for(i=0; i< PMU_size(); i++) vnumnp[i] = info->numnp[i]; */
/*    #endif */
/*      SonFath(info->ncvec, vnumnp, *ifath, *nsons, &SONFATH, */
/*              &PMU_size(), &info->nshgTot, &maxshg); */
/*      free(vnumnp); */
  }//SONFATH > 0
}

#ifdef __cplusplus
}
#endif

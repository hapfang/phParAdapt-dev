/* collect the periodic boundary condition array by finding
the degrees of freedom associated with periodic mesh entities */

#include <stdio.h>
#ifndef SIM
#include "MSops.h"
#else
#include "MeshSim.h"
#endif
#include "func.h"

/**********************************************************************/
/* compute the periodic boundary condition array. Periodic boundary   */
/* conditions which reside on other processors are treated as a       */
/* communication, and thus are not included in iper.                  */
/**********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId MYCTID;
extern pMeshDataId NDOFID;
extern pMeshDataId IPERID;

// called in writeEnsaFiles.cc
void attachPeriodicBC(pMesh mesh, globalInfo *info, int *iper)
{
  pEntity ent, em;
  pMatch ment;
  VIter viter;  EIter eiter;  FIter fiter;
  int i, nd;
  int myctint, myctint2, iperint;
  double* iperptr;

  viter = M_vertexIter(mesh);

  printf("\nentering attachPeriodicBC\n"); 


  // "IPER" is set setupGlobalTasks.c
  // based on 
  while (ent = (pEntity) VIter_next(viter))  {
    if ((EN_getDataPtr(ent, IPERID, (void**)&em)) ||  
           ( isPeriodic(ent, &ment) && isOnThisProc(ment) && (em = Match_ent(ment)))) {                             
        EN_getDataInt(em,  MYCTID, &myctint);
        EN_getDataInt(ent, MYCTID, &myctint2);
        iper[myctint2] = myctint;
    }
  }
  VIter_delete(viter);

  if (info->edgeson)  {
    eiter = M_edgeIter(mesh);
    while (ent = (pEntity) EIter_next(eiter))  {
      if (EN_getDataInt(ent, NDOFID, &nd)) {
        if ((EN_getDataPtr(ent, IPERID, (void**)&em)) || (isPeriodic(ent, &ment)
                   && isOnThisProc(ment) && (em = Match_ent(ment))))  {
          for (i=0; i < nd; i++) {
            EN_getDataInt(em,  MYCTID, &myctint);
            EN_getDataInt(ent, MYCTID, &myctint2);
            iper[myctint2+i] = myctint + i;
          }
        }
      }
    }
    EIter_delete(eiter);
  }

  if (info->faceson)  {
    fiter = M_faceIter(mesh);
    while (ent = (pEntity) FIter_next(fiter))  {
      if (EN_getDataInt(ent, NDOFID, &nd))  { 
        if ((EN_getDataPtr(ent, IPERID, (void**)&em)) || (isPeriodic(ent, &ment)
             && isOnThisProc(ment) && (em = Match_ent(ment))))  { 
          for (i=0; i < nd; i++) {
            EN_getDataInt(em,  MYCTID, &myctint);
            EN_getDataInt(ent, MYCTID, &myctint2);
            iper[myctint2+i] = myctint + i;
          }
        }
      }
    }
    FIter_delete(fiter);
  }
}

#ifdef __cplusplus
}
#endif

#ifndef H_Solution
#define H_solution

#include <map>
#include "MeshSim.h"

using namespace std;

void 
attachSoln( char* solfile, 
            pMesh mesh,
            pMeshDataId phSol,
            int ndof= 5,
            int P=1 );

int
attachVPSoln( char *solfile, 
              pMesh mesh,
              pMeshDataId phSol,
              pMeshDataId mark,
              pMeshDataId modes,
              int ndof= 5);
void
solution( pEntity ent,
          pMeshDataId phSol,
          double** sol );

void
numberofmodes( pEntity ent,
               pMeshDataId modes,
               int* sol );

void 
attachSolution( char* solfile, 
                pMesh mesh, 
                map<pEntity, double *>& data, 
                int ndof, 
                int P );

#endif

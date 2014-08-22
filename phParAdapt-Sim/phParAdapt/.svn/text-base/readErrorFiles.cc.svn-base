#include <stdio.h>
#include <fstream>
#include "MeshSimInternal.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#endif
#include "MeshSimAdapt.h"
#include <unistd.h>
#include <string>
#include <strings.h>
#include <fstream>
#include <iostream>
#include "phastaIO.h"
#include <stdlib.h>
#include "mesh_interface.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////
// function that reads in erro files in restart format
// using phastaIO library
////////////////////////////////////////////////////////////////////////////////////////////
void
readErrorFiles(double* nodalErrors, int stepNumber){

    char* iformat = "binary";
    int readStepNum;
    int restart;
    char rfname[40];
    sprintf(rfname,"error.%d.%d",stepNumber,PMU_rank()+1);
    
    openfile_( rfname, "read",  &restart );
    
    // contains: nshg,numVars (here: 10 different "error-types"),lstep
    int iarray[4];
    int isize = 3;
    
    readheader_( &restart, "solution", iarray,
                 &isize, "double", iformat );
    
    
    int nshg = iarray[0]; // linears: numVertices   
    int numVar  = iarray[1];//(here: 10 different "error-types")
    isize = iarray[0]*iarray[1];
    readStepNum = iarray[2];// stepnumber

    if(readStepNum != stepNumber ){
        if(PMU_rank() == 0){

            cerr<<"\nerror in readErrorFiles.cc: stepNumber retrieved differs\n"
                <<"from the one requested\n";
            
        }
        SimPartitionedMesh_stop();
        exit(1);
    }

    double* qlocal = new double [ isize ];
    
    readdatablock_( &restart, "solution", qlocal, &isize,
                    "double" , iformat );
    
    
    // loop num Variables (here: 10 different "error-types")
    for(int i=0; i < numVar; i++)
        // loop num shapefuns
        for(int j=0; j < nshg; j++){
            // first nshg elements of nodalErrors are nodal 
            // values of "1st" error-type (=residual)
            // these values are local to the partition and must be 
            // re-mapped in case a global mesh is employed
            nodalErrors[i*nshg+j] = qlocal[i*iarray[0]+j];
        }
    delete [] qlocal;
    
    closefile_(&restart, "read");
}

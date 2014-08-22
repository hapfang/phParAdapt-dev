#include <stdio.h>
#include <fstream>
#include "MeshSimInternal.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#endif
#include "MeshSimAdapt.h"
#include <unistd.h>
#include <string.h>
#include <strings.h>


// process the errors read in
// here just pick the 1st (and 2nd+3rd) out of 10 possible contributions:
// the heat eqn-residual
// this functions will be further parametrized
// setOfNodalErrors: contains different error measures (10)
// nodalErrors: setOfNodalErrors are mapped onto nodalErrors
void
processError(pMesh mesh, double* setOfNodalErrors, double* nodalErrors, int numOfDifferentErrorMeasures)
{
    int numErrorMeasures = 10;


    double etmp;
    double wght[10];

    // setting the weights manually
    for(int i=0; i<numErrorMeasures; i++){
        wght[i]=0;
        if(i<3)wght[i]=1.0;
    }

    for (int i=0; i < M_numVertices(mesh); i++) {
        etmp = 0.0;
        for ( int j=0; j < numErrorMeasures; j++ ){      
             // combing the error indicators with appropriate   weights
            etmp += wght[j]*setOfNodalErrors[i+j*M_numVertices(mesh)];
        }           
        // the different error measures are combined into
        // ONE single nodal value
        nodalErrors[i] = etmp;
    }
}
    

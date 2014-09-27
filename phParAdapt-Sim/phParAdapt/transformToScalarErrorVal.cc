#include <stdio.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"
#include <stdlib.h>
//  using std::cout;
//  using std::endl;
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif


extern pMeshDataId errorIndicatorID;
extern pMeshDataId phasta_solution;

// nvar is number of error indicators
// as provided by phasta
void
transformToScalarErrorVal(pMesh mesh, int nvar, int option)
{
    // loop over vertices
    pVertex v;
    VIter vIter=M_vertexIter(mesh);

    while(v = VIter_next(vIter)) {

        double *nodalErrorSet;
        double* scalarValue = new double;
	double *nodalSolutionSet;
        if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&nodalErrorSet)){
            
            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(v);
            exit(0);
        }
//Now we have to do the same with solution set
        if(!EN_getDataPtr((pEntity)v, phasta_solution,(void**)&nodalSolutionSet)){

            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(v);
            exit(0);
        }

        *scalarValue = processErrorAG(nodalErrorSet,nodalSolutionSet,nvar,option);
        
        delete [] nodalErrorSet;
        delete [] nodalSolutionSet;
        EN_deleteData((pEntity)v,errorIndicatorID );

        EN_attachDataPtr( (pEntity)v, errorIndicatorID, (void *)
		       scalarValue);

#ifdef DEBUG_NOT_DONE_NOW
        double *testScalarValue;
        if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&testScalarValue)){
            
            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(v);
            exit(0);
        }
//         else{
//             cout<<"\nin transformToScalarErrorVal: attached a scalar"
//                 <<"errorIndicatorID value : "<<testScalarValue[0]<<"\n";
//         }
#endif//DEBUG        


    }
    VIter_delete(vIter);

}

#ifdef __cplusplus
}
#endif
	

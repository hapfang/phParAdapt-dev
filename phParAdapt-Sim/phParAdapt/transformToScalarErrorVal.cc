#include <stdio.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <math.h>
#include "phParAdapt.h"
#include <stdlib.h>
//  using std::cout;
//  using std::endl;
using namespace std;
#ifdef SIM
#include "SimMeshTools.h"
#include "SimAdvMeshing.h"
#endif


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
    double *nodalSolutionSet_other;
    double coordvcur[3]; // note not set if not taking this option 
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
// add  option to compute max pressure gradient from incident edges
	  V_coord(v, coordvcur );
        if(0) {
/*  */
	  double *Test;
          int numEdges = V_numEdges(v);
          pEdge edge;
          pVertex vother;
          double gradp=0;
          for (int i=0; i < numEdges; i++) {
            edge = V_edge(v,i);
            vother=E_otherVertex(edge, v);
            if(!EN_getDataPtr((pEntity)vother, phasta_solution,(void**)&nodalSolutionSet_other)){
            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(vother);
            exit(0);
            }
// max p hard to relate to paraview            gradp=max(gradp,fabs((nodalSolutionSet[0]-nodalSolutionSet_other[0])/E_length(edge)));
            gradp+=fabs((nodalSolutionSet[0]-nodalSolutionSet_other[0])/E_length(edge));
         }
         nodalErrorSet[3]=gradp/numEdges;  // overwrite 4th error field 
       delete [] Test;
/*  */  
     }
            

        *scalarValue = processErrorAG(nodalErrorSet,nodalSolutionSet,nvar,option,coordvcur);
        
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
//whyfail       delete [] nodalSolutionSet_other;

}

#ifdef __cplusplus
}
#endif
	

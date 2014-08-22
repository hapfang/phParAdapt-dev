#include "phParAdapt.h"
#include "phReadWrite.h"
#include "phastaIO.h"
#include <time.h>
#include <iostream>
#include <unistd.h>
#include <string>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalHessianID;
extern int timeStepNumber;
extern char* oformat;     

/////////////////////////////////
// write out the restart files //
// which contain the hessians  //
// the timestep timeStepNumber already is used for a restart
// therefore the number 1111
/////////////////////////////////
void
writeRestartHessians(pMesh mesh )
{
  

  double* nodalHessian;
  pVertex vertex;
  VIter vIter = M_vertexIter( mesh );
  // variables written out
  int nshg_fine = M_numVertices(mesh);
  double* q = new double [ nshg_fine * 6 ];    
  int vcounter=0;
  
  while ( vertex = VIter_next( vIter ) ) {
    if(EN_getDataPtr( (pEntity)vertex, nodalHessianID,
		      (void**)&nodalHessian)) {
      int j=0;
      for( int i=0; i<6; i++ ) {	
	// store the hessians into q vector
	q[vcounter+j] = nodalHessian[i]; 
	
	j=j+nshg_fine;
      }
      vcounter = vcounter+1;  
      
    }
    // there is no hessian on this vertex
    else{
      cerr<<"\nerror in writeRestartHessians: encountered\n"
	  <<"a vertex that carries no hessian\n";
      exit(0);
    }    
  }
  VIter_delete(vIter); 

  char filename[28];
  sprintf(filename,"restart.1111.%d" ,PMU_rank()+1 );

  writeArrayToFile(filename,"hessian",oformat,"write",
                    nshg_fine,6,1111,q);


  

}

#ifdef __cplusplus
}
#endif

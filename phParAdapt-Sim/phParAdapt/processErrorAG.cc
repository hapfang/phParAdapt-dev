#include <iostream>
#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"
#include <math.h>
//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern double* wght;
extern double epsilon;

using namespace std;
// processing error indicator values
// provided by phasta according to your needs
double
processErrorAG(double* nodalErrorSet, double* nodalSolutionSet, int nvar, int option)
{
    double scalarVal = 0;
    double phi;

// option =2 - base error on proximity to interface
    if (option == 2) {
      phi = sqrt(nodalSolutionSet[6]*nodalSolutionSet[6]);
//    cout << "phi = "<<phi<<"\n";
      if (phi < epsilon) {
        scalarVal = 1.0 - 0.5*(1.0-phi/epsilon+1.0/3.141593*sin(3.141593*phi/epsilon)); 
      } else {
        scalarVal = 0.001;
      }
    }
    else if (option == 5) {
       scalarVal += sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]);
    }
    else if (option == 10) {
       scalarVal = nodalErrorSet[3] * log( sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]) + 1E-10 );
    }
    else if (option == 11) {
       double rms_mag = sqrt(nodalErrorSet[6]*nodalErrorSet[6]+nodalErrorSet[7]*nodalErrorSet[7]+nodalErrorSet[8]*nodalErrorSet[8]);
//       scalarVal = rms_mag * log( sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]) + 1E-10 ) / (1e-2 + 100*nodalErrorSet[3]);
// lastused       scalarVal = log( rms_mag*sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]) + 1E-10 ) / (1e-0 + 10*nodalErrorSet[3])+500*nodalErrorSet[4];
//       double pde_mag=sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]);
       double pde_mag=sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]);
       scalarVal =  rms_mag *log(pde_mag + 1E-10 );
    }
    
    else {
    // using linear weights provided by input file adapt.inp
    for(int i=0; i<nvar;i++){
//        cout << nodalErrorSet[i] << " "  << wght[i] << endl;
      }
    }
//
    return scalarVal;
}
        
        



#ifdef __cplusplus
}
#endif
	

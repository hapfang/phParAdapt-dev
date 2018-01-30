#include <iostream>
#include "func.h"
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
processErrorAG(double* nodalErrorSet, double* nodalSolutionSet, int nvar, int option, double* coord)
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
       double rms_mag = 0.0; // HACK KEJ The following will go out of bounds when doing pre-adapt
// since the solution has only 6 fields and error has been pointed at solution for pre-adapt
//  sqrt(nodalErrorSet[6]*nodalErrorSet[6]+nodalErrorSet[7]*nodalErrorSet[7]+nodalErrorSet[8]*nodalErrorSet[8]);
//       scalarVal = rms_mag * log( sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]) + 1E-10 ) / (1e-2 + 100*nodalErrorSet[3]);
// lastused       scalarVal = log( rms_mag*sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]) + 1E-10 ) / (1e-0 + 10*nodalErrorSet[3])+500*nodalErrorSet[4];
//       double pde_mag=sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]);
       double pde_mag=sqrt(nodalErrorSet[0]*nodalErrorSet[0]+nodalErrorSet[1]*nodalErrorSet[1]+nodalErrorSet[2]*nodalErrorSet[2]);
// Previous       scalarVal =  rms_mag *log(pde_mag + 1E-10 );
//      double coord[3];
//      //double plane;
//      V_coord(vertex,coord); 
/*      double ptCheck[3]; // This is a point that should not be ref 
      ptCheck[0]=0.99; //38.025 11.625 4.975
      ptCheck[1]=0.40;
      ptCheck[2]=0.05;
      int stop;
      double distTop = sqrt(dist(coord, ptCheck));
      if(distTop < 5e-2) {
         stop=1;
      }  
*/
     if(0) {   // this was what was used for CRM
       double pgrad_p[3];
       double volInv=1.0/nodalErrorSet[0];
       pgrad_p[0]=nodalErrorSet[3]*volInv;
       pgrad_p[1]=nodalErrorSet[4]*volInv;
       pgrad_p[2]=nodalErrorSet[5]*volInv;
       double maggradp=sqrt(pgrad_p[0]*pgrad_p[0]
                           +pgrad_p[1]*pgrad_p[1]
                           +pgrad_p[2]*pgrad_p[2]);
//       scalarVal =  30.0*nodalSolutionSet[5];
       scalarVal =  90.0*nodalSolutionSet[5];
// flat plate hack       if(0.8*(coord[0]-0.4)-0.45*(coord[1]-17.0) > 0){
       if(1) { //0.8*(coord[0]-40.4)-0.45*(coord[1]-17.0) > 0) {
       scalarVal = max(scalarVal,log10(maggradp + 1E-10 ));
//       scalarVal *= log10(maggradp + 1E-10 );
//Check         scalarVal=max(scalarVal,1.0);  // beyond this plane pick up pgrad without evisc weight but only when it is bigger than the product
       }
    } else {

// Boeing CalTech   Using Dwal and ybar requires us to HACK STRONGLY the creation of the only two vectors this code gets passed-- 
// nodalErrorSet and nodalSolution set.   This is done in  adapt (search for Bad hack)  THIS BREAKS PGRAD ABOVE WHICH DOES SAME HACK IN PHASTA

// pasted from ParaView 
// EVbar*2e7+0.5e-1*pde-res_Z*dwal = 2e5
     scalarVal =  nodalErrorSet[4]*2e7+0.05*nodalErrorSet[2]*nodalErrorSet[3];
// pasted from ParaView with value 0.01 for HLCRM-16degrees Simmmetrix Coarse-Mixed
// mag(rms-vel)*log10(mag(pde-res)+1.0e-10)*EVbar*0.4 + 0.2*EVbari
// Previous Riccardo CRM work      scalarVal =  (rms_mag *log(pde_mag + 1E-10 )*0.4 +0.2)*nodalSolutionSet[5];
//       scalarVal =  nodalSolutionSet[5];
    }
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
	

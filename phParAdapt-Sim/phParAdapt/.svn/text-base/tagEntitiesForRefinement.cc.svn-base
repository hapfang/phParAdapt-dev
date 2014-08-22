#include "phParAdapt.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;

// tag the entities for refinement (for isotropic refinement)
// factor is used to evaluate the threshold for refinement
// as of now do not use hmin and hmax 
void tagEntitiesForRefinement(pMesh mesh,
			      pMSAdapt simAdapter,
			      double factor,
			      double hmax, 
			      double hmin,
			      int option)
{
//  option =0;
    
  double totalError = 0.;
  double maxError = 0.;
  double minError = 0.;
  double threshold = 0.;

  int edgesTagged = applyMarkingStrategy(mesh,simAdapter,factor,hmax,hmin,
					 totalError,maxError,minError,
					 threshold,option);

  cout<<" Info. on adaptation parameters are: "<<endl;
  cout<<" Total Error : "<<totalError<<endl;
  cout<<" Max. Error  : "<<maxError<<endl;
  cout<<" Min. Error  : "<<minError<<endl;
  cout<<" Threshold   : "<<threshold<<endl;
  cout<<" factor : "<<factor<<endl;
  cout<<" hmax   : "<<hmax<<endl;
  cout<<" hmin   : "<<hmin<<endl;
  cout<<"\nNumber of edges tagged to be refined : "<<edgesTagged<<"\n"<<endl;

  std::ofstream adaptSimLog("phAdapt.log");
  adaptSimLog<<"Strategy chosen for adaptation is tag driven"<<endl;
  adaptSimLog<<"(i.e., isotropic refinement)"<<endl;
  adaptSimLog<<"Info. on adaptation parameters are: "<<endl;
  adaptSimLog<<"Total Error : "<<totalError<<endl;
  adaptSimLog<<"Max. Error  : "<<maxError<<endl;
  adaptSimLog<<"Min. Error  : "<<minError<<endl;
  adaptSimLog<<"Threshold   : "<<threshold<<endl;
  adaptSimLog<<"factor : "<<factor<<endl;
  adaptSimLog<<"hmax   : "<<hmax<<endl;
  adaptSimLog<<"hmin   : "<<hmin<<endl;
  adaptSimLog<<"\nNumber of edges tagged to be refined : "<<edgesTagged<<endl;
  adaptSimLog.close();
}

#ifdef __cplusplus
}
#endif

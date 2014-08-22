#ifndef H_EssentialBC
#define H_EssentialBC

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "MeshSimInternal.h"      // for gType in essentialBC
#include "phParAdapt.h"
/* compatibility */
#ifdef SIM
#include "MeshSim.h"
#endif
#include "ModelTypes.h"
#include "BoundaryCondition.h"

using namespace std;
extern map<pair<string,pGEntity>,int> GEntityDataContainerInt;
extern string AXCL;

class EssentialBC : public BoundaryCondition {
public:
  EssentialBC (pGFace gface);
  EssentialBC (pGEdge gedge);
  EssentialBC (pGVertex gvert);
  int eval (pVertex, double *);
  int evalE (pEdge, double *, int);      // for hierarchic basis
  int evalF (pFace, double *, int);      // for hierarchic basis
  void takeBCfromIC(double *BC, double *qTot, int nshgTot, int GDOFnum);
private:
  gType gtype;                            // Gvertex, Gedge or Gface
  pGFace gf;     // could store in a union since either face, edge or vertex
  pGEdge ge;     //   but too much bother for insignificant memory savings
  pGVertex gv;
  int zxcl;      // Flag to do Axisymmetric Centerline

  void update_inherit ();     // to update dontinherit
  void update_thermo(const double*, double*, int*);  // thermo BC's
  void update_velo(const double*, double*, int*);    // velo BC's
  void update_scalar(const double*, double*, int*);  // scalar BC's
  void update_axisym(const double*, double*, int*);  // axisym centerline
  void pteval (const double*, double *, int *);
  int isZA()
  {
    GEntity* gty;
    switch (gtype){
    case Gface:
      gty = (pGEntity) gf;
      break;
    case Gedge:
      gty = (pGEntity) ge;
      break;
    case Gvertex:
      gty = (pGEntity) gv;
      break;
    default:
      cerr << " Entity does not have type "<<endl;
      exit(1);
    }
    int tmp;
    return GEN_dataI(gty,"AXCL",&tmp);
  }
};

#endif

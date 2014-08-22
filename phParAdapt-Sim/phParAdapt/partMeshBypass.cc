#ifdef PARALLEL

/* quick and dirty partition */
#include <iostream>
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"
#include <fstream>
#ifndef SIM
#include "modeler.h"
#endif

using namespace std;

extern int numelX;    // number of elements in x direction
extern int Idirection;
extern int wGraph;

double centroid_x(pRegion region,int Idir);

extern "C" void partMeshBypass(pParMesh pmesh, pPartitioner up, int chkwgt,void *)
{
  pRegion  region;
  int i;
  int Idir = Idirection-1;
  double * maxX;    /* bounding max X of a partition */
  double * minX;    /* bounding min X of a partition */
  int deltaNL;
  double min[3], max[3];

  // this is better than before coz it will be same for all procs, plus of
  // course, it's more efficient
  pGModel model = M_model(pmesh);
  GM_bounds(model, min, max);
  double xmin = min[Idir];
  double xmax = max[Idir];

  int npL = PMU_size();
  int numelL = numelX;

  double deltax = (xmax -xmin)/numelX;

  //Generate partition stencil....

  maxX = new double [PMU_size()];
  minX = new double [PMU_size()];

  deltaNL = ((numelL+1)/npL);

  maxX[0]=xmin+deltax*(deltaNL);

  for(i =1; i< npL-1;i++) {
    maxX[i]=maxX[i-1]+deltax*deltaNL;
  }

  maxX[PMU_size() -1] = xmax;   /* just for sureity */

  minX[0] = xmin;

  for(i=1; i<npL;i++) {
    minX[i] = maxX[i-1];
  }

  pMesh mesh = PM_mesh(pmesh, 0);
  RIter rIter = M_regionIter(mesh);
  ofstream out;
  int cgid = PMU_gid(PMU_rank(), 0);

  if (wGraph == 1) {
    char fn[20];
    sprintf(fn, "partit.out.%d", PMU_rank()+1);
    out.open(fn);
  }
  while(region =  RIter_next(rIter)) {
    double cx = centroid_x(region,Idir);
    for (i=0; i < PMU_size(); i++) {
      int newgid = PMU_gid(i, 0);
      if (minX[i] <= cx && cx < maxX[i])  {

          
        if(PMU_rank()==0){

            cout<<"Partitioner_addRegion() deprecated: exiting\n";
        }
        exit(0);
        // taken out 4/5/05 JM: would have to be setup differently
        // according to MeshSim5.4 Manual
        //Partitioner_addRegion(up, region, newgid, cgid);


        if (wGraph == 1)
          out << i << "\n";
        break;   // oww
      }
    }
  }
  delete[] maxX;  // owww
  delete[] minX;
}

double centroid_x(pRegion region,int dir) {

  pPList v_list;
  double cenx=0;
  double xyz[3];
  v_list = return_vertices(region);
  int n = PList_size(v_list);

  for(int i=0; i< n;i++) {
    V_coord((pVertex)PList_item(v_list,i),xyz);
    cenx = cenx + xyz[dir];
  }
  PList_delete(v_list);
  cenx /= n;
  return cenx;
}

#endif

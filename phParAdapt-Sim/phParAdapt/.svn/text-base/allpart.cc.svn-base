#ifdef PARALLEL
/**************************************************************************
 * does multidirection cuts for partitioning.
 *
 *************************************************************************
 * Format of the input file (partition.in)
 *
 *  nx ny nz
 *  x0 x1 ... xk xk+1 ... ...xnx xn_max
 *  y0 y1 ... yk yk+1 ... ...yny yn_max
 *  y0 z1 ... zk zk+1 ... ...znz zn_max
 *
 * 1. nx , ny and nz NEED NOT BE EQUAL
 *
 * 2. Here nx,ny and nz are the number of slices in each direction
 *    so should be set to zero for directions in which no partitioning
 *    is intended
 *
 * 3. i0 and in_max (i = x,y,z) form the maximum and minimum in the i direction
 *    they should always be specified for every direction. if ni is "0" then
 *    i0 and in_max will be the only values on the (i+1)th line.
 **************************************************************************/

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
// #ifndef SIM
// #include "MSopsInternal.h"
// #else
#include "MeshSimInternal.h"
// #endif
#include "func.h"
#include "ccfunc.h"
#include "parallel.h"
#include "mesh_interface.h"
/* global variables */

using namespace std;

extern int wGraph;

void centroid(pRegion region, double* cen);
void centroid_mod(pRegion region);

typedef struct{
  double min[3];
  double max[3];
} chunk;

int is_Inside(pRegion region, chunk* cube,int inpart);


extern "C" void allpart(pParMesh pmesh, pPartitioner up, int chkwgt, void *)

{
  // if this routine is invoked , there should be a file called
  // partition.in in the current working directory.
    

  double **cut;
  int nx,ny,nz;
  int i,j,k,m;

  FILE* ipf = fopen("partition.in","r");
  if(ferror(ipf)){
    cerr << " Error occured when trying to read [partition.in]"<<endl;
    cerr << " please check if the file exists and is in correct format"<<endl;
    exit(1);
  }
  fscanf(ipf,"%d %d %d", &nx, &ny, &nz);
  cut = new double* [3];
  cut[0]=new double [nx+2];
  cut[1]=new double [ny+2];
  cut[2]=new double [nz+2];

  for(i=0; i< nx+2;i++) fscanf(ipf,"%lf",&cut[0][i]);
  for(i=0; i< ny+2;i++) fscanf(ipf,"%lf",&cut[1][i]);
  for(i=0; i< nz+2;i++) fscanf(ipf,"%lf",&cut[2][i]);

  fclose(ipf);

  //  cut[i][j] gives the jth cut location in the  ith direction.
  // also nx, ny, nz give number of cuts in x y and z;

  pRegion region;
  int npart = (nx+1)*(ny+1)*(nz+1);
  chunk *cube;
  cube = new chunk [ npart];

  // creating the partition templates using the maximun and minimum
  // co-ordinate values for each direction.

  for(k=0; k<= nz; k++){
    for(j=0; j<= ny; j++){
      for(i=0; i<= nx; i++){
       m = k*(nx+1)*(ny+1)+j*(nx+1)+i;
       cube[m].min[0]=cut[0][i];
       cube[m].min[1]=cut[1][j];
       cube[m].min[2]=cut[2][k];
       cube[m].max[0]=cut[0][i+1];
       cube[m].max[1]=cut[1][j+1];
       cube[m].max[2]=cut[2][k+1];
      }
    }
  }

  // the actual partitioning, ie identifying which partition each region
  // belongs to.

  pMesh mesh = PM_mesh(pmesh, 0);
  int cgid = PMU_gid(PMU_rank(), 0);
  RIter rIter = M_regionIter(mesh);
  ofstream out;

  if (wGraph == 1)  {
    char fn[20];
    sprintf(fn, "partit.out.%d", PMU_rank()+1);
    out.open(fn);
  }
  while(region = RIter_next(rIter)) {
    int pid = is_Inside(region,cube,npart);
    int newgid = PMU_gid(pid, 0);

    if(PMU_rank()==0){

        cout<<"Partitioner_addRegion() deprecated: exiting\n";
    }
    exit(0);
    // taken out 4/5/05 JM: would have to be setup differently
    // according to MeshSim5.4 Manual
            //Partitioner_addRegion(up, region, newgid, cgid);
    out << pid << "\n";
  }
  RIter_delete(rIter);
  delete[] cut[0]; delete[] cut[1]; delete[] cut[2];
  delete[] cut;
  delete[] cube;
}

void centroid(pRegion region, double* cen) {
  pPList v_list;
  double* xyz= new double[3];
  int Idir;
  for(Idir =0; Idir < 3; Idir++){
    xyz[Idir]=0.0;
    cen[Idir]=0.0;
  }

  int nen;


      v_list = R_vertices(region,1);


  nen = PList_size(v_list);
  for(int i=0; i< nen;i++) {
    V_coord((pVertex)PList_item(v_list,i),xyz);
    for(Idir=0; Idir<3; Idir++)
      cen[Idir] = cen[Idir] + xyz[Idir];
  }
  delete [] xyz;
  PList_delete(v_list);
  for(Idir=0; Idir<3; Idir++) cen[Idir] = cen[Idir]/nen;

}
void centroid_mod(pRegion region) {
  pPList v_list;

  double* xyz= new double[3];
  double* cen= new double[3];
  int Idir;
  for(Idir =0; Idir < 3; Idir++){
    xyz[Idir]=0.0;
    cen[Idir]=0.0;
  }
  int nen;


  v_list = R_vertices(region,1);



  nen = PList_size(v_list);
  for(int i=0; i< nen;i++) {
    V_coord((pVertex)PList_item(v_list,i),xyz);
    //debug
    cout << i<<": "<<xyz[0]<<"  "<<xyz[1]<<"  "<<xyz[2]<<endl;
    //debug
    for(Idir=0; Idir<3; Idir++)
      cen[Idir] = cen[Idir] + xyz[Idir];
  }
  delete [] xyz;
  delete [] cen;
  PList_delete(v_list);
  // debug
  cout << "centroid : "<<cen[0]<<"  "<<cen[1]<<"  "<<cen[2]<<endl;
  //debug
  for(Idir=0; Idir<3; Idir++) cen[Idir] = cen[Idir]/nen;
  // debug
  cout << "centroid : "<<cen[0]<<"  "<<cen[1]<<"  "<<cen[2]<<"  "<<nen<<endl;
  //debug

}

int is_Inside(pRegion region, chunk* cube, int npart){

   double* xyz=new double [3];
   int pid=362436;
   static int count=0;
   centroid(region,xyz);

   for(int i=0; i< npart; i++){
     if((xyz[0] >= cube[i].min[0]) && (xyz[0] <= cube[i].max[0])){
       if((xyz[1] >= cube[i].min[1]) && (xyz[1] <= cube[i].max[1])){
         if((xyz[2] >= cube[i].min[2]) && (xyz[2] <= cube[i].max[2])){
           pid = i;
           break;
         }
       }
     }
   }

   if (pid > npart) {
     cout << "Error in isInside, pid: "<< pid<<" "<< count++<<endl;
     cout << xyz[0]<<" "<< xyz[1]<<" "<<xyz[2]<<endl;
     getchar();
     centroid_mod(region);
     getchar();
   }
   delete [] xyz;
   return pid;
}

#endif

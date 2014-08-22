#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "mpi.h"

#include "phParAdapt.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern int BCInflowFaceTag;

pMeshDataId vertID;
pMeshDataId inflowOwnedVertID;

void BCInflowFaceInfo(pGModel model, pParMesh pmesh, pMesh mesh) {

  vertID =  MD_newMeshDataId( "vertex id.");
  inflowOwnedVertID =  MD_newMeshDataId( "inflow vertex");

  int nodesLoc = BCInflowFaceNodesInfo(model,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  commuBCInflowFaceInfo(pmesh,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  int facesLoc = BCInflowFaceConnectInfo(model,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  cout<<"\n ["<<PMU_rank()<<"] "<<nodesLoc<<" : local owned nodes found (BCInflowFaceInfo)\n";
  cout<<" ["<<PMU_rank()<<"] "<<facesLoc<<" : local owned faces found (BCInflowFaceInfo)\n";

  int nodesTot, facesTot;
  MPI_Allreduce(&nodesLoc, &nodesTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&facesLoc, &facesTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  // combine all the procs
  if(PMU_rank()==0)
    BCInflowFaceGlobalInfoInVtk(nodesTot, facesTot);

  int id;
  pVertex vertex;
  VIter vit = M_vertexIter(mesh);
  while(vertex = VIter_next(vit)) {
    if(EN_getDataInt((pEntity)vertex,vertID,&id)) {
      EN_deleteData((pEntity)vertex,vertID);
      EN_deleteData((pEntity)vertex,inflowOwnedVertID);
    }
  }
  VIter_delete(vit);

  MD_deleteMeshDataId(vertID);
  MD_deleteMeshDataId(inflowOwnedVertID);

}

int BCInflowFaceNodesInfo(pGModel model, pMesh mesh) {

  char filename[256];
  sprintf(filename,"inflow_mesh_points.%i.dat",PMU_gid(PMU_rank(),0)+1);

  ofstream fout(filename);  

  // need to put this in global var. list
  int tag = BCInflowFaceTag;
  int id, gtag, count = 0;
  // we can used classified vertex iter on model face with tag="tag"
  // but in case tag was not specified correctly it would 
  // crash the code (after adapt and before preprocessing)
  pVertex vertex;
  pFace face;
  double xyz[3];
  FIter fit = M_faceIter(mesh);
  while(face = FIter_next(fit)) {
    gtag = GEN_tag(F_whatIn(face));
    // if the model face' tag is  the  one specified
    if(gtag == tag) {
      pPList fVerts = F_vertices(face,0);
      for(int i=0; i<PList_size(fVerts); i++) {
	vertex = (pVertex)PList_item(fVerts,i);
	if(EN_isOwnerProc((pEntity)vertex)) {
          if(!EN_getDataInt((pEntity)vertex,vertID,&id)) {
	    count++;
	    EN_attachDataInt((pEntity)vertex,vertID,count);
	    EN_attachDataInt((pEntity)vertex,inflowOwnedVertID,1);
	    V_coord(vertex,xyz);
	    fout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n";
	  }
	}
      }
      PList_delete(fVerts);
    }
  }
  FIter_delete(fit);

  fout.close();

  return count;
}

int BCInflowFaceConnectInfo(pGModel model, pMesh mesh) {

  char filename[256];
  sprintf(filename,"inflow_mesh_connect.%i.dat",PMU_gid(PMU_rank(),0)+1);

  ofstream fout(filename);  

  // need to put this in global var. list
  int tag = BCInflowFaceTag;
  int id, gtag, count = 0, offsetVtk = 1;
  pVertex vertex;
  pFace face;
  FIter fit = M_faceIter(mesh);
  while(face = FIter_next(fit)) {
    if(!EN_isOwnerProc((pEntity)face))
      continue;

    gtag = GEN_tag(F_whatIn(face));
    // if the model face' tag is  the  one specified
    if(gtag == tag) {
      pPList fVerts = F_vertices(face,0);
      for(int i=0; i<PList_size(fVerts); i++) {
	vertex = (pVertex)PList_item(fVerts,i);
	if(!EN_getDataInt((pEntity)vertex,vertID,&id)) {
	  cout<<"\nError in BCInflowFaceConnectInfo : vertex id. not attached\n"<<endl;
	  exit(0);
	}
	fout<<id-offsetVtk<<" ";
      }
      PList_delete(fVerts);
      fout<<"\n";
      count++;
    }
  }
  FIter_delete(fit);

  fout.close();

  return count;
}

void BCInflowFaceGlobalInfoInVtk(int nodesTot, int facesTot) {

  ofstream fout("inflow_mesh_face.vtk");

  // write the file
  fout<<"# vtk DataFile Version 3.0\n";
  fout<<"vtk output\n";
  fout<<"ASCII\n";
  fout<<"DATASET POLYDATA\n";
  fout<<"POINTS "<<nodesTot<<" float\n";

  double tmpDbl;
  for(int i=0; i<PMU_size(); i++) {
    char in_filename[128];
    sprintf(in_filename,"inflow_mesh_points.%i.dat",i+1);

    ifstream fin(in_filename);

    if(!fin) {
      cout<<"\n BCInflowFaceGlobalInfoInVtk() : error opening file : "<<in_filename<<"\n"<<endl;
      exit(0);
    }

    while(!fin.eof()) {
      fin>>tmpDbl;
      if(fin.eof())
	break;

      fout<<tmpDbl<<" ";
      for(int j=0; j<2; j++) {
	fin>>tmpDbl;
	fout<<tmpDbl<<" ";
      }
      fout<<"\n";
    }

    fin.close();
    char rm_file[256];
    sprintf(rm_file,"rm %s",in_filename);
    system(rm_file);
  }

  fout<<"\n";
  fout<<"POLYGONS "<<facesTot<<" "<<4*facesTot<<"\n";
  fout<<"\n";

  int tmp;
  for(int i=0; i<PMU_size(); i++) {
    char in_filename[128];
    sprintf(in_filename,"inflow_mesh_connect.%i.dat",i+1);

    ifstream fin(in_filename);

    if(!fin) {
      cout<<"\n BCInflowFaceGlobalInfoInVtk() : error opening file : "<<in_filename<<"\n"<<endl;
      exit(0);
    }

    while(!fin.eof()) {
      fin>>tmp;
      if(fin.eof())
	break;

      fout<<"3 ";
      fout<<tmp<<" ";
      // already read the first one and then read the next two
      // (assuming triangular faces)
      for(int j=0; j<2; j++) {
	fin>>tmp;
	fout<<tmp<<" ";
      }
      fout<<"\n";
    }

    fin.close();
    char rm_file[256];
    sprintf(rm_file,"rm %s",in_filename);
    system(rm_file);
  }

  // THIS IS NOT USED/REQUIRED FOR BCT-GENERATION !!!
  // dump final mapping
  fout<<"\n";
  fout<<"CELL_DATA "<<facesTot<<"\n";
  fout<<"POINT_DATA "<<nodesTot<<"\n";
  fout<<"SCALARS scalars float\n";
  fout<<"LOOKUP_TABLE default\n";

    // loop over unique vertex list
  for(int i=0; i<nodesTot; i++) {
    fout<<i+1<<"\n";
  }

  fout.close();

  cout<<"\n "<<nodesTot<<" : total nodes found (BCInflowFaceInfo)\n";
  cout<<" "<<facesTot<<" : total faces found (BCInflowFaceInfo)\n";
}
                
#ifdef __cplusplus
}
#endif    



    
    

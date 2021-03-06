#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "mesh_interface.h"
#include "phParAdapt.h"
#ifdef FMDB
#include "ParUtil.h"
#include "AOMD_Internals.h"
#endif
using namespace std; 

void exportPVTK(const char* fname, int numPtn, int numdof)
{
  char realName[256];
  sprintf(realName,"%s.pvtu",fname); 

   FILE* outFile = fopen(realName, "w");  
   if (!outFile)
       cerr<<"unable to create pvtu File\n";

  fprintf(outFile, "<?xml version=\"1.0\"?>\n"); 
  fprintf(outFile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(outFile, "<PUnstructuredGrid GhostLevel=\"0\">\n");

  fprintf(outFile, "<PPointData>\n");
  fprintf(outFile, "<PDataArray type=\"Float32\" Name=\"field\" format=\"ascii\" NumberOfComponents=\"%d\"/>\n", numdof);
  fprintf(outFile, "</PPointData>\n");
  fprintf(outFile, "<PPoints>\n<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n </PPoints>\n"); 

  fprintf(outFile, "<PCells>\n");
  fprintf(outFile, "<PDataArray type=\"Int32\"  Name=\"connectivity\" format=\"ascii\"/>\n");
  fprintf(outFile, "</PCells>\n");

  for (int pid=0; pid<numPtn; pid++)
    {
      fprintf(outFile, "<Piece Source=\"%s%d.vtu\"/>\n", fname,pid); 
    }

  fprintf(outFile, "</PUnstructuredGrid>\n</VTKFile>\n");
   fclose(outFile);    
}


void exportVTK(const char* fname, pMesh mesh, pMeshDataId field, int numdof) 
{

   //return;

   FILE* outFile = fopen(fname, "w");  
   if (!outFile)
       cerr<<"unable to create pvtu File\n";

   int numVtc=M_numVertices(mesh);
   int numRgn=M_numRegions(mesh); 
   
   char tmpTagName[256]; 
   sprintf(tmpTagName, "tmp%s", fname);
//   printf("%s\n",tmpTagName);
   pMeshDataId tagId = MD_newMeshDataId(tmpTagName);
   
   int mypid = PMU_rank(); 

  fprintf(outFile, "<?xml version=\"1.0\"?>\n"); 
  fprintf(outFile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(outFile, "<UnstructuredGrid>\n"); 
  
  fprintf(outFile, "<Piece NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n", numVtc, numRgn);
  fprintf(outFile, "<Points>\n "); 
  fprintf(outFile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  pEntity ent;
  VIter viter = M_vertexIter(mesh);
  int i=0, id; 
  while( ent = (pEntity)VIter_next(viter)) {
    double loc[3];  
    EN_attachDataInt(ent, tagId, i);
    //  fprintf(outFile, " \n%d", EN_id(ent)); 
    V_coord((pVertex)ent,loc);  
    fprintf(outFile, "  %.15f  %.15f  %.15f  \n", loc[0], loc[1], loc[2]);  
    i++;
  }
  VIter_reset(viter);

  char dof[16];
  sprintf(dof,"%d",numdof);
 
  // Header for vtk
  fprintf(outFile, "\n</DataArray>\n"); 
  fprintf(outFile, "</Points>\n"); 

  fprintf(outFile, "<PointData>\n");
  fprintf(outFile, "<DataArray type=\"Float32\" Name=\"field\" NumberOfComponents=\"%d\" format=\"ascii\">\n", numdof);
  double *solution;
  while( ent = (pEntity)VIter_next(viter)) {
      EN_getDataPtr(ent, field, (void**)&solution);
      for(int i=0;i<numdof;i++)
          fprintf(outFile, "  %.15f ", solution[i]);
      fprintf(outFile, "\n");
  }
  
  VIter_reset(viter);
  fprintf(outFile, "\n</DataArray>\n"); 
  fprintf(outFile, "</PointData>\n");
 
  // Header for vtk
  fprintf(outFile, "<Cells>\n");  
  fprintf(outFile, "<DataArray type=\"Int32\"  Name=\"connectivity\" format=\"ascii\">\n");  

  RIter riter = M_regionIter(mesh);
  //int numV[numRgn];
  int *numV = NULL;
  //numV = new int [numRgn];
  numV = new int [numRgn+1];
  pRegion region;
  numV[0]=0; 
  //int type[numRgn];
  int * type = NULL;
  //type = new int[numRgn];
  type = new int[numRgn+1];
  int iNum=0;
#ifdef FMDB 
  while(ent = (pEntity)RIter_next(riter) ) {
    void *iter=0;  
    iNum++; 
    type[iNum]=M_GetElementType(ent);
    pEntity vt;
    pPList vtxs = R_vertices((pRegion)ent,1); 
    numV[iNum]=PList_size(vtxs)+numV[iNum-1];
    while ((vt = (pEntity)PList_next(vtxs,&iter))) 
      {	   
	int id; 
	EN_getDataInt(vt, tagId, &id);
	fprintf(outFile, " %d", id); 
      }
    PList_delete(vtxs);
  }
  RIter_delete(riter);
#endif
#ifdef SIM
  while(ent = (pEntity)RIter_next(riter) ) {
    void *iter=0;
    iNum++;
    pEntity vt;
    pPList vtxs = R_vertices((pRegion)ent,1);
    numV[iNum]=PList_size(vtxs)+numV[iNum-1];
    while ((vt = (pEntity)PList_next(vtxs,&iter)))
      {
   int id;
   EN_getDataInt(vt, tagId, &id);
   fprintf(outFile, " %d", id);
      }
    PList_delete(vtxs);
  }
  RIter_reset(riter);
  iNum=0;
  while(region = RIter_next(riter)){
     iNum++;
     type[iNum]=R_topoType(region);
  }
  RIter_delete(riter);
#endif

  fprintf(outFile, "\n</DataArray>\n");
  fprintf(outFile, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        
  for (int i=1; i<numRgn+1; i++)
      fprintf(outFile, " %d ", numV[i]);  


  fprintf(outFile, "\n</DataArray>\n"); 
  fprintf(outFile, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"); 

#ifdef FMDB
  for (int i=1; i<numRgn+1; i++)
//  for (int i=0; i<numRgn; i++)
    switch(type[i]) {
    case TET :
      fprintf(outFile, " 10 ");
      break; 
    case PYRAMID  :
      fprintf(outFile, " 14 ");
      break;
    case PRISM :
      fprintf(outFile, " 13 ");
      break;
    case HEX :
      fprintf(outFile, " 12  ");
      break; 

    default:
      cout<<"region topology ["<<i<<"] NOT supported"<<endl;
    exit(0);
    }
#endif  
#ifdef SIM
  for (int i=1; i<numRgn+1; i++)
//  for (int i=0; i<numRgn; i++)
     switch(type[i]) {
    case Rtet :
      fprintf(outFile, " 10 ");
      break;
    case Rpyramid  :
      fprintf(outFile, " 14 ");
      break;
    case Rwedge  :
      fprintf(outFile, " 13 ");
      break;
    case Rhex  :
      fprintf(outFile, " 12  ");
      break;

    default:
      cout<<"region topology ["<<i<<"] NOT supported"<<endl;
    exit(0);
    }
#endif

  // Header for vtk
  fprintf(outFile, "\n</DataArray>\n");  
  fprintf(outFile, "</Cells>\n"); 
/*
  fprintf(outFile, "<CellData Scalars=\"mypid\">\n");
  fprintf(outFile, "<DataArray type=\"Int32\" Name=\"mypid\" format=\"ascii\">\n"); 

  for (int i=1; i<numRgn+1; i++)
    fprintf(outFile, " %d ", mypid); 
  fprintf(outFile, "\n</DataArray>\n");
  fprintf(outFile, "</CellData>\n");  
*/
  fprintf(outFile, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");

  while( ent = (pEntity)VIter_next(viter)) {
     if(EN_getDataInt(ent, tagId, &id)){
       EN_deleteData(ent, tagId);
     }
  }
  VIter_delete(viter);
  
 fclose(outFile);  

 delete[] type;
 delete[] numV;

}

// fname should be given without extension. 
void M_writeVTKFile (pMesh mesh,const char * fname, pMeshDataId field, int numdof)
{

///**********************************
///   writing XML VTK files
///**********************************

 int numPtn = PMU_size(); 
 int mypid = PMU_rank(); 
 int dim = 3;
// int dim=M_getMaxDim(mesh);

 if(dim<3)
   {
     printf("2D mesh NOT supported\n");
     exit(0);
   }
 if(numPtn>1)
   {
     // piece files 
     for (int pid=0; pid<numPtn; pid++)
       {
	 char realName[256];
	 sprintf(realName,"%s%d.vtu",fname,pid);
	 if (mypid==pid) 
	   exportVTK(realName, mesh, field, numdof);  
       }  
     
     // PUnstructuredGrid file
     if (mypid == 0)
       { 
	 exportPVTK(fname,numPtn, numdof); 
       }
   }
 else
   {
     char realName[256];
     sprintf(realName,"%s.vtu",fname);
     exportVTK(realName, mesh, field, numdof);  
   }

}


#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "phParAdapt.h"
#include "MeshSimInternal.h"
#include "SimMeshTools.h"
#include "func.h"

extern double rhoinp;
extern double muinp;
extern pMeshDataId wallStressID;

using namespace std;

void CalcYPlus(pMesh mesh) {
  
   VIter vIter  = M_vertexIter(mesh);
   pVertex vert;
   FILE* fyplus;

   fyplus = fopen("yplus.dat", "w");

   int numWallVert = 0;
   double yplusSum = 0.0;

   while (vert = VIter_next(vIter))  { 

      if(V_isOriginatingNode(vert)) {
         
        double node1[3], node2[3];
        pVertex vert2;
        pPList VGrowth = V_growthCurveVertices(vert);

        vert2 = (pVertex)PList_item(VGrowth, 1);

        V_coord(vert, node1);
        V_coord(vert2, node2);
        
        double yDist = sqrt(dist(node1, node2));

        double *WallStress;
        EN_getDataPtr(vert, wallStressID, (void**)&WallStress);

        double uTau = sqrt(WallStress[0]/rhoinp);
        double nu = muinp/rhoinp;
        double yplus = yDist*uTau/nu;

        fprintf(fyplus, "%lf\n", yplus);

        yplusSum = yplusSum + yplus;

        PList_delete(VGrowth);
        numWallVert ++; 
      }
   }
   fclose(fyplus);
   double yplusAvg = yplusSum/numWallVert;
   printf("Average yplus for First cell: %lf\n", yplusAvg); 

   VIter_delete(vIter); 
}

void VolumeDebug(pMesh mesh) {
   
   char file[250];
   sprintf(file, "negVol%d",PMU_rank());
   pRegion region; 
   vector<pRegion> RList;
   RIter rIter = M_regionIter(mesh);
   int haveNeg = 0;
   while( region = RIter_next(rIter)) {
      double vol;
      vol = R_VolumeSim(region);
/*
      int rTopoType = region->getType();
      double weight;
      FMDB_Ent_GetWeight(region, &weight);

       pPList vertices;
       int size;
       vertices = R_vertices(region, 1);
       size = PList_size(vertices);
      
      printf("Type: %d, weight: %f, numNodes: %d \n",rTopoType, weight, size);
*/      
      if(vol<1e-17) {
         R_info(region);
         RList.push_back(region);
         haveNeg = 1;
      }
   }
//   if(haveNeg == 1) 
///      write_3Dpart_vtu(mesh_instance, file, RList);
   RIter_delete(rIter);
}
/*
void NormalizedEdgeLength(pMesh pmesh) {
 
 EIter eit=M_edgeIter(pmesh);
 pEdge edge;
 ofstream nfile;
 
 nfile.open("normEdgeLength.dat");
 while( edge=EIter_next(eit) ) {
    double normE, normESq;

    normE = E_normalizedLength(edge);
    normESq = E_normalizedLengthSq(edge);
    
    nfile << normE << " " << normESq << endl; 
 }
 nfile.close();
}
*/
 //to calculate minimum and maximum edge lengths
void EdgeLength(pMesh pmesh) {
 EIter eit=M_edgeIter(pmesh);
 pEdge edge;
 int numEdges = 0;
 double local_uniformSize = 0., local_len[2] = {1.0e10, 0.};
 while( edge=EIter_next(eit) ) {
   numEdges++;
   double len = sqrt(E_lengthSq(edge));
   local_uniformSize += len;
   if(len<local_len[0])
      local_len[0] = len;
   if(local_len[1]<len)
      local_len[1] = len;
 }
EIter_delete(eit);
printf("min & max edge length: %.16e %.16e\n",local_len[0],local_len[1]);
}
//end
                                                                              


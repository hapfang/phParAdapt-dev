#include <cstdlib>
#include <cstring>
#include "phParAdapt.h"
#include <iostream>
#include <fstream>
#include "func.h"
#include "SimAdvMeshing.h"
#include <math.h>
#include <assert.h>

extern pMeshDataId wallStressID;
extern pMeshDataId ybarID;
extern pMeshDataId isOrgNodeID;
extern double rhoinp;
extern double muinp;

//the variable nuYrho is actually nu*Y1plus*sqrt(rho)
double DispThickness(pVertex vertex) {
   double node1[3], node2[3], Dely, uMag;
   double* nodalData;
   pVertex v;
   vector<pVertex> GC;
//   BL_getGrowthCurveNodes(vertex, GC);
   int GCSize =  GC.size();
   double DispThick = 0.0;
   double uFrStrm;  //for M6wing
   V_coord(vertex, node1);
   
   v=GC[GCSize-1];
   EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalData);
   
   uFrStrm = sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]); 

   for(int i=1; i<GCSize; i++) { 
      v=GC[i];
      if(!EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalData)){
         cout<<"\nerror in WallStress: no data attached to  vertex\n";
         V_info(v);
         exit(0);
      }
      
      V_coord(v, node2);
      Dely = sqrt(dist(node1, node2));
      
      uMag = sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]);
      
      DispThick += (1-uMag/uFrStrm)*Dely;
      node1[0]=node2[0]; node1[1]=node2[1]; node1[2]=node2[2];
      double* Dist = new double[1];

//      EN_deleteData((pEntity)vertex, wallDistID, (void *) Wa
      
   }
//   printf("Total BL height: %lf\n", DispThick);
   return DispThick;
}

      
void SpaldingLaw(double ydist, double u, double& uTau, double& yplus) {
   double kup, kup2, kup3, kappa;
   
   double mu=muinp; 
   double rho=rhoinp;
   double err=1e-6;
   double rat = 1.0;
   double efac=0.1108;
   double uplus, yrmi, f, dfds;
   double nu; nu = mu/rho;

   uTau=0.04;
   yrmi=ydist/nu;
   kappa=0.35;

   int i;

   i=0;
   while(fabs(rat)>err || i<500) {
      yplus = yrmi*uTau;
      uplus = u/uTau;
      kup = kappa*uplus;
      kup2 = kup*kup;
      kup3 = kup2*kup;

      f = uplus-yplus+efac*(exp(kup)-1.0-kup-kup2*0.5-kup3/6);
      dfds = uplus+yplus+efac*(exp(kup)*kup-kup-kup2-kup3*0.5);
      rat = f*uTau/dfds;
      uTau = uTau + rat;
      i++;
   }
   if(fabs(rat)>err)
     printf("failed to converge on uTau %lf\n");
}

void WallStress(pMesh mesh) {
      
      double uTau;
      pPList BaseVtxList; 
      pVertex BaseVert;
      double* WallStress = new double[1];
      int isOrg = 0;

      double nu = muinp/rhoinp;
      double rho = rhoinp;

      VIter vIter  = M_vertexIter(mesh);
      pVertex vert;
      BaseVtxList = PList_new();

      while (vert = VIter_next(vIter))  {
       
       WallStress[0] = 0.0;  
       EN_attachDataPtr((pEntity)vert, wallStressID,(void *)WallStress);
       EN_attachDataInt((pEntity)vert, isOrgNodeID, isOrg);
       if(EN_isBLEntity(vert)) {
            
          double node1[3], node2[3];
          pPList VGrowth; 
          pVertex v;
          int GCLoopSize;

          VGrowth = V_growthCurveVertices(vert);
          int Size = PList_size(VGrowth);
            
          BaseVert = (pVertex)PList_item(VGrowth, 0);
//          PList_delete(VGrowth);
          int FoundVtx = PList_inList(BaseVtxList, BaseVert);
          if(!FoundVtx) {
             PList_append(BaseVtxList, BaseVert);
             int BaseOrg = 1; 
            
            double uMag, VinpMag, yDist, Vw;
            double Vel[3], wvec[3], wnorm[3], Vinpl[3];
            uTau = 0.0;

            V_coord(BaseVert, node1);
            GCLoopSize = 4;
            for(int i=1; i<GCLoopSize; i++) {
              if(Size>GCLoopSize) { 
                double* nodalData;
                v = (pVertex)PList_item(VGrowth, i);
                if(!EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalData)){
                  cout<<"\nerror in WallStress: no data attached to  vertex\n";
                  V_info(v);
                  exit(0);
               }
               
               V_coord(v, node2); 
               yDist = sqrt(dist(node1, node2));
               uMag = sqrt(nodalData[0]*nodalData[0]+nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]);
/*
//calculate the inplane component of velocity for wall shear stress
               Vel[0]=nodalData[0]; Vel[1]=nodalData[1]; Vel[2]=nodalData[2];
//calculate wall normal vector
               diffVt(node2, node1, wvec);
               normVt(wvec, wnorm);
               Vw = dotProd(wnorm, Vel);

               for(int j=0; j<3; j++) {
                  Vinpl[j] = Vel[j] - Vw*wnorm[j];
               }

//we are gonna use in plane velocity magnitude for wall shear stress          
//         VinpMag =sqrt(Vinpl[0]*Vinpl[0]+Vinpl[1]*Vinpl[1]+Vinpl[2]*Vinpl[2]);
*/
               double yplus, utau;
//               SpaldingLaw(yDist,uMag, utau, yplus); 
//               utau = sqrt(nu*VinpMag/yDist);
               utau = sqrt(nu*uMag/yDist);
               uTau += utau;
//              uTau = 1.785e-5*uMag/yDist;
              }
             }
             PList_delete(VGrowth);
            
             uTau = uTau/(GCLoopSize-1);
             if(Size == 1)  uTau = 0.0;

/*//           following only for compressible flows  
             double *nodalData;
             EN_getDataPtr((pEntity)BaseVert, ybarID,(void**)&nodalData);
             rho = nodalData[3]/(287.05*nodalData[4]);
*///           end

             WallStress[0] = uTau*uTau*rho;
             double *WallCopy = new double[1];
             WallCopy[0] = WallStress[0];
//            printf("Wall Stress %lf\n",WallStress[0]);

            EN_modifyDataPtr((pEntity)BaseVert, wallStressID,(void *)WallCopy);
            EN_modifyDataInt((pEntity)BaseVert, isOrgNodeID, BaseOrg);

           }  // ifFoundVtx
         }   //isBLEntity
       }      //VIter 
#ifdef DEBUG      
       M_writeVTKFile(mesh, "WallStress", wallStressID, 1);
//       M_writeVTKFile(mesh, "OrgNodes", isOrgNodeID, 1); 
#endif      
}


#include <cstdlib>
#include <cstring>
#include "math.h"
#include "phParAdapt.h"
#include <iostream>
#include "SimAdvMeshing.h"
#include "func.h"
#include <iostream>

extern double rhoinp;
extern double muinp;

using namespace std;

//Y1plus value for RANS(spallart Almaras = 1.0-5.0, 1.0 is better
//Y1plus value for k-eps and k-omega models 0.1-0.5
//Y1plus for Wall model 30-50  (refer literature for best suitable value)

extern pMeshDataId wallStressID;
extern pMeshDataId wallDistID;
extern pMeshDataId nodalVorticityID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId localGradientID;
extern pMeshDataId localPatchVolID;
extern pMeshDataId isOrgNodeID;
extern pMeshDataId SeparatedID;
extern pMeshDataId ybarID;
extern pMeshDataId ShockID;

void TotalBLHeight(pMesh mesh) {
     
     pVertex vertex, v;
     pPList VGrowth;
     double uFrStrm, uMag;
     double node1[3], node2[3];
//     double wvec[3], wnorm[3], Vortinpl[3];
     double* nodalData;
     double* Sep = new double[1];
     Sep[0] = 0.0;
//     double threshold = 0.001;
//     double threshold = 1.0;

     VIter vit=M_vertexIter(mesh);
     while(vertex = VIter_next(vit)){
        double* Dist = new double[3];
//Dist[0] is first point height, Dist[1] is total BL height, Dist[2] is displacement thickness       
        Dist[0] = 0.0; Dist[1] = 0.0; Dist[2] = 0.0;
        EN_attachDataPtr((pEntity)vertex, wallDistID, (void *)Dist);
        EN_attachDataPtr((pEntity)vertex, SeparatedID, (void *)Sep);

       if(V_isOriginatingNode(vertex)) {
//         int Gtype; 
         VGrowth = V_growthCurveVertices(vertex);
         int Size = PList_size(VGrowth);

         double height = 0.0;
//         height = BLHeightVel(VGrowth);
         height = BLHeightVort(vertex, VGrowth);
      
      v = (pVertex)PList_item(VGrowth, Size-1);
      EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalData);
      uFrStrm = sqrt(nodalData[0]*nodalData[0]+nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]);           

      double DispThick = 0.0;
      double Dely; 
      for(int i=1; i<Size; i++) {
        v = (pVertex)PList_item(VGrowth, i);
         EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalData);

         V_coord(v, node2);
         Dely = sqrt(dist(node1, node2));
         
         uMag = sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]);
         
         DispThick += (1-uMag/uFrStrm)*Dely;
//         DispThick += (1-uMag/uFrStrm)*Dely;
         node1[0]=node2[0]; node1[1]=node2[1]; node1[2]=node2[2];
//         printf("Disp thick %d: %lf\n", i, DispThick);
      } //loop of pPlist

       double SepH = 0.0;
       SepH = SeparatedBL(vertex);
//       double interH = GrowthCurveIntersect(vertex);

       if(SepH > 0.0) {
          height = SepH;
       }
//         cout << height << endl; 

/*
       if(interH > 0.0 && height > interH) { 
           height = interH;
//           cout << "coming in inter" << endl;
       }
*/
       double* DistCopy = new double[3];
       DistCopy[0] = 0.0; DistCopy[1] = height; DistCopy[2] = DispThick;
//       cout << DistCopy[1] << endl;
//       EN_deleteData((pEntity)vertex, wallDistID);
       EN_modifyDataPtr((pEntity)vertex, wallDistID, (void *)DistCopy);
       PList_delete(VGrowth);
      } //ifOriginatingNode
     } //while VIter

     VIter_delete(vit); 
#ifdef DEBUG     
     M_writeVTKFile(mesh, "SeparatedBL", SeparatedID, 1);
//     M_writeVTKFile(mesh, "DistWall", wallDistID, 3);
#endif     
     MD_deleteMeshDataId(SeparatedID);
//     delete [] Sep;

//       printf("height:%lf\n",height);
//       return height;
// check if the vertex lies on model boundary and if owned by the proc    
}

void DistOffWall(pMesh mesh, pMSAdapt simAdapter){

   pVertex vertex;
   pVertex v;
   pPList VGrowth;
   
   int nL;
   int nMaxLayers = 50;
   double GrowthFactor = 1.25;
   double rho = rhoinp;
   double mu = muinp;
   double nu;
   double yplus = 1.0;
   nu = mu/rho;

   VIter vit=M_vertexIter(mesh);
   while(vertex = VIter_next(vit)){
      double WallDist, WallStress; 
      double height, uTau;
      double* Dist;
      double* Stress;
      height = 0.0;
      EN_getDataPtr((pEntity)vertex, wallStressID,(void **)&Stress);  
      
      if(V_isOriginatingNode(vertex)) {

         EN_getDataPtr((pEntity)vertex, wallDistID, (void **)&Dist);
/* for shear stress from phasta
      uTau=0.0;   
      for(int i=0; i<6; i++){
         uTau = uTau + Tau[i]*Tau[i];
      }  //assembling the dot product of the wall shear stress vector
*/
      WallStress = Stress[0];
      uTau = sqrt(WallStress/rho);
     
      WallDist = yplus*nu/uTau; //nu*Y1plus*sqrt(rho)/uTau incompressible
/*
////For compressible flows:
      double *nodalData;
      EN_getDataPtr((pEntity)vertex, ybarID,(void**)&nodalData);

      rho = nodalData[3]/(287.05*nodalData[4]);
//      printf("rho: %lf\n",rho);
      uTau = sqrt(WallStress/rho);
      nu = mu/rho;    
      WallDist = yplus*nu/uTau; //nu*Y1plus*sqrt(rho)/uTau incompressible
*/
      if(uTau<1e-9)
         WallDist = 0.0;

//insert a limiting condition when uTau tends to be zero, since that would
//make WallDist fairly large and we can not accept that.
//calculation reqd to compare this distance to other 2 in plane directions
//WallDist value should ALWAYS be less than any other delta value
//to have acceptable aspect ratio.
      
      double* H = new double[2];
//      double power;
      H[0] = WallDist;
      H[1] = Dist[1];
/*      
      double powSum = 0.0;
      for(int i=1; i<nL; i++){
         power = pow(GrowthFactor, i-1);
         powSum = powSum + power;
      }
      H[1] = H[0]*powSum;
*/
// /*      
      double T, sum;
      sum = 0.0;
      nL = nMaxLayers;
      height = Dist[1];
      for(int i=0; i<nMaxLayers; i++) {
         sum = pow(GrowthFactor, i-1);
         T = WallDist*sum;
         if(T>height) {
            nL = i;
            i = nMaxLayers;
         }   
      }
// */      
      H[2] = Dist[2];
//      nL = 20;      
//      printf("number of layers:%d\n",nL);
      MSA_setBLNormalSize(simAdapter, vertex, 1, nL, H);
      
      EN_deleteData((pEntity)vertex, wallDistID); 
      EN_attachDataPtr((pEntity)vertex, wallDistID, (void*)H);
      } //isOriginatingNode
   }   //VIter
   VIter_delete(vit);
#ifdef DEBUG   
   M_writeVTKFile(mesh, "wallDist", wallDistID, 3);
#endif   
   MD_deleteMeshDataId(wallDistID);

}

double SeparatedBL(pVertex vertex) {

   double* Sep = new double[1];

   double* nodalData;
   double dir[3];
   double magV;
   double height = 0.0;
   double xyz1[3], xyz2[3];
   pVertex v;
   pPList VGrowth;

   VGrowth = V_growthCurveVertices(vertex);
   int Size = PList_size(VGrowth);
   
   V_coord(vertex, xyz1);

   v = (pVertex)PList_item(VGrowth, 1);
   EN_getDataPtr((pEntity)v, ybarID, (void**)&nodalData); 

   magV = sqrt(nodalData[0]*nodalData[0]+nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]);

  dir[0] = nodalData[0]/magV;
  dir[1] = nodalData[1]/magV;
  dir[2] = nodalData[2]/magV;
  
  int ibreak=0;
  for(int i=2; i<Size; i++) {
     
    if(ibreak==0) {
      double* nodalDataOther;
      double dirOther[3];
      double magVOther;

      double dot;
      v = (pVertex)PList_item(VGrowth, i);
      if(!EN_getDataPtr((pEntity)v, ybarID, (void**)&nodalDataOther)) {
         cout << "error in Separated BL" << endl;
      }

      magVOther = sqrt(nodalDataOther[0]*nodalDataOther[0]+nodalDataOther[1]*nodalDataOther[1]+nodalDataOther[2]*nodalDataOther[2]);

      dirOther[0] = nodalDataOther[0]/magVOther;
      dirOther[1] = nodalDataOther[1]/magVOther;
      dirOther[2] = nodalDataOther[2]/magVOther;

      dot = dotProd(dir, dirOther);

      V_coord(v, xyz2);
      
     if(dot < -0.1) { 
//       cout << "separation!" << endl;
        Sep[0] = 1.0; 
        EN_modifyDataPtr((pEntity)vertex, SeparatedID, (void *)Sep);
        height = 2.0*sqrt(dist(xyz1,xyz2));
        ibreak=1;
     }

//     if 
    } //ibreak
  }  //pPlist
  PList_delete(VGrowth);
  return height;
  
}

double GrowthCurveIntersect(pVertex vertex) {
    
    pVertex vL;
    pPList VGrowth;
    double xyz1[3], xyz2[3];
    double distBase, distTop, distGC1, distGC2;
    double h = 0.0;
    double hOld = 0.0;
    int numE;
    VGrowth = V_growthCurveVertices(vertex);
    int Size = PList_size(VGrowth);
    
    vL = (pVertex)PList_item(VGrowth, Size-1); 

    V_coord(vertex, xyz1);
    V_coord(vL, xyz2);

    numE = V_numEdges(vL);
    for(int i=0; i<numE; i++) {

       pEdge lEdge;
       pVertex vOther;

       lEdge = V_edge(vL, i);
       vOther = E_otherVertex(lEdge, vL);

       if(PList_inList(VGrowth, vOther) || (!EN_isBLEntity(lEdge))) {
         continue;
       } else {

          pPList VGrowthOther;
          double xyzOther1[3], xyzOther2[3];
          VGrowthOther = V_growthCurveVertices(vOther);
          int SizeOther = PList_size(VGrowthOther);

          pVertex vOtherBase;
          vOtherBase = (pVertex)PList_item(VGrowthOther, 0);

          V_coord(vOtherBase, xyzOther1);
          V_coord(vOther, xyzOther2);

          distBase = sqrt(dist(xyz1, xyzOther1));
          distTop = sqrt(dist(xyz2, xyzOther2));
          distGC1 = sqrt(dist(xyz1, xyz2));
          distGC2 = sqrt(dist(xyzOther1, xyzOther2));

          if(distBase > distTop) {
               h = distGC1*distBase/(distBase-distTop);   
//               cout << "intersection in " << h << endl;
               if(hOld > h) 
                   h = hOld;
          }


       } //if
    hOld = h;
    } //loop on edges
/*
    double HConstraint = BLConstraint(vertex, distGC1);
    if(HConstraint >= 0.0) {
       h = HConstraint;
    }
*/

  PList_delete(VGrowth);
  return h;
}

double BLHeightVel(pPList VGrowth) {

     pVertex vertex, v; 
     int Size = PList_size(VGrowth);
     double* nodalData;
     double* nodalOtherData;
     double node1[3], node2[3];
     double threshold = 0.0001;
     double VelGrad, Dely, height, VelGradBase;
     double VelGradOld = 0.0;
     double Speed, SpeedOther;
     
     vertex=(pVertex)PList_item(VGrowth, 0);
     EN_getDataPtr((pEntity)vertex, ybarID,(void**)&nodalData);
     V_coord(vertex, node1);
     Speed = sqrt(nodalData[0]*nodalData[0]+nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]);
     
     double hOld = 0.0;
     for(int i=1; i<Size; i++) {

        v = (pVertex)PList_item(VGrowth, i);
        EN_getDataPtr((pEntity)v, ybarID,(void**)&nodalOtherData);
        V_coord(v, node2);
        Dely = sqrt(dist(node1, node2));
        SpeedOther = sqrt(nodalOtherData[0]*nodalOtherData[0]+nodalOtherData[1]*nodalOtherData[1]+nodalOtherData[2]*nodalOtherData[2]);

       VelGrad = (SpeedOther-Speed)/Dely;

       if(i==1)
          VelGradBase = VelGrad;

        Speed = SpeedOther;
        node1[0] = node2[0]; node1[1] = node2[1]; node1[2] = node2[2];
       
         if(VelGrad <= threshold*VelGradBase) {                     
//interpolate between heightOld and height to get the actual height
           height = hOld + Dely*(VelGradBase*threshold-VelGrad)/(VelGrad-VelGradOld); 
//           printf("heights: %lf %lf %lf\n",hOld, height, h);
           i = Size;
         }
      VelGradOld = VelGrad;
      hOld =  hOld + Dely; 

     } //loop on PList

     if(VelGrad>threshold*VelGradBase) {
          height = hOld*VelGrad/(threshold*VelGradBase);

//      if(fabs(VelGradBase) < 1e-8)  
//         cout << "NAN!" << endl;
//        cout << "height of BL " << hOld << endl;

          if(height>2*hOld)
               height = 2*hOld;
          if(height<hOld/3)
               height = hOld/3;


     }

//     cout << "VelGrad " << VelGrad << endl;
     return height;
} 

double BLHeightVort(pVertex vertex, pPList VGrowth) {
     
     pVertex v; 
     double node1[3], node2[3];
     double uTau, baseVort, height, Vorti, Vortw, VortiOld, hOld, h;
     int Size = PList_size(VGrowth);
     double* nodalVort;

     double threshold = 0.0002;

     if(EN_getDataPtr((pEntity)vertex, nodalVorticityID,
             (void**)&nodalVort)==NULL) {
        cout<<"\nerror in DistOffWall: no data attached to vertex\n";
        exit(0);
      }

     baseVort = sqrt(nodalVort[0]*nodalVort[0]+nodalVort[1]*nodalVort[1]+nodalVort[2]*nodalVort[2]);     
//     baseVort = 10000000;
     V_coord(vertex, node1);
/*     
     v = (pVertex)PList_item(VGrowth, 1);
     V_coord(v, node2);

//calculate the inplane component of velocity for wall shear stress
     Vort[0]=nodalVort[0]; Vort[1]=nodalVort[1]; Vort[2]=nodalVort[2];
     //calculate wall normal vector
     diffVt(node2, node1, wvec);
     normVt(wvec, wnorm);
     Vortw = dotProd(wnorm, Vort);
     
     for(int j=0; j<3; j++) {
        Vortinpl[j] = Vort[j] - Vortw*wnorm[j];
     }
     
     //we are gonna use in plane velocity magnitude for wall shear stress      

     baseVort =sqrt(Vortinpl[0]*Vortinpl[0]+Vortinpl[1]*Vortinpl[1]+Vortinpl[2]*Vortinpl[2]);
*/     
//total height based on vorticity
      VortiOld = baseVort;
      hOld = 0.0;
      for(int i=1; i<Size; i++) {
         double* nodalOtherVort;
         v = (pVertex)PList_item(VGrowth, i);
         if(EN_getDataPtr((pEntity)v, nodalVorticityID,
                   (void**)&nodalOtherVort)==NULL) {
             cout<<"\nerror in DistOffWall: no data attached to vertex\n";
             exit(0);
         }
         
         Vorti = sqrt(nodalOtherVort[0]*nodalOtherVort[0]+nodalOtherVort[1]*nodalOtherVort[1]+nodalOtherVort[2]*nodalOtherVort[2]);         
         V_coord(v, node2);
/*         
//calculate the inplane component of velocity for wall shear stress
     Vort[0]=nodalOtherVort[0]; Vort[1]=nodalOtherVort[1]; Vort[2]=nodalOtherVort[2];
     //calculate wall normal vector
     diffVt(node2, node1, wvec);
     normVt(wvec, wnorm);
     Vortw = dotProd(wnorm, Vort);
     
     for(int j=0; j<3; j++) {
        Vortinpl[j] = Vort[j] - Vortw*wnorm[j];
     }
     
     //we are gonna use in plane velocity magnitude for wall shear stress      

     Vorti =sqrt(Vortinpl[0]*Vortinpl[0]+Vortinpl[1]*Vortinpl[1]+Vortinpl[2]*Vortinpl[2]);
*/
         if(Vorti<=threshold*baseVort) {                     
           h = sqrt(dist(node1, node2));
//interpolate between heightOld and height to get the actual height
           height = hOld + (h-hOld)*(baseVort*threshold-VortiOld)/(Vorti-VortiOld); 
//           printf("heights: %lf %lf %lf\n",hOld, height, h);
           i = Size;
         }
      VortiOld = Vorti;
      hOld =  sqrt(dist(node1, node2));   
      }
         
       if(Vorti>threshold*baseVort) {
          height = sqrt(dist(node1, node2))*Vorti/(threshold*baseVort);
          if(height>2*sqrt(dist(node1, node2)))
              height = 2*sqrt(dist(node1, node2));
       }
       if(height>sqrt(dist(node1, node2))*2)
          height = sqrt(dist(node1, node2)*2);

       if(height<sqrt(dist(node1, node2))/2)
          height = sqrt(dist(node1, node2)/2);

  if(baseVort<1e-3) height = 1e-5;
  return height;       
}

double BLConstraint(pVertex vert, double distGC) {
 
  double height = -1.0;
  double xyz[3];

  V_coord(vert, xyz);

// following for M6wing
 /*          
          if(xyz[1]>0.18 && xyz[0]>1.06) {
            height = distGC;
//            cout << distGC1 << endl;
          }
 */

 // following for Delery bump
 /* 
          if((0.295 > xyz[0] > 0.282)) {
//             h = distGC1;
             height = 0.004;
//             cout << distGC1 << endl;
          }
*/        
 
 //following for NACA0012
 /* 
   if(xyz[0] > 0.243 && xyz[0] < 0.249) {
     height = 0.01;
   }
   else if(xyz[0] > -0.168 && xyz[0] < -0.159) {
      height = 0.001;
   }
 */
  return height;  
}

void DetectShock(pParMesh pmesh, pMesh mesh) {
   
   VIter vIter;
   pVertex vertex;  

   int fieldIndexForGrad = 3;
   double* nodalGrad;
   double* dist;
   double* Shock = new double[1];
   double *distNew = new double[3];
   double* ShockModify = new double[1];
   Shock[0]=0.0;

   VortgradientsFromPatch(mesh, fieldIndexForGrad);
   commuGradientsFromPatch(pmesh, mesh);

   vIter = M_vertexIter(mesh);
   while(vertex = VIter_next(vIter)) {

      EN_deleteData((pEntity)vertex, localGradientID);
      EN_deleteData((pEntity)vertex, localPatchVolID);
      
      EN_attachDataPtr((pEntity)vertex, ShockID, (void *)Shock);

      if(V_isOriginatingNode(vertex)) {

         EN_getDataPtr((pEntity)vertex, nodalGradientID,(void**)&nodalGrad);

         if(fabs(nodalGrad[0]) > 1.0e6) {
//            cout << "found shock" << endl;
            EN_getDataPtr((pEntity)vertex, wallDistID,(void**)&dist);
            
            dist[1] = 0.0002; 
            distNew[0] = dist[0]; distNew[1] = dist[1]; distNew[2] = dist[2];
            EN_modifyDataPtr((pEntity)vertex, wallDistID, (void *)distNew);
            
            ShockModify[0] = 1.0;
            EN_modifyDataPtr((pEntity)vertex, ShockID, (void *)ShockModify);
//            delete [] ShockModify;
         }   
      }

      EN_deleteData((pEntity)vertex, nodalGradientID);

   }
   VIter_delete(vIter);
#ifdef DEBUG 
//  M_writeVTKFile(mesh, "PresGrad", nodalGradientID, 3);
//  M_writeVTKFile(mesh, "DetectShock", ShockID, 1);
#endif   
  MD_deleteMeshDataId(ShockID);
}

int V_isOriginatingNode(pVertex vertex) {
    int isOrg;
    EN_getDataInt((pEntity)vertex, isOrgNodeID, &isOrg);
    return isOrg;
}

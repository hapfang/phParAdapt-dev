#include <cstdio>
#include "phParAdapt.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#ifndef ibm
#define ludcmp ludcmp_
#define lubksb lubksb_
#endif 
#include "Eigen.h"
#ifdef SIM
#include "MeshSim.h"
#include "MeshTypes.h"
#include "SimMeshTools.h"
#endif
#include <math.h>
#ifdef FMDB
  #include "mEntity.h"
#endif

extern "C" double XYZ_volume(dArray*);
extern int adaptOption;
extern double rhoinp;

using namespace std;

extern pMeshDataId localPatchVolID;
extern pMeshDataId localGradientID;
extern pMeshDataId ybarID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalVorticityID;
// recover gradient at a vertex using patch of elements around it
//  attaches:
//  Sum of weighted (local) gradients to vertices via localGradientID
//  Sum of (local) patch Volumes to vertices via localPatchVolID
// this is a LOCAL operation

void Vorticity(pParMesh pmesh, pMesh mesh) {
   
  VIter vIter;
  vIter = M_vertexIter(mesh);
  pVertex vertex;

  int nshg = M_numVertices(mesh);
//  double velGrad[3][nshg*3];
  double **velGrad;
  velGrad = new double*[3];
  double* nodalGrad;
  int vCount; 
//Calculate gradients of u, v and w for vorticity calculation
   for(int i=0; i<3; i++) { 
      velGrad[i] = new double[nshg*3]; 
      VortgradientsFromPatch(mesh, i);
      commuGradientsFromPatch(pmesh, mesh);
      
      vCount = 0;
      while(vertex = VIter_next(vIter)) {
         EN_deleteData((pEntity)vertex, localGradientID);
         EN_deleteData((pEntity)vertex, localPatchVolID);

         if(EN_getDataPtr((pEntity)vertex, nodalGradientID,
                       (void**)&nodalGrad)==NULL){
           cout<<"\nerror in Vorticty: no data attached to vertex\n";
           exit(0);
         }
         
         for(int j=0; j<3; j++) {
            velGrad[i][vCount*3+j] = nodalGrad[j];
         }
         
         delete [] nodalGrad;
         EN_deleteData((pEntity)vertex, nodalGradientID);
         vCount++;
      }
   VIter_reset(vIter);
   }
      
   vCount = 0;

   while(vertex = VIter_next(vIter)) {
      double* vorticity = new double[3];
      vorticity[0] = velGrad[2][vCount*3+1] - velGrad[1][vCount*3+2];
      vorticity[1] = velGrad[0][vCount*3+2] - velGrad[2][vCount*3];
      vorticity[2] = velGrad[1][vCount*3] - velGrad[0][vCount*3+1];

      EN_attachDataPtr((pEntity)vertex, nodalVorticityID, (void *)vorticity);
      vCount++;
   }

   delete [] velGrad;
   VIter_delete(vIter);
#ifdef DEBUG   
   M_writeVTKFile(mesh, "Vorticity", nodalVorticityID, 3);
#endif   
}


void
VortgradientsFromPatch(pMesh mesh, int fieldIndexForGrad)
{
  VIter vIter;
  vIter = M_vertexIter(mesh);
  pVertex vertex;
  pRegion region;

  double elemGradient[3];
  int count = 0;
  
  // loop over vertices and get patch of elements
  while(vertex = VIter_next(vIter)) {
      pPList elementPatch;
      elementPatch = V_regions(vertex);
      
      double* patchVolume = new double[1];
      patchVolume[0] = 0;
      
      double* patchGradient = new double[3];
      for(int i=0; i<3; i++){
          patchGradient[i]= 0.0  ;
      }

#if  ( defined  DEBUG ) 
//       printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
//       printf("[%i]:in gradientsFromPatch: vertex:",PMU_rank());
//       double c[3];
//       V_coord(vertex,c);
//       printf("coords %f %f %f\n",c[0],c[1],c[2]);
//       printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif      

      // loop over elements and recover gradient at vertex   
      for (int i=0; i<PList_size(elementPatch); i++) {

	region = (pRegion) PList_item(elementPatch,i);
#ifdef SIM   
	double vol = R_volume(region);
#endif   
#ifdef FMDB   
	double vol = R_Volume2(region);
#endif   
	*patchVolume+=vol;

        // retrieve element grad. for each element in the patch at vertex
        VortelementGradient(region, elemGradient, fieldIndexForGrad);

#if  ( defined  DEBUG )
//           printf("\n[%i]: ELEM grad :\n %f\n %f\n %f\n",PMU_rank(),elemGradient[0],elemGradient[1],elemGradient[2]);
//           pPList regionVerts;
//           regionVerts = R_vertices( (pRegion)PList_item(elementPatch,i),0);
//           double co[3];
//           printf("\nand region for vertex %i:\n",count);
//           for(int k=0;k<4;k++){
//               V_coord((pVertex)PList_item(regionVerts,k),co);
//               printf(" %f %f %f\n",co[0],co[1],co[2]);
//           }
//           PList_delete(regionVerts);   
#endif      


          // use element grad. (each component) to recover gradient at vertex
          for(int j=0; j<3; j++) {
              patchGradient[j]= patchGradient[j] + 
                  elemGradient[j]*vol;
          }
      }//loop over element patch ends
      
      
#if  ( defined  DEBUG )
//       printf("\n[%i]:in gradientsFromPatch  local vals:\n %f\n %f\n %f\n",PMU_rank(),patchGradient[0],patchGradient[1],patchGradient[2]);
//       double ci[3];
//       printf("\nand local vol %f:\n",patchVolume[0]);
//       V_coord(vertex,ci);
//       printf(" for vertex %i   %f %f %f\n",count,ci[0],ci[1],ci[2]);
#endif       

      // attach the recovered gradient to the vertex
      EN_attachDataPtr( (pEntity)vertex, localGradientID, (void *)
                        patchGradient );
      
      EN_attachDataPtr( (pEntity)vertex, localPatchVolID, (void *)
                        patchVolume);
      
      PList_delete(elementPatch);
      count++;
  }
  VIter_delete(vIter);

}

// reconstruct the element gradient
void
VortelementGradient(pRegion region, double* elementGradient, int fieldIndexForGrad)
{
  pVertex v;
  double* nodalData;
  
  int four = 4;
  for(int i=0; i<3; i++) {
    elementGradient[i] = 0.;
  }

  double matrix[16];
  double fieldVector[4];

  pPList regionVerts;
  // regionVerts are deleted at the end of the routine
  regionVerts = R_vertices(region, 1);

#ifdef FMDB
  int rTopoType = region->getType();
  switch(rTopoType) {
  case mEntity::TET : {
#endif
#ifdef SIM
  int rTopoType = R_topoType(region);
  switch(rTopoType) {
  case Rtet : {
#endif
    // build the linear system
    buildSystem(region,matrix);

    // get the field vals
    for (int i=0; i<PList_size(regionVerts); i++) {
      v = (pVertex)PList_item(regionVerts,i);
      if(EN_getDataPtr((pEntity)v, ybarID,
                       (void**)&nodalData)==NULL){
        cout<<"\nerror in elementGradient: no data attached to vertex (tet element)\n";
        exit(0);
      }
    fieldVector[i]=nodalData[fieldIndexForGrad];

    }

    int indx[4];
    double fnumber;
    ludcmp(matrix , &four, &four, indx, &fnumber);

    // fieldVector is going to be overridden and will
    // contain the solution
    lubksb(matrix , &four, &four, indx, fieldVector );

    for(int i=0; i<3; i++) {
      // the poly's structure is
      // a0 + a1*x + a2*y + a3*z
      elementGradient[i] = fieldVector[i+1];

#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;

#ifdef FMDB
  case mEntity::PRISM : {
    double volume;
    volume = R_Volume2(region);
#endif
#ifdef SIM
  case Rwedge : {
    double volume;
#endif
    double volTet, x_y_z[4][3];
    // fit three tets in a prism
    // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
    int tetVerts[3][4] = {0,1,2,5, 0,1,5,4, 0,4,5,3};

    double v01[3],v02[3],v03[3];
    double normal[3];

    // three tets in a prism
    // compute gradient on prism as vol. averaged gradients on three child tets
    for(int iTet=0; iTet<3; iTet++) {
      for(int iVert=0; iVert<4; iVert++) {
        v = (pVertex)PList_item(regionVerts,tetVerts[iTet][iVert]);

        if(EN_getDataPtr((pEntity)v, ybarID,
                         (void**)&nodalData)==NULL) {
          cout<<"\nerror in elementGradient: no data attached to vertex (wedge element)\n";
          exit(0);
        }
        fieldVector[iVert]=nodalData[fieldIndexForGrad];
        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;

#ifdef FMDB
//       #ifdef DEBUG
      // volume of wedge should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped (and if base tri. is thin)
      // --- probably can just check the volume of the wedge/prism
      if(volume<0.&&volTet<0.) {
        cout<<"\nError in elementGradient()..."<<endl;
        cout<<"negative volume of tet ["<<iTet<<"] in prism/wedge of negative volume"<<endl;

        R_info(region);
        exit(0);
      }
//       #endif
#endif

      // volume of wedge should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped (and if base tri. is thin)
      volTet = fabs(volTet);

      // build the linear system
      buildSystemXYZ(x_y_z,matrix);

      int indx[4];
      double fnumber;
      ludcmp(matrix , &four, &four, indx, &fnumber);

      // fieldVector is going to be overridden and will
      // contain the solution
      lubksb(matrix , &four, &four, indx, fieldVector );

      for(int i=0; i<3; i++) {
        // the poly's structure is
        // a0 + a1*x + a2*y + a3*z
//         elementGradient[i] += volTet*fieldVector[i+1];
        elementGradient[i] += fieldVector[i+1];
        if(isnan(fieldVector[i])) {
            printf("ISNAN in elementGrad prism %d\n",i);                       
        }
      }

    } // loop over three tets fitted into prism

    // gradient for prism/wedge element is
    // vol. weighted average of gradients on three fitted child tets
    for(int i=0; i<3; i++) {
      // elementGradient[i] = elementGradient[i]/volume;
      elementGradient[i] = elementGradient[i]/3;
#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;

#ifdef FMDB
  case mEntity::PYRAMID : {
    double volume;
    volume = R_Volume2(region);
#endif
#ifdef SIM
  case Rpyramid : {
    double volume;
#endif
    double volTet, x_y_z[4][3];
    // fit two tets in a pyramid
    // (diagonal edge of a pyramid is : V0-V2)
    int tetVerts[2][4] = {{0,1,2,4}, {0,2,3,4}};

    double v01[3],v02[3],v03[3];
    double normal[3];

    // two tets in a pyramid
    // compute gradient on pyramid as vol. averaged gradients on two child tets
    for(int iTet=0; iTet<2; iTet++) {
      for(int iVert=0; iVert<4; iVert++) {
        v = (pVertex)PList_item(regionVerts,tetVerts[iTet][iVert]);

        if(EN_getDataPtr((pEntity)v, ybarID,
                         (void**)&nodalData)==NULL) {
          cout<<"\nerror in elementGradient: no data attached to vertex (wedge element)\n";
          exit(0);
        }
        fieldVector[iVert]=nodalData[fieldIndexForGrad];
        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;

#ifdef FMDB
//       #ifdef DEBUG
      // volume of pyramid should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped
      // --- probably can just check the volume of the pyramid
      if(volume<0.&&volTet<0.) {
        cout<<"\nError in elementGradient()..."<<endl;
        cout<<"negative volume of tet ["<<iTet<<"] in pyramid of negative volume"<<endl;

        R_info(region);
        exit(0);
      }
//       #endif
#endif

      // volume of pyramid should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped
      volTet = fabs(volTet);

      // build the linear system
      buildSystemXYZ(x_y_z,matrix);

      int indx[4];
      double fnumber;
      ludcmp(matrix , &four, &four, indx, &fnumber);

      // fieldVector is going to be overridden and will
      // contain the solution
      lubksb(matrix , &four, &four, indx, fieldVector );

      for(int i=0; i<3; i++) {
        // the poly's structure is
        // a0 + a1*x + a2*y + a3*z
        // elementGradient[i] += volTet*fieldVector[i+1];
        elementGradient[i] += fieldVector[i+1];
      }

    } // loop over two tets fitted into prism

    // gradient for pyramid element is
    // vol. weighted average of gradients on two fitted child tets
    for(int i=0; i<3; i++) {
      // elementGradient[i] = elementGradient[i]/volume;
      elementGradient[i] = elementGradient[i]/2;
#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;
  default: {
    cout<<"\nError in elementGradient()..."<<endl;
    cout<<"element topology ["<<rTopoType<<"] CANNOT be handelled"<<endl;
    exit(0);
  }
  }

  PList_delete(regionVerts);
}


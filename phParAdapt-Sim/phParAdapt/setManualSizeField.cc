#include "phParAdapt.h"
#include "Eigen.h"
#include <iostream>
#include <fstream>
#include "attachData.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;
extern pMeshDataId nodalSizeID;
extern pMeshDataId OrgSizeID;
extern pMeshDataId nodalDirectionID;

extern pProgress prog;
void 
setManualSizeField( pParMesh pmesh, 
         pMesh mesh,
		   pMSAdapt simAdapter, 
     	   int strategy,
         int option   ) {
  double dir[3][3];
  char option2[28];
  ofstream fout;
  
  switch(strategy) {
  case -1:
    {
      sprintf(option2,"constant");

      dir[0][0]=1.0;//x-direction
      dir[1][1]=1.0;//y-direction
      dir[2][2]=1.0;//z-direction
    
      dir[0][1]=dir[0][2]=0.;
      dir[1][0]=dir[1][2]=0.;
      dir[2][0]=dir[2][1]=0.;

      pVertex vertex;
      VIter vit=M_vertexIter(mesh);
      while(vertex=VIter_next(vit)) {  
	MSA_setAnisoVertexSize(simAdapter, 
			       vertex,
			       dir);
      }
      VIter_delete(vit);
    }
    break;
  case -2:
    {
      sprintf(option2,"cylindrical");

      double norm,tol=1.e-8;
    
      double sizeR=0.001;
      double sizeTheta=0.001;
      double sizeZ=0.01;

      pVertex vertex;
      VIter vit=M_vertexIter(mesh);
      while(vertex=VIter_next(vit)) {  
	double xyz[3];
	V_coord(vertex,xyz);
	norm=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	if( norm>tol ) 
	  {
	    dir[1][1]=sizeR*(xyz[1]/norm);
	    dir[1][2]=sizeR*(xyz[2]/norm);
	    dir[1][0]=0.;
	    dir[2][1]=-sizeTheta*(xyz[2]/norm);
	    dir[2][2]=sizeTheta*(xyz[1]/norm);
	    dir[2][0]=0.;
	    dir[0][0]=0.;
	    dir[0][1]=0.;
	    dir[0][0]=sizeZ;
	  }
	else
	  {
	    dir[1][1]=sizeR;
	    dir[0][1]=0.;
	    dir[0][2]=0.;
	    dir[1][0]=0.;
	    dir[2][2]=sizeTheta;
	    dir[1][2]=0.;
	    dir[2][0]=0.;
	    dir[2][1]=0.;
	    dir[0][0]=sizeZ;
	  }

	MSA_setAnisoVertexSize(simAdapter, 
			       vertex,
			       dir);
    double* dir2 = new double[9];
    dir2[0]=dir[0][0];dir2[1]=dir[0][1];dir2[2]=dir[0][2];
    dir2[3]=dir[1][0];dir2[4]=dir[1][1];dir2[5]=dir[1][2];
    dir2[6]=dir[2][0];dir2[7]=dir[2][1];dir2[8]=dir[2][2];

    EN_attachDataPtr( (pEntity)vertex, nodalDirectionID, (void *)
                        dir2);
 
// /* compute the size field
        double* h2 = new double[3];
        for(int k=0; k<3; k++) {
           h2[k]=0.0;
        }
        for (int k=0;k<3;k++){
          for (int j=0 ;j<3; j++) {
            h2[k] +=(dir[k][j])*(dir[k][j]);
          }
          h2[k] = sqrt(h2[k]);
        }
        EN_attachDataPtr( (pEntity)vertex, nodalSizeID, (void *)
                        h2);
      }
      VIter_delete(vit);
    }
#ifdef DEBUG    
//    M_writeVTKFile(mesh, "nodalSize", nodalSizeID, 3);
//    M_writeVTKFile(mesh, "nodalDirection",nodalDirectionID, 9);
#endif    
    break;
  case -3:
    {
    
    char sizefile[256];
    sprintf(sizefile,"OrgSize%d.dat", PMU_rank());
    fout.open(sizefile);

    int nshg = M_numVertices(mesh);
    fout<< nshg <<"\n"; 
    pVertex vertex;
    VIter vit=M_vertexIter(mesh);
    int vCount = 0;
    while(vertex=VIter_next(vit)) { 
       
      double* OrgSize = new double;
      double newSize;
      double OrgAnisoSize[3][3];
      int iSize;

      double* ybar;
      if(option==-1) {
         if(!EN_getDataPtr((pEntity)vertex, errorIndicatorID, (void **)&ybar)) {
               printf("error in setManualSize, no data attached\n");
               exit(0);
         }
      }
      
      vCount++;
      iSize = V_size(vertex, OrgSize, OrgAnisoSize);
      double* Orgh = new double[3];
      for(int k=0; k<3; k++) {
          Orgh[k]=0.0;
      }
      
      if (iSize == 1){
         for(int k=0; k<3; k++) {
            Orgh[k]= *OrgSize;
         }  
            dir[0][1]=dir[0][2]=0.;
            dir[1][0]=dir[1][2]=0.;
            dir[2][1]=dir[2][0]=0.;
//          printf("ybar: %lf\n", ybar[5]);
//         if(option==-1 && (ybar[5] > 002 || ybar[5] < -0.01)) {
//            dir[0][0] = dir[1][1] = dir[2][2] = *OrgSize/4;
         
//            newSize = *OrgSize/4;
//         } else {

            dir[0][0] = dir[1][1] = dir[2][2] = *OrgSize/2;
            newSize = *OrgSize;
         }
     
      if (iSize == 2){
         for (int k=0;k<3;k++){
            for (int j=0 ;j<3; j++) {
               Orgh[k] +=(OrgAnisoSize[k][j])*(OrgAnisoSize[k][j]);
            }
            Orgh[k] = sqrt(Orgh[k]);

//            if(option==-1 && (ybar[5] > 0.02 || ybar[5] < -0.01)) {
//               for (int l=0;l<3; l++){
//                 dir[k][l] = OrgAnisoSize[k][l]/4;
//               }   
//            } else {
               for (int l=0;l<3; l++){
                  dir[k][l] = OrgAnisoSize[k][l]/2;
               }
         }
         newSize = Orgh[0]/2;
      }
      
//         MSA_setAnisoVertexSize(simAdapter, vertex, dir);
       MSA_setVertexSize(simAdapter, vertex, newSize);

         fout<<Orgh[0]<<" "<<Orgh[1]<<" "<<Orgh[2]<<"\n";
//         fout<<dir[0][0]<<" "<<dir[0][1]<<" "<<dir[0][2]<<" ";
//         fout<<dir[1][0]<<" "<<dir[1][1]<<" "<<dir[1][2]<<" ";
//         fout<<dir[2][0]<<" "<<dir[2][1]<<" "<<dir[2][2]<<"\n";
         
// /* compute the size field
         double* h2 = new double[3];
         for(int k=0; k<3; k++) {
            h2[k]=0.0;
         }
         for (int k=0;k<3;k++){
            for (int j=0 ;j<3; j++) {
               h2[k] +=(dir[k][j])*(dir[k][j]);
            }
            h2[k] = sqrt(h2[k]);
         }

//         double* newSize2 = new double;
//         *newSize2 = newSize;

         EN_attachDataPtr( (pEntity)vertex, nodalSizeID, (void *)
               h2);

         EN_attachDataPtr( (pEntity)vertex, OrgSizeID, (void *)
              Orgh);    
    }
    VIter_delete(vit);
    fout.close();
#ifdef DEBUG    
//    M_writeVTKFile(mesh, "nodalSize", nodalSizeID, 3);
//    M_writeVTKFile(mesh, "OrgSize", OrgSizeID, 3);
//    M_write(mesh, "mesh_size.sms", prog);
#endif    
    }
    break;   
  default :
    cout<<"check strategy [in setManualSizefield(...)]"<<endl;
    exit(-1);
    break; 
  }

  ofstream adaptSimLog("phAdapt.log");
  adaptSimLog<<"Strategy chosen for adaptation is size-field driven"<<endl;
  adaptSimLog<<"Mesh size-field is set manually"<<endl;
  adaptSimLog<<"Size-field option : "<<option<<endl;
  adaptSimLog.close();
}

#ifdef __cplusplus
}
#endif

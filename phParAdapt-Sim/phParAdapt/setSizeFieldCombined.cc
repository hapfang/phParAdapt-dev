#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef SIM
#include "SimMeshTools.h"
#endif
#include "MeshSimAdapt.h"

#include "phParAdapt.h"
#include "Eigen.h"
#include "attachData.h"
#include "phReadWrite.h"

#include "mpi.h"


#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;
extern pMeshDataId nodalSizeID;
extern pMeshDataId nodalDirectionID;
extern pMeshDataId oldMeshSizeID;
extern pMeshDataId OrgSizeID;
extern pMeshDataId QualitySizeID;

extern pMeshDataId localGradientID; 
extern pMeshDataId localPatchVolID;
extern pMeshDataId localHessianID;
extern pMeshDataId numSurroundNodesID;
extern pMeshDataId locMaxInterpolErrorID;
extern pMeshDataId globMaxInterpolErrorID;

extern int nErrorVars;
extern int timeStepNumber;
extern int isSizeLimit;

extern pMeshDataId meshSizeID;

extern int preLBforAdaptivity;

void setSizeFieldUsingCombined(pParMesh pmesh,
                               pMesh mesh,
			       pMSAdapt simAdapter,
			       double factor,
			       double hmax,
			       double hmin,
			       int option)
{

  meshSizeID = MD_newMeshDataId("mesh size ID");
  oldMeshSizeID = MD_newMeshDataId("isotropic old mesh size");

  int nshg=M_numVertices(mesh);  // only true for linear elements
  double tol=1.e-12;
  double totalError = 0.;
  double threshold  = 0.;
  double sumOfError = 0.;
  
  double T[3][3];
  double maxEigenval=0;
  int i,j;
  
  pVertex vertex;
  double eloc;  	  // local error at a vertex
  double etot=0.;	  // total error for all vertices
  double etotLoc=0.;	  // total error for all vertices local to partition
  double emean; 	  // emean = etot / nv
  double elocmax=0.;	  // max local error
  double elocmin=1.e20;   // min local error

  pProgress prog;

  char adaptSimLogFile[128];
  sprintf(adaptSimLogFile,"phAdapt.%i.log",PMU_gid(PMU_rank(),0)+1);
  // ofstream adaptSimLog(adaptSimLogFile);
  FILE *adaptSimLog=fopen(adaptSimLogFile,"w");

  // to get an idea refer Appendix A in Li's thesis
  // attach the max local interpolation error LOCAL
  // to each partition via locMaxInterpolErrorID;
  maxLocalPartLocError(mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // attach the max local interpolation error GLOBAL
  // to the parallel mesh via globMaxInterpolErrorID;
  commuMaxLocalPartLocError(pmesh,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // struct Hessian contains decomposed values
  // mesh sizes and directional information
  Hessian *hess = new Hessian[nshg];

  if (option == 4) {

     getGlobalErrorInfo(mesh,totalError,sumOfError);

     computeOldMeshSize(mesh,option);
     MPI_Barrier(MPI_COMM_WORLD);

     commuOldMeshSize(pmesh,mesh);
     MPI_Barrier(MPI_COMM_WORLD);
  }

  VIter vit=M_vertexIter(mesh);
  i=0;
  while(vertex=VIter_next(vit)) { 
    double xyz[3];
    V_coord(vertex,xyz);

    // get hessian of appropriate variable at a `vertex'
    if(option==1 || option==2 || option==9 || option == 4) 
      V_Hessian(vertex,T);
    else{

//      V_AnalyticHessian(vertex,T,option);
        printf("in setSizeFieldUsingHessians: only option=1 works at the moment \n");
        exit(0);
    } 

    // `T' is the actual Hessian 
    // of any solution variable (e.g., T, u, v, w etc.)
    // `hess' has the information related to Mesh-Metric (eventually)
    // used for mapping the ellipsoid into a unit sphere 
    
    // compute eigen values and eigen vectors
    int nbEigen = eigen(T,hess[i].dir,hess[i].h);

    if(nbEigen == 1) {
      printf(" [%d] ",PMU_rank());
      printf("WARNING : nbEigen is 1\n");
      printf(" Vertex Number : %d\n",i);
      printf(" xyz : %f %f %f \n",xyz[0],xyz[1],xyz[2]);
      printf("       %f %f %f %f %f %f\n\n",T[0][0], 
             T[0][1],T[0][2],T[1][1],T[1][2],T[2][2]);
    }

    hess[i].h[0] = ABS(hess[i].h[0]);
    hess[i].h[1] = ABS(hess[i].h[1]);
    hess[i].h[2] = ABS(hess[i].h[2]);
    
    if( hess[i].h[0] < tol ) {
      printf(" [%d] ",PMU_rank());
      printf("Error: zero maximum eigenvalue !!!\n");
      printf(" xyz : %f %f %f \n",xyz[0],xyz[1],xyz[2]);
      printf("       %f %f %f\n", hess[i].h[0],
             hess[i].h[1],hess[i].h[2]);
      printf("       %f %f %f %f %f %f\n\n",T[0][0], 
             T[0][1],T[0][2],T[1][1],T[1][2],T[2][2]);

//       adaptSimLog<<" ["<<PMU_rank()<<"]"<<endl;
//       adaptSimLog<<"Error: zero maximum eigenvalue !!!"<<endl;
//       adaptSimLog<<" xyz : "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n"<<endl;
//       adaptSimLog<<hess[i].h[0]<<" "<<hess[i].h[1]<<" "<<hess[i].h[2]<<endl;
//       adaptSimLog<<T[0][0]<<" "<<T[0][1]<<" "<<T[0][2]<<" ";
//       adaptSimLog<<T[1][1]<<" "<<T[1][2]<<" "<<T[2][2]<<endl;
     printf("For the above point, setting unit vectors as directions and hmax as sizes. If the point is not farfield, this is a serious error!\n\n");
     hess[i].dir[0][0]=hess[i].dir[1][1]=hess[i].dir[2][2]=1.0;
     hess[i].dir[0][1]=hess[i].dir[0][2]=hess[i].dir[1][0]= 0.0;
     hess[i].dir[1][2]=hess[i].dir[2][0]=hess[i].dir[2][1]= 0.0;
//      continue;
    }

    // estimate relative interpolation error
    // needed for scaling metric field (mesh size field)

    // retrieve the local interpolation error valid globally via 
    // globMaxInterpolErrorID
    double* globMaxInterpolError;
    if(!EN_getDataPtr((pEntity)vertex,globMaxInterpolErrorID ,
		      (void**)&globMaxInterpolError)){
      printf("\nerror in setSizeFieldUsingHEssians: no data attached to vertex\n");
                    exit(0);
    }
    eloc=globMaxInterpolError[0];


    // summing up contribs over verts
    // ONLY if vertex is OWNED
    // otherwise we would have double contributions
    if(EN_isOwnerProc((pEntity)vertex)){
        etotLoc += eloc;
    }
    if( eloc>elocmax )  elocmax=eloc;
    if( eloc<elocmin )  elocmin=eloc;

    ++i;
  }//while(vertex
  printf("Info [%d]: Reading hessian... done...\n\n",PMU_rank());

  MPI_Barrier(MPI_COMM_WORLD);
  // communicate the total
  MPI_Allreduce(&etotLoc, &etot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  // get the total number of verts, overal procs/parts
  int nshgTOT = getNSHGTOT(mesh);
  // cout<<"\n["<<PMU_rank()<<"] nshgTOT : "<<nshgTOT<<"\n"<<endl;
  
  emean =  etot / nshgTOT;

  if(PMU_rank()==0) {
    printf("\n Info on relative interpolation error: ");
    printf("   total: %f\n",etot);
    printf("   mean:  %f\n",emean);
    printf("   max local: %f\n",elocmax);
    printf("   min local: %f\n",elocmin);
  }

  eloc = emean*factor;

  if(PMU_rank()==0) {
    printf("\n towards uniform local error distribution of %f\n", eloc);
    printf("   with max edge length=%f; min edge length=%f\n\n",hmax,hmin);
  }

//   adaptSimLog<<"MeshSim   Library build ID : "<<SimMeshing_buildID()<<endl;
//   adaptSimLog<<"PMesh     Library build ID : "<<SimPMesh_buildID()<<endl;
//   adaptSimLog<<"MeshTools Library build ID : "<<SimMeshTools_buildID()<<endl;
//   adaptSimLog<<endl;
  fprintf(adaptSimLog,"MeshSim   Library build ID : %s\n",SimMeshing_buildID());
  fprintf(adaptSimLog,"PMesh     Library build ID : %s\n",SimPartitionedMesh_buildID());
  fprintf(adaptSimLog,"MeshTools Library build ID : %s\n",SimMeshTools_buildID());
  fprintf(adaptSimLog,"\n");

  if(option==1 || option==2) {
    // adaptSimLog<<"Strategy chosen is anisotropic adaptation, i.e., size-field driven"<<endl;
    fprintf(adaptSimLog,"Strategy chosen is anisotropic adaptation, i.e., size-field driven\n");
  }
  if(option==9) {
    // adaptSimLog<<"Strategy chosen is isotropic adaptation (based on local hmin of mesh metric tensor from hessians), i.e., size-field driven"<<endl;
    fprintf(adaptSimLog,"Strategy chosen is isotropic adaptation (based on local hmin of mesh metric tensor from hessians), i.e., size-field driven\n");
  }
  if(option==1 || option==2 || option==9) {
    // adaptSimLog<<"Using numerical/computed hessian... done..."<<endl;
    fprintf(adaptSimLog,"Using numerical/computed hessian... done...\n");

  }
  else {
    // adaptSimLog<<"Using analytic hessian... done..."<<endl;
    fprintf(adaptSimLog,"Using analytic hessian... done...\n");
  }

//   adaptSimLog<<"Info on relative interpolation error :"<<endl;
//   adaptSimLog<<"total : "<<etot<<endl;
//   adaptSimLog<<"mean : "<<emean<<endl;
//   adaptSimLog<<"factor : "<<factor<<endl;
//   adaptSimLog<<"min. local : "<<elocmin<<endl;
//   adaptSimLog<<"max. local : "<<elocmax<<endl;
//   adaptSimLog<<"towards uniform local error distribution of "<<eloc<<endl;
//   adaptSimLog<<"with min. edge length : "<<hmin<<endl;
//   adaptSimLog<<"with max. edge length : "<<hmax<<endl;

  fprintf(adaptSimLog,"Info on relative interpolation error :\n");
  fprintf(adaptSimLog,"total : %f\n",etot);
  fprintf(adaptSimLog,"mean : %f\n",emean);
  fprintf(adaptSimLog,"factor : %f\n",factor);
  fprintf(adaptSimLog,"min. local : %f\n",elocmin);
  fprintf(adaptSimLog,"max. local : %f\n",elocmax);
  fprintf(adaptSimLog,"towards uniform local error distribution of %f\n",eloc);
  fprintf(adaptSimLog,"with min. edge length : %f\n",hmin);
  fprintf(adaptSimLog,"with max. edge length : %f\n",hmax);

  // scale the hessian
  // double errorSqRoot = hmin*sqrt(maxEigenval);
  VIter_reset(vit);
  i=0;

  int foundHmin = 0;
  int foundHmax = 0;
  int hminCount = 0;
  int hmaxCount = 0;
  int bothHminHmaxCount = 0;
  double tol2 = 0.01*hmax, tol3 = 0.01*hmin;

  while ( vertex=VIter_next(vit)) {

    foundHmin = 0;
    foundHmax = 0;

    for( j=0; j<3; j++ ) {
      if( hess[i].h[j] < tol )
        hess[i].h[j] = hmax;
      else {
        hess[i].h[j] = sqrt(eloc/hess[i].h[j]);
        if( hess[i].h[j] > hmax )
          hess[i].h[j] = hmax;
        if( hess[i].h[j] < hmin )
          hess[i].h[j] = hmin;
      }
    }


    if(option == 4) { //combined strategy

        double *nodalValue;
        if(!EN_getDataPtr((pEntity)vertex,errorIndicatorID,(void**)&nodalValue)) {
           cout<<"\nerror in setIsotropicSizeField(...) : no error data attached to vertex\n";
           exit(0);
         }

        threshold = factor*totalError;
        // based on mesh optimal criteria and error convergence rate
        // implemented for linear cases
        // refer L.-Y. Li, Comm. Num. Meth. Eng.,
        // Vol. 11, 857-868 & 911-915 (1995)
        double adaptFactorIso = threshold/(pow(*nodalValue,0.25)*sqrt(sumOfError));

        double MaxRefineFactor = hmin;
        double MaxCoarsenFactor = hmax;

        double *oldSize;
        if(!EN_getDataPtr((pEntity)vertex,oldMeshSizeID,(void**)&oldSize)) {
           cout<<"\nerror in setIsotropicSizeField(...) : no old mesh size data attached";
           exit(0);
        }
        double newSize = adaptFactorIso*(*oldSize);
//        if (newSize > MaxCoarsenFactor) newSize = MaxCoarsenFactor;
//        if (newSize < MaxRefineFactor) newSize = MaxRefineFactor;
        
        //Now that we know the isotropic size based on the error field
        //we will scale it with the size field from the hessian
        //currently using iso size for the max size value and then scale
        //accordingly
        
        hess[i].h[1] = newSize*hess[i].h[1]/hess[i].h[0];
        hess[i].h[2] = newSize*hess[i].h[2]/hess[i].h[0];
        hess[i].h[0] = newSize;

         for(int k=0; k<3; k++) {
           if (hess[i].h[k] > hmax) hess[i].h[k] = hmax;
           if (hess[i].h[k] < hmin) hess[i].h[k] = hmin;
         }   
       }  // if option = 4

      double localHmin = hmax;
      double* h2 = new double[3];
      double* dir2 = new double[9];

      for(j=0; j<3; j++) {
      // for option 9
        if(localHmin > hess[i].h[j])
       	localHmin = hess[i].h[j];

        if(ABS(hess[i].h[j]-hmax) <= tol2)
       	foundHmax = 1;
        if(ABS(hess[i].h[j]-hmin) <= tol3)
       	foundHmin = 1;


       h2[j]=hess[i].h[j];
       for(int k=0; k<3; k++) {
          dir2[j*3+k] = hess[i].dir[j][k];
       } 
     } //loop on j
 
      if(foundHmin)
        hminCount++;
      if(foundHmax)
        hmaxCount++;
      if(foundHmin && foundHmax)
        bothHminHmaxCount++;

      EN_attachDataPtr( (pEntity)vertex, nodalSizeID, (void *)h2);
      EN_attachDataPtr( (pEntity)vertex, nodalDirectionID, (void *)
                                 dir2);

      ++i;
    } // Vertex loop
    VIter_reset(vit);
    i=0;

#ifdef DEBUG    
//  M_writeVTKFile(mesh,"nodalSize",nodalSizeID,3);
#endif    
//  SizePropogate(mesh);

   int numSmooth=3;
   for (int k=0; k<50; k++){
     SmoothSize(mesh,numSmooth); //Size field smoothing similar to hessians    
     commuSmoothSize(pmesh, mesh,numSmooth);
   }
   
   if(isSizeLimit)
      SizeLimit(mesh);
/*  
   for (int k=0; k<5; k++) {
      SmoothDir(mesh);
      commuSmoothDir(pmesh, mesh); 
   }
*/

  M_writeVTKFile(mesh,"nodalSize",nodalSizeID,3);
  while ( vertex=VIter_next(vit)) {

    // Modify mesh-metric to take special care
    // like in parallel plates (i.e, w=f(y)) 
    // if a user wants different hmax in x-direction
    // other examples can be situations where
    // user don't want to have one mesh edge 
    // connecting two geometric model edges (like no dofs)
    // ModifyMetric(vertex,hess[i].dir,hess[i].h);

    // set the data in MeshSim format
    // M'= R*Lambda   (R:eigenvecs, Lambda:eigenVals (diag)

    double* h2 = new double[3];
    double* dir2 = new double[9];

    EN_getDataPtr((pEntity)vertex,nodalSizeID ,
                      (void**)&h2);
    EN_getDataPtr((pEntity)vertex,nodalDirectionID ,
                      (void**)&dir2); 

    for (int iRow=0; iRow<3; iRow++) {

       if(h2[iRow] < hmin) h2[iRow] = hmin;
       if(h2[iRow] > hmax) h2[iRow] = hmax;

       hess[i].h[iRow] = h2[iRow];
//      hess[i].h[iRow] = localHmin;
       for(int iDir=0; iDir<3; iDir++) {
         hess[i].dir[iRow][iDir]=dir2[iRow*3+iDir]*hess[i].h[iRow];
       }
    }

    EN_modifyDataPtr( (pEntity)vertex, nodalSizeID, (void *)h2);

#if  ( defined  DEBUG )
    // double c[3];
    // V_coord(vertex,c);
    // printf("%i %f %f %f ",EN_id((pEntity)vertex),c[0],c[1],c[2]);
    // printf("%f %f %f ",hess[i].dir[0][0],hess[i].dir[0][1],hess[i].dir[0][2]);
    // printf("%f %f %f ",hess[i].dir[1][0],hess[i].dir[1][1],hess[i].dir[1][2]);
    // printf("%f %f %f\n",hess[i].dir[2][0],hess[i].dir[2][1],hess[i].dir[2][2]);
#endif  

   //Here we calculate the current size of the original mesh to compare
   //with our size field request. 
/*   
   double* OrgSize = new double;
   double OrgAnisoSize[3][3];
   int iSize;

   iSize = V_size(vertex, OrgSize, OrgAnisoSize);
   double* Orgh = new double[3];

   for(int k=0; k<3; k++) {
        Orgh[k]=0.0;
   }

   if (iSize == 1){
    for(int k=0; k<3; k++) {
        Orgh[k]= *OrgSize;
    }
   }

   if (iSize == 2){
     for (int k=0;k<3;k++){
        for (int j=0 ;j<3; j++) {
           Orgh[k] +=(OrgAnisoSize[k][j])*(OrgAnisoSize[k][j]);
        }
        Orgh[k] = sqrt(Orgh[k]);
      }
   }
   EN_attachDataPtr( (pEntity)vertex, OrgSizeID, (void *)
                      Orgh);
*/    
    if(!preLBforAdaptivity) {
      if(option==1 || option==2 || option==4) {
      	MSA_setAnisoVertexSize(simAdapter,
			       vertex,
			       hess[i].dir);

// /* compute the size field        
    double* h2 = new double[3];
    for(int k=0; k<3; k++) {
        h2[k]=0.0;
    }
    for (int k=0;k<3;k++){
       for (int j=0 ;j<3; j++) {
           h2[k] +=(hess[i].dir[k][j])*(hess[i].dir[k][j]);
       }
       h2[k] = sqrt(h2[k]);
     }
     EN_attachDataPtr( (pEntity)vertex, nodalSizeID, (void *)
                       h2);
    }
/*      
      else if(option==9) {
	MSA_setVertexSize(simAdapter,
			  vertex,
			  localHmin);
      }
*/      
    }
    else {
      if(option==1  || option==2) {
	double *metric = new double[9];

	for(int iRow=0; iRow<3; iRow++)
	  for(int iDir=0; iDir<3; iDir++)
	    metric[iRow*3+iDir] = hess[i].dir[iRow][iDir];

	EN_attachDataPtr((pEntity)vertex,meshSizeID,metric);
      }
/*      
      else if(option==9) {
	double *size = new double;
	*size = localHmin;

	EN_attachDataPtr((pEntity)vertex,meshSizeID,size);
      }
*/
    }

    ++i;
  }
  VIter_delete(vit);
  
//  SizeQuality(mesh);
#ifdef DEBUG  
//  M_writeVTKFile(mesh,"nodalSize",nodalSizeID,3);
//  M_writeVTKFile(mesh, "nodalDirection", nodalDirectionID, 9);
#endif  
//  if(PMU_size()==1) {
//     pMesh meshMerge;
//     meshMerge = M_createFromParMesh(pmesh, 3, prog);
//     M_write(mesh, "mesh_size.sms", 0, prog);
//     M_release(meshMerge);
//  } else {
     PM_write(pmesh, "mesh_size.sms", sthreadNone, prog);
//  }

  printf(" [%d] Nodes with hmin into effect : %d (%4.2f%%)\n",PMU_rank(),hminCount,((double)hminCount/nshg)*100);
  printf(" [%d] Nodes with hmax into effect : %d (%4.2f%%)\n",PMU_rank(),hmaxCount,((double)hmaxCount/nshg)*100);
  printf(" [%d] Nodes with both hmin/hmax into effect : %d (%4.2f%%)\n",PMU_rank(),bothHminHmaxCount,((double)bothHminHmaxCount/nshg)*100);
  printf("\n");

//   adaptSimLog<<"Nodes with hmin into effect : "<<hminCount<<endl;
//   adaptSimLog<<"Nodes with hmax into effect : "<<hmaxCount<<endl;
//   adaptSimLog<<"Nodes with both hmin/hmax into effect : "<<bothHminHmaxCount<<endl;
//   adaptSimLog.close();

  fprintf(adaptSimLog,"Nodes with hmin into effect : %d (%4.2f%%)\n",hminCount,((double)hminCount/nshg)*100);
  fprintf(adaptSimLog,"Nodes with hmax into effect : %d (%4.2f%%)\n",hmaxCount,((double)hmaxCount/nshg)*100);
  fprintf(adaptSimLog,"Nodes with both hmin/hmax into effect : %d (%4.2f%%)\n",bothHminHmaxCount,((double)bothHminHmaxCount/nshg)*100);
  fclose(adaptSimLog);

  // assuming one partition per proc
  int gid = PMU_gid(PMU_rank(),0);

  // hess is a Hessian struct
  // 3rd arg is timestep
  writeMEDITSizeField(hess,mesh,timeStepNumber,gid);

  delete [] hess;

  if(option==1 || option==2 || option==9) {
    cleanAttachedData(mesh,nodalGradientID,0);
    cleanAttachedData(mesh,nodalHessianID,0);
    cleanAttachedData(mesh,localPatchVolID,0);
    cleanAttachedData(mesh,localGradientID,0);
    cleanAttachedData(mesh,localHessianID,0);
    cleanAttachedData(mesh,numSurroundNodesID,0);
    cleanAttachedData(mesh,locMaxInterpolErrorID,0);
    cleanAttachedData(mesh,globMaxInterpolErrorID,0);  
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(preLBforAdaptivity) {
#ifndef ibm
    printf("\n[%2d] memory usage before partioning: %d (KB)\n\n",PMU_rank(),phParAdaptProcSize());
#endif
    int partLB4AptOption;
    if(option==1 || option==2)
       partLB4AptOption = 1;
    else if(option==9)
      partLB4AptOption = 9;
    partitionMeshToLoadBalanceForAdaptivity(pmesh,mesh,partLB4AptOption,nErrorVars);

#ifndef ibm
    printf("\n[%2d] memory usage after partioning: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
    pMesh LBmesh = PM_mesh(pmesh,0);
    mesh = LBmesh;

    VIter vIter = M_vertexIter(LBmesh);
    while(vertex = VIter_next(vIter)) {

      double *metric;
      if(!EN_getDataPtr((pEntity)vertex,meshSizeID,(void**)&metric)) {
	printf("\nerror in setSizeFieldUsingHEssians: no data attached with meshSizeID to vertex\n");
	exit(0);
      }

      if(option==1 || option==2) {
	double simMetric[3][3];

	for(int iRow=0; iRow<3; iRow++)
	  for(int iDir=0; iDir<3; iDir++)
	    simMetric[iRow][iDir] = metric[iRow*3+iDir];

	MSA_setAnisoVertexSize(simAdapter,vertex,simMetric);
      }
      else if(option==9) {

	MSA_setVertexSize(simAdapter,vertex,*metric);
      }
    }
    VIter_delete(vIter);

    if(option==1 || option==2) {
      cleanAttachedData(LBmesh,meshSizeID,0);
    }
    else if(option==9) {
      cleanAttachedData(LBmesh,meshSizeID,0,0);
    }
  }
  MD_deleteMeshDataId(meshSizeID);

  MPI_Barrier(MPI_COMM_WORLD);
}

void SizeQuality(pMesh mesh){
   
  pVertex vertex;
  pVertex otherVert;
  pEdge edge;
   double* h = new double[3];  
   double* h2 = new double[3]; 
   VIter vit=M_vertexIter(mesh);
   while(vertex=VIter_next(vit)) {
      EN_getDataPtr((pEntity)vertex, nodalSizeID, (void**)&h);
      int numEdges = V_numEdges(vertex);

      double* averageValue = new double[3];
      for (int j=0; j<3; j++){
        averageValue[j] = 0.0;
      }

      for (int i=0; i < numEdges; i++) {
         edge = V_edge(vertex,i);
         otherVert = E_otherVertex(edge,vertex);
         EN_getDataPtr((pEntity)otherVert, nodalSizeID, (void**)&h2);
         
// Calculating the average size variation at a vertex
//         printf("%f %f \n",h[0],h2[0]);
         for(int j=0; j<3; j++){
           averageValue[j] = averageValue[j] + fabs(h[j]-h2[j]);
         }
         
      }

      for(int j=0; j<3; j++){
        averageValue[j] = averageValue[j]/numEdges;
      }
//      printf("%f \n",averageValue[0]);      
      EN_attachDataPtr((pEntity)vertex, QualitySizeID, (void *)averageValue);
  }

#ifdef DEBUG   
//  M_writeVTKFile(mesh,"QualitySize", QualitySizeID, 3);
#endif   
}

#ifdef __cplusplus
}
#endif

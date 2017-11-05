#include "phParAdapt.h"
#include "func.h"
#include "Eigen.h"
#include <iostream>
#include <fstream>
#include "attachData.h"
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#ifdef SIM
#include "SimMeshTools.h"
#endif

using namespace std;

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern double* wght;
extern int nErrorVars;

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;
extern pMeshDataId nodalSizeID;

extern pMeshDataId numSurroundNodesID;
extern pMeshDataId oldMeshSizeID;

extern pMeshDataId meshSizeID;

extern int preLBforAdaptivity;
///Joe
extern double epsilon;
extern double max_size;
extern double min_size;
extern int size_flag;
extern int isSizeLimit;
extern int MaxLimitFact;
extern int MinLimitFact;

void setIsotropicSizeField(pParMesh pmesh,
			   pMesh mesh,
			   pMSAdapt simAdapter,
			   double factor,
			   double hmax, 
			   double hmin,
			   int option)
{
  // for 3D problems
  int dim = 3;
  // assuming linear basis
  int poly_order = 1;

  double totalError = 0.;
  double threshold  = 0.;
  double sumOfError = 0.;

  // get error info. from all partitions
  getGlobalErrorInfo(mesh,totalError,sumOfError);

  oldMeshSizeID = MD_newMeshDataId("isotropic old mesh size");

  computeOldMeshSize(mesh,option);
  MPI_Barrier(MPI_COMM_WORLD);

  commuOldMeshSize(pmesh,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  char size_file[256], adaptfactor_file[256];

  // make sure there is only one partition per proc.
  if(PM_numParts(pmesh) > 1) {
      if(PMU_rank() == 0) {
	cout<<"\nError in setIsotropicSizeField()... "<<endl;
	cout<<"only one part per proc at the moment allowed \n"<<endl;
      }
      SimPartitionedMesh_stop();
      exit(1);
  }

#ifdef DEBUG
  // assuming 1 part. on each proc.
  sprintf(size_file,"isoSize.%d.dat",PMU_rank()+1);
  sprintf(adaptfactor_file,"adaptFactor.%d.dat",PMU_rank()+1);

  ofstream sizes(size_file);
  ofstream adaptFactorFile(adaptfactor_file);  
#endif

  if(preLBforAdaptivity)
    meshSizeID = MD_newMeshDataId("mesh size ID");

  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
    double *nodalValue;
    if(!EN_getDataPtr((pEntity)vertex,errorIndicatorID,(void**)&nodalValue)) {
      cout<<"\nerror in setIsotropicSizeField(...) : no error data attached to vertex"<<endl;
      exit(0);
    }

    threshold = factor*totalError;
    // based on mesh optimal criteria and error convergence rate
    // implemented for linear cases
    // refer L.-Y. Li, Comm. Num. Meth. Eng.,
    // Vol. 11, 857-868 & 911-915 (1995)
    double adaptFactor = threshold/(pow(fabs(*nodalValue),0.25)*sqrt(sumOfError));

    double MaxRefineFactor = hmin;
    double MaxCoarsenFactor = hmax;

    // these factors set cut-off levels for adaptation
    // can set such that no coarsening (or refinement) occurs
/*    if (MaxRefineFactor > 1.e-10  && adaptFactor < 1.0/MaxRefineFactor) 
      adaptFactor = 1.0/MaxRefineFactor;

    if (MaxCoarsenFactor > 1.e-10 && adaptFactor > MaxCoarsenFactor) 
      adaptFactor = MaxCoarsenFactor;
*/
    double *oldSize;
    if(!EN_getDataPtr((pEntity)vertex,oldMeshSizeID,(void**)&oldSize)) {
      cout<<"\nerror in setIsotropicSizeField(...) : no old mesh size data attached to vertex"<<endl;
      exit(0);
    }
    
    double newSize;
    if (option == 10){
      if (*nodalValue > factor) newSize = *oldSize/3;
      else newSize = *oldSize;
    }
    else if (option == 11 ){
      //double coord[3];
      //double plane;
      //V_coord(vertex,coord); 
      //plane = -0.516616558513076*coord[0]+0.787121033907457*coord[1]+0.336968558548957*coord[2]+0.761858682657;
      //if (*nodalValue > factor && plane >= 0.0) newSize = *oldSize/2;
      //else  newSize = *oldSize;
      //
      if (*nodalValue > factor) newSize = *oldSize/2;
      else  newSize = *oldSize;
    }
    else {   
      newSize = adaptFactor*(*oldSize);
    }


    if (newSize > MaxCoarsenFactor) {
      // If the newSize is smaller than max theshold, set the newSize back to the theshold
      newSize = MaxCoarsenFactor;
      // But preserve the old size that was already larger than the threshold so that it does not get refined
      if (*oldSize > MaxCoarsenFactor)  newSize = *oldSize;
    }

    if (newSize < MaxRefineFactor) {
      // If the newSize is smaller than min theshold, set the newSize back to the theshold
      newSize = MaxRefineFactor;
      // But preserve the old size that was already smaller than the theshold so that it does not get coarsened!
      if ( *oldSize < MaxRefineFactor) newSize = *oldSize;        
    }

    if(isSizeLimit) {
        if(newSize > *oldSize*MaxLimitFact) newSize = *oldSize*MaxLimitFact;
        if(newSize < *oldSize/MinLimitFact) newSize = *oldSize/MinLimitFact;
    }

//from Joe: setting sizefield based upon level set at nodes
    if (size_flag != 0) {
//Now we have to do the same with solution set
      double *nodalSolutionSet;
      double phi;
      if(!EN_getDataPtr((pEntity)vertex, phasta_solution,(void**)&nodalSolutionSet)){
         cout<<"\nerror in setIsotropicSizeField(...) : no data attached to vertex\n";
         V_info(vertex);
         exit(0);
      }
      phi = sqrt(nodalSolutionSet[1]*nodalSolutionSet[1]);
      if (phi < epsilon) {
        newSize = min_size;
      }
      else if (phi < 3*epsilon) {
        newSize = (min_size + max_size) / 2.0;
      }
      else {
        newSize = max_size;
      } 
//      V_coord(vertex,xyz);
//      if (sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])> epsilon) {
    }
    double dir[3][3];
    if(option == 3) {
       dir[0][0] = newSize; dir[1][1] = newSize; dir[2][2] = newSize;
       dir[0][1] = dir[0][2] = 0; 
       dir[1][0] = dir[1][2] = 0;
       dir[2][0] = dir[2][1] = 0;
       MSA_setAnisoVertexSize(simAdapter, vertex, dir);
    }   
    double* sizeField = new double;
    if(option != 3) {
     if(!preLBforAdaptivity) {
//KEDAR: MSA_setVertexSize delete from here to call SmoothSize
      *sizeField = newSize;
//      printf("size: %lf\n", newSize);
     }
     else {
      double *size = new double;
      *size = newSize;
      EN_attachDataPtr((pEntity)vertex,meshSizeID,size);
     }
    }
    EN_attachDataPtr((pEntity)vertex,nodalSizeID,(void *)sizeField);

#ifdef DEBUG
//    sizes<<newSize<<"\n";
//    adaptFactorFile<<adaptFactor<<"\n";
#endif

  }
  VIter_reset(vIter);

  M_writeVTKFile(mesh,"nodalSizeB",nodalSizeID,3);

//KEDAR: SmoothSize moved outside the vertex iterator.
//SetVertexSize needs to be called after again 
  int numSmooth=1;
  for (int k=0; k<10; k++){
     SmoothSize(mesh,numSmooth); //Size field smoothing similar to hessians    
     commuSmoothSize(pmesh, mesh,numSmooth);
      if(PMU_rank()==0) {
         cout<<"Size Field Smoothing iteration      : "<<k<<endl;
      }
  }
  M_writeVTKFile(mesh,"nodalSizeA",nodalSizeID,3);

//now set sizes with simmetrix

  double* OrgSize = new double;
  double OrgAnisoSize[3][3];
  int iSize;
  double *oldSize;
  double* h = new double;
  double ratThresh=0.9; // not certain of the best number here as smoothing was applied to the original size

  while ( vertex=VIter_next(vIter)) {
    EN_getDataPtr((pEntity)vertex,oldMeshSizeID,(void**)&oldSize);
    EN_getDataPtr((pEntity)vertex,nodalSizeID ,
                      (void**)&h);
    double sizeRat;
    sizeRat= h[0]/(*oldSize);
    if(sizeRat <= ratThresh){
      iSize = V_size(vertex, OrgSize, OrgAnisoSize);
//      if (iSize == 1){   // this vertex is isotropic
//         MSA_setVertexSize(simAdapter, 
//                        vertex,
//                        h[0]);
//         }
//      if (iSize == 2){
          for (int k=0 ;k<3; k++) {
            for (int j=0 ;j<3; j++) {
               OrgAnisoSize[k][j]=0.0;
            }
          }
//          OrgAnisoSize[0][2]= 0.025;
//          OrgAnisoSize[1][1]= 0.1; // <-now try....worked -> 0.025;
//          OrgAnisoSize[2][0]= 0.1;
          OrgAnisoSize[0][1]= 0.025;
          OrgAnisoSize[1][2]= 0.025; 
          OrgAnisoSize[2][0]= 0.1;
          MSA_setAnisoVertexSize(simAdapter, 
                        vertex,
                        OrgAnisoSize);

//      }
    }
  }
  VIter_delete(vIter);
  delete [] h;

#ifdef DEBUG  
//  M_writeVTKFile(mesh, "IsotropicSize", nodalSizeID, 1);
#endif
#ifdef DEBUG
  sizes.close();
  adaptFactorFile.close();
#endif
   pProgress prog;
   /*if(PMU_size()==1) {
      pMesh meshMerge;
      cout << "\n converting pParMesh to pMesh here " << endl;
      meshMerge = M_createFromParMesh(pmesh, 3, prog);
      M_write(meshMerge, "mesh_size.sms", 0, prog);
      M_release(meshMerge);
   } else {*/
     PM_write(pmesh, "mesh_size.sms", prog);
   //}

  // data (single double ptr.) is attached to vertices
  cleanAttachedData(mesh,numSurroundNodesID,0,0);
  cleanAttachedData(mesh,oldMeshSizeID,0,0);
  MD_deleteMeshDataId(oldMeshSizeID);

  if(PMU_rank()==0) {
    cout<<"\nInfo. on adaptation parameters are: "<<endl;
    cout<<"Total Error      : "<<totalError<<endl;
    cout<<"factor           : "<<factor<<endl;
    cout<<"Weights          : ";
    for(int iEVar=0; iEVar<nErrorVars; iEVar++)
      cout<<wght[iEVar]<<" ";
    cout<<endl;
    cout<<"Threshold        : "<<threshold<<endl;
    cout<<"sumOfError       : "<<sumOfError<<endl;
    cout<<"MaxCoarsenFactor : "<<hmax<<endl;
    cout<<"MaxRefineFactor  : "<<hmin<<"\n"<<endl;
  }

  char log_file[256];

  // make sure there is only one partition per proc.
  if(PM_numParts(pmesh) > 1) {
      if(PMU_rank() == 0) {
	cout<<"\nError in setIsotropicSizeField()... "<<endl;
	cout<<"only one part per proc at the moment allowed \n"<<endl;
      }
      SimPartitionedMesh_stop();
      exit(1);
  }

#ifdef DEBUG
  sprintf(log_file,"phAdapt.%d.log",PMU_rank()+1);

  // assuming 1 part. on each proc.
  ofstream adaptSimLog(log_file);

  // as of now log file is identical for all procs
  adaptSimLog<<"Strategy chosen is size-field driven for isotropic adaptation"<<endl;  
  adaptSimLog<<"Info. on adaptation parameters are: "<<endl;
  adaptSimLog<<"Total Error      : "<<totalError<<endl;
  adaptSimLog<<"factor           : "<<factor<<endl;  
  adaptSimLog<<"Weights          : ";
  for(int iEVar=0; iEVar<nErrorVars; iEVar++)
    adaptSimLog<<wght[iEVar]<<" ";
  adaptSimLog<<endl;
  adaptSimLog<<"Threshold        : "<<threshold<<endl;
  adaptSimLog<<"sumOfError       : "<<sumOfError<<endl;
  adaptSimLog<<"MaxCoarsenFactor : "<<hmax<<endl;
  adaptSimLog<<"MaxRefineFactor  : "<<hmin<<endl;
  adaptSimLog.close();
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  if(preLBforAdaptivity) {
#ifndef ibm
    printf("\n[%2d] memory usage before partioning: %d (KB)\n\n",PMU_rank(),phParAdaptProcSize());
#endif
    // 10 EIs transformed to one scalar
    int nErrorVarsMod = 1;

    // 3rd arg. tells tensorial order of the size field 
    // (i.e., scalar or 3x3 mesh metric)
    // value "9" is for scalar size field
    partitionMeshToLoadBalanceForAdaptivity(pmesh,mesh,9,nErrorVarsMod);

#ifndef ibm
    printf("\n[%2d] memory usage after partioning: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
    pMesh LBmesh = PM_mesh(pmesh,0);
    mesh = LBmesh;

    VIter vIter = M_vertexIter(LBmesh);
    while(vertex = VIter_next(vIter)) {

      double *size;
      if(!EN_getDataPtr((pEntity)vertex,meshSizeID,(void**)&size)) {
	printf("\nerror in setSizeFieldUsingHEssians: no data attached with meshSizeID to vertex\n");
	exit(0);
      }

      MSA_setVertexSize(simAdapter,vertex,*size);
    }
    VIter_delete(vIter);

    cleanAttachedData(LBmesh,meshSizeID,0,0);
    MD_deleteMeshDataId(meshSizeID);
  }

}

void getGlobalErrorInfo(pMesh mesh, double& totalError, double& sumOfError)
{

  double totalErrorLoc = 0., sumOfErrorLoc = 0.;  

  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  // get error info. on each partition
  // vertices on partition bdry. will contribute only to owner proc.
  while( vertex = VIter_next(vIter) ) {
    if( EN_isOwnerProc((pEntity)vertex) ) {
      double *nodalValue;
      if(!EN_getDataPtr((pEntity)vertex,errorIndicatorID,(void**)&nodalValue)) {
	cout<<"\nerror in setIsotropicSizeField(...) : no error data attached to vertex"<<endl;
	exit(0);
      }
      
      // eta_domain^2 = sum_k (eta_k^2)
      totalErrorLoc += *nodalValue;
//      cout << "Sum of Err: " << totalErrorLoc << endl;
      // eta_sum = sum_k (eta_k)
      sumOfErrorLoc += sqrt(fabs(*nodalValue));
//      cout << "Sum of Err: " << sumOfErrorLoc << endl;
    }
  }
  VIter_delete(vIter);

  MPI_Barrier(MPI_COMM_WORLD);
  // communicate the error info. over patitions
  MPI_Allreduce(&totalErrorLoc, &totalError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sumOfErrorLoc, &sumOfError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // eta_domain = sqrt( sum_k (eta_k^2) )
  totalError = sqrt(totalError);
}

void computeOldMeshSize(pMesh mesh, int option)
{

  // option is not being used currently
  // can be utilized later to compute old size field 
  // based on different criteria, like average edge length or
  // average element inner radius etc.

  pVertex vertex;
  pEdge edge;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
    double *oldSize = new double;
    int *numSurrNodes = new int;
    *oldSize = 0.;
    *numSurrNodes = 0;

    if(!EN_isOnPartBdry((pEntity)vertex)) {
      // if vertex is not on partition bdry. 
      // treat is as in serial case
            
      int numEdges = V_numEdges(vertex);
      *numSurrNodes = numEdges;
      // old size at a location (vertex)
      // is mean value of the lengths 
      // of edges around (seems ok for isotropic meshes)
      // (can have different choices)
      for (int i=0; i < numEdges; i++) {
	edge = V_edge(vertex,i);
	*oldSize += E_length(edge);
      }
    }
    else {
      int numEdges = V_numEdges(vertex);
      // old size at a location (vertex)
      // is mean value of the lengths 
      // of edges around (seems ok for isotropic meshes)
      // (can have different choices)
      for (int i=0; i < numEdges; i++) {
	edge = V_edge(vertex,i);

	// edge length must be taken into consideration ONLY if it is owned by the proc.
	if(EN_isOwnerProc((pEntity)edge)) {
	  *numSurrNodes += 1;
	  *oldSize += E_length(edge);
	}
      }
    }

    EN_attachDataPtr((pVertex)vertex,oldMeshSizeID,(void *)oldSize);
    EN_deleteData( (pEntity)vertex, numSurroundNodesID);
    EN_attachDataPtr((pVertex)vertex,numSurroundNodesID,(void *)numSurrNodes);
  }
  VIter_delete(vIter);
}


#ifdef __cplusplus
}
#endif

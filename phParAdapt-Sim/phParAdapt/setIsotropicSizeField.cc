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
#include <strings.h>
#ifdef SIM
#include "SimMeshTools.h"
#include "SimAdvMeshing.h"
#endif

using namespace std;

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern double* wght;
extern int nErrorVars;
extern int AnisoSimmetrix;
extern double sizeRatio;
extern double ratioThresh;
extern int numSmooth;

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

void setIsotropicSizeField(pGModel model,
                           pParMesh pmesh,
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
// before this point oldSize was computed based on averaging sizes of edges around a
// node but Simmetrix has a function for that now.  Let's test to see if it is better
//   int itype=2; // we want anisotropic
// double anisosize[3][3];
//   itest=V_estimateSize(vertex, itype, NULL, anisosize);
   int itest,itype=1; // we want isotropic
   double size;  
   itest=V_estimateSize(vertex, itype, &size, NULL);
   *oldSize=size; // replace oldSize with Simmetrix computation
    
    double newSize;
    if (option == 10){
      if (*nodalValue > factor) newSize = *oldSize/3;
      else newSize = *oldSize;
    }
    else if (option == 11 ){
      double coord[3];
      //double plane;
      V_coord(vertex,coord); 
/*      double ptCheck[3]; // This is a point that should not be ref 38 11.5 4.93675
      ptCheck[0]=0.99; //38.025 11.625 4.975
      ptCheck[1]=0.40;
      ptCheck[2]=0.05;
      int stop;
      double distTop = sqrt(dist(coord, ptCheck));
      if(distTop < 5e-2) {
         stop=1;
      }  
*/
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


// Note, the logic below will block refinement of cells that are larger than input hmax
// which is far from the users intent when they choose a max element size
//  for now I am blocking it for 
   if(0) {
    if (newSize > MaxCoarsenFactor) {
      // If the newSize is smaller than max theshold, set the newSize back to the theshold
      newSize = MaxCoarsenFactor;
      // But preserve the old size that was already larger than the threshold so that it does not get refined
      if (*oldSize > MaxCoarsenFactor)  newSize = *oldSize;
    }
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
//  int numSmooth=1;
//  int iSmooth=1;
  if (numSmooth>0) {
     for (int k=0; k<numSmooth; k++){
        SmoothSize(mesh,numSmooth); //Size field smoothing similar to hessians    
        commuSmoothSize(pmesh, mesh,numSmooth);
        if(PMU_rank()==0) {
           cout<<"Size Field Smoothing iteration      : "<<k<<endl;
        }
     }
  }
  M_writeVTKFile(mesh,"nodalSizeA",nodalSizeID,3);

  double ratThresh=ratioThresh; // 0.75; // not certain of the best number here as smoothing was applied to the original size
  int Isotrop;
  double sizeRat=1.0;
  double coordvcur[3];
  int iSize,icountVertsNotBL,icountVertsOnFaces,icountIsotrop,icountAnisotrop;
  icountVertsNotBL=0;
  icountVertsOnFaces=0;
  icountIsotrop=0;
  icountAnisotrop=0;
  double* h = new double;
  double *oldSize;

// START OF WRITING TAGS
  char tfile[255];
  bzero( (void*)tfile, 255 );
  sprintf(tfile,"edgeids.%d", PMU_rank());
  FILE * itf;
  itf=fopen (tfile,"wt");
//  openfile_(tfile,"write",&itf);
  int taggedInt;
//CWS  Std:vector <int> vals;
//  int icountTags=0;
//END OF WRITING TAGS
/* from CWS   

int icountTags=0;
Std:vector<int> vals;
inside of the Loop that finds the write ints {
  vals.push_back(anInt)
  icountTags++
}

Then when it is time to write after the loop is complete:
           fprintf(itf,"%d \n",icountTags);
followed by this command which writes all the ints collected.
fwrite(&Vals[0], sizeof(int), Vals.size(), filehandle)

Use fread to read them:
http://www.cplusplus.com/reference/cstdio/fread/

The vector knows it's size and will resize itself as needed.
*/
//now set sizes with simmetrix outside of the BL
  while ( vertex=VIter_next(vIter)) {
      double coord[3];
      //double plane;
      V_coord(vertex,coord); 
/*      double ptCheck[3]; // This is a point that should not be ref 38 11.5 4.93675
      ptCheck[0]=0.99; //38.025 11.625 4.975
      ptCheck[1]=0.40;
      ptCheck[2]=0.05;
      int stop;
      double distTop = sqrt(dist(coord, ptCheck));
      if(distTop < 5e-2) {
         stop=1;
      }  
*/
   if(!EN_isBLEntity((pEntity)vertex)) { // true if this is a NOT BL entity
     EN_getDataPtr((pEntity)vertex,oldMeshSizeID,(void**)&oldSize);
     EN_getDataPtr((pEntity)vertex,nodalSizeID ,
                      (void**)&h);
     sizeRat= h[0]/(*oldSize);
     icountVertsNotBL++;
     if(sizeRat<ratThresh){   // only set if size ratio under threshold
         sizeRat=sizeRatio; // flat plate suggests that Simm overestimates sizes this and the 1.32 below correct to expected size but not sure if universal
       icountIsotrop++;
       if(AnisoSimmetrix==4){  // setRefineLevel way
         int numEdges = V_numEdges(vertex);
         pEdge edge;
         int level=1;
         for (int i=0; i < numEdges; i++) {
	   edge = V_edge(vertex,i);
           MSA_setRefineLevel(simAdapter, (pEntity)edge, level);	
// START INT PER LINE
//CWS           vals.push_back(EN_id((pEntity)edge));
//           icountTags++;
           taggedInt=EN_id( (pEntity)edge);
           fprintf(itf,"%d \n",taggedInt);
// END INT PER LINE

           if(0){ 
// recurse one level out
           pEdge edgeO;
           pVertex vother;
           vother=E_otherVertex(edge, vertex);	
           int numEdgesO = V_numEdges(vother);
           for (int j=0; j < numEdgesO; j++) {
              edgeO = V_edge(vother,j);
              MSA_setRefineLevel(simAdapter, (pEntity)edgeO, level);	
           }
           }
         }
       } else{
       if(1) { // gaps in BL's are being missed in the else
         int itype=2; // we want anisotropic
         double size;  // failed so I passed null
         double anisosize[3][3];
         int itest=V_estimateSize(vertex, itype, NULL, anisosize);
         if(itest) {
         anisosize[0][0]*= sizeRat;
         anisosize[0][1]*= sizeRat;
         anisosize[0][2]*= sizeRat;
         anisosize[1][0]*= sizeRat;
         anisosize[1][1]*= sizeRat;
         anisosize[1][2]*= sizeRat;
         anisosize[2][0]*= sizeRat;
         anisosize[2][1]*= sizeRat;
         anisosize[2][2]*= sizeRat;
         MSA_setAnisoVertexSize(simAdapter, 
                            vertex,
                           anisosize);
        } else {
           cout<<" V_estimateSize returned 0: "<<endl;
        }  
       }  else {
         h[0]=*oldSize*sizeRatio;  // lets user control how much to reduce, not just 1/2
         MSA_setVertexSize(simAdapter, 
			vertex,
			h[0]);
      } // else from if(1)
      } // else from setSize way
     } // sizeRat < ratThresh
    } // not a bl entity
  } // vertex iterator
  VIter_delete(vIter);
//  fprintf(itf,"%d \n",icountTags);
//CWS  fwrite(&vals[0], sizeof(int), vals.size(),itf);
//  fclose (itf);

//now set sizes with simmetrix inside of the BL
  int maxE,minE,midE;
  double* OrgSize = new double;
  double OrgAnisoSize[3][3];


  if(PMU_rank()==0) {
    cout << "Starting aniso Adapt  with AnisoSimmetrix="<<AnisoSimmetrix << endl;
    cout << "ratioThresh="<<ratioThresh << endl;
    cout << "sizeRatio="<<sizeRatio << endl;
    cout << "numSmooth="<<numSmooth << endl;
  }
  
  pGFace gface;
  GFIter gfIter=GM_faceIter(model);
  pGRegion gregion;
  GRIter grIter=GM_regionIter(model);
  while ( gregion=GRIter_next(grIter)) {   // only way I know to find region
  while ( gface=GFIter_next(gfIter)) {
    VIter vofIter  = M_classifiedVertexIter(mesh, (pGEntity)gface, 1); // 1 gives closure 0 not....I think we want closure to pick up verts on model edges and verts
    while ( vertex=VIter_next(vofIter)) {
//  above 5 lines swap to a vertices on face iteration while below is all    
// reuse of reset iterator while ( vertex=VIter_next(vIter)) {

    V_coord(vertex, coordvcur ); 

// adding the ability to scan the growth curve to find the highest error.
   int itest,itesto,fromSide=0; // 1-> face normal points outside region, 0-> opposite
   pGEntity into;
/// copied from above for format  pVertex vertex;
   pEntity seed;
   sizeRat=1;
   icountVertsOnFaces++;
//
/*      double ptCheck[3]; // This is a point that should not be ref 38 11.5 4.93675
      ptCheck[0]=0.99; //38.025 11.625 4.975
      ptCheck[1]=0.40;
      ptCheck[2]=0.05;
      int stop;
      double distTop = sqrt(dist(coordvcur, ptCheck));
      if(distTop < 5e-2) {
         stop=1;
      }  
*/
//
   if(EN_isBLEntity((pEntity)vertex)) { // true if this is a BL entitGy
     if(BL_isBaseEntity((pEntity)vertex,(pGEntity)gface)) { // current vertex is base
       itest=BL_stackSeedEntity((pEntity)vertex,(pGEntity)gface,fromSide,(pGEntity)gregion,&seed);
       if(itest<1) { // try other side?
         fromSide=1;
         itesto=BL_stackSeedEntity((pEntity)vertex,(pGEntity)gface,fromSide,(pGEntity)gregion,&seed);
         if(itesto<1) { 
         sizeRat=101; // code for both sides fail
//  This will fail for all faces that BL's grow "up", e.g., symmetry planes 
// disable warning until we figure out how to block these faces.
//         cout << "WARNING: BL_stackSeedEntity failed" <<     endl;
         }
       }
       if(sizeRat<100) { // one side succeeded so traverse the stack
         pPList verStack=PList_new();
         pPList edgeStack=PList_new();
         itest=BL_growthVerticesAndEdges((pEdge)seed,verStack,edgeStack);
         void *iter = 0; // must initialize to 0
         double sizeRatMin=1.0;
         pVertex ent;
         pVertex entlast;
         while(ent = (pVertex)PList_next(verStack, &iter)){
             // process each item in list
           EN_getDataPtr((pEntity)ent,oldMeshSizeID,(void**)&oldSize);
           EN_getDataPtr((pEntity)ent,nodalSizeID ,
                      (void**)&h);
           sizeRat= h[0]/(*oldSize);
           sizeRatMin=min(sizeRat,sizeRatMin);
           entlast=ent;
          }
          sizeRat=sizeRatMin;
          PList_delete(verStack);
          PList_delete(edgeStack);
          if(AnisoSimmetrix > 0) {  // turn this on or off with 0 or 1
// inserted codeblock to use Simmetrix size
          if(sizeRat<ratThresh) {
            int itype=2; // we want anisotropic
            double size;  // failed so I passed null
            double anisosize[3][3];
            itest=V_estimateSize(vertex, itype, NULL, anisosize);
            if(itest==1) {
              if(AnisoSimmetrix==4){  // setRefineLevel way
                int numEdges = V_numEdges(vertex);
                pEdge edge;
                int level=1;
                for (int i=0; i < numEdges; i++) {
	          edge = V_edge(vertex,i);
                  MSA_setRefineLevel(simAdapter, (pEntity)edge, level);
// Test tag driven adapation following Saurabh's advice
                  taggedInt=EN_id( (pEntity)edge); 
                  fprintf(itf,"%d \n",taggedInt);	
// recurse one level out
                  if(0) {
                  pEdge edgeO;
                  pVertex vother;
                  vother=E_otherVertex(edge, vertex);	
                  int numEdgesO = V_numEdges(vother);
                  for (int j=0; j < numEdgesO; j++) {
                     edgeO = V_edge(vother,j);
                     MSA_setRefineLevel(simAdapter, (pEntity)edgeO, level);	
                  }
                 }
                }
                numEdges = V_numEdges(entlast);
                for (int i=0; i < numEdges; i++) {
	          edge = V_edge(entlast,i);
                  MSA_setRefineLevel(simAdapter, (pEntity)edge, level);
// Test tag driven adapation following Saurabh's advice
                  taggedInt=EN_id( (pEntity)edge);
                  fprintf(itf,"%d \n",taggedInt);	
// recurse one level out
                  if(0) {
                  pEdge edgeO;
                  pVertex vother;
                  vother=E_otherVertex(edge, entlast);	
                  int numEdgesO = V_numEdges(vother);
                  for (int j=0; j < numEdgesO; j++) {
                     edgeO = V_edge(vother,j);
                     MSA_setRefineLevel(simAdapter, (pEntity)edgeO, level);	
                  }
                 }
                }
              } else { // set size way
              sizeRat=sizeRatio; // flat plate suggests that Simm overestimates sizes this and the 1.32 below correct to expected size but not sure if universal
// saurabh says we should be careful to limit the smaller sizes from becoming larger than the size we are shrinking.  To do that we calculate the length of all three size vectors, and then, if their size will be larger than the future size of l2 then create a scale factor that makes them the same size viz:
              double l0,l1,l2;
              l1=sqrt(anisosize[1][0]*anisosize[1][0]
                     +anisosize[1][1]*anisosize[1][1]
                     +anisosize[1][2]*anisosize[1][2]);
              l2=sqrt(anisosize[2][0]*anisosize[2][0]
                     +anisosize[2][1]*anisosize[2][1]
                     +anisosize[2][2]*anisosize[2][2]);
              l0=sqrt(anisosize[0][0]*anisosize[0][0]
                     +anisosize[0][1]*anisosize[0][1]
                     +anisosize[0][2]*anisosize[0][2]);
              double sizeRat1=1.0; // /1.32;
              if(AnisoSimmetrix>1) sizeRat1=sizeRatio;
              double sizeRat0=sizeRat1; // /1.32;
              if(sizeRat1*l1>sizeRat*l2) sizeRat1*=l2*sizeRat/l1;
              if(sizeRat0*l0>sizeRat*l2) sizeRat0*=l2*sizeRat/l0;
              if(AnisoSimmetrix>2) sizeRat0=sizeRatio;
              anisosize[2][0]*= sizeRat;
              anisosize[2][1]*= sizeRat;
              anisosize[2][2]*= sizeRat;
              anisosize[1][0]*= sizeRat1;
              anisosize[1][1]*= sizeRat1;
              anisosize[1][2]*= sizeRat1;
              anisosize[0][0]*= sizeRat0;
              anisosize[0][1]*= sizeRat0;
              anisosize[0][2]*= sizeRat0;
              MSA_setAnisoVertexSize(simAdapter, 
                            vertex,
                           anisosize);
// also set on the top of the blstack
              itest=V_estimateSize(entlast, itype, NULL, anisosize);
              anisosize[2][0]*= sizeRat;
              anisosize[2][1]*= sizeRat;
              anisosize[2][2]*= sizeRat;
              anisosize[1][0]*= sizeRat1;
              anisosize[1][1]*= sizeRat1;
              anisosize[1][2]*= sizeRat1;
              anisosize[0][0]*= sizeRat0;
              anisosize[0][1]*= sizeRat0;
              anisosize[0][2]*= sizeRat0;
              MSA_setAnisoVertexSize(simAdapter, 
                            entlast,
                           anisosize);
              icountAnisotrop++;
             }
            } else {
              cout << "WARNING: V_estimateSize returned 0 " <<     endl;
            }
          }
          sizeRat=1.0; // this is a cheat to allow us to keep other option but skip since sizeRat > ratThresh on line below
       }   // end of inserted codeblock


       } // end of stack found
      } // not a base sizeRat left at 1.0 so no ref
    } // end of isBL


 // end stack scan for worst

// 4 lines below are what we did when we checked every vertex
//    EN_getDataPtr((pEntity)vertex,oldMeshSizeID,(void**)&oldSize);
//    EN_getDataPtr((pEntity)vertex,nodalSizeID ,
//                      (void**)&h);
// sizeRat= h[0]/(*oldSize);

    if(sizeRat <= ratThresh){

// begin of computation of current anisotropic size

     int numEdges = V_numEdges(vertex);
     double NormdotProdTable[numEdges][numEdges];
     double edgesIonV[numEdges][3];
//step 1: compute edge lengths to determine if anisotropic size needed
     double edgeL,edgemax,edgemin;
     double eLength [3];
     int ifirst, isecond, ithird;
     edgemax=0.0;
     edgemin=1.0e6;
     pEdge edge;
     for (int i=0; i < numEdges; i++) {
	edge = V_edge(vertex,i);
	edgeL= E_length(edge);
 	edgemax=max(edgemax,edgeL);
	edgemin=min(edgemin,edgeL);
     }
     if (0) { // always do else now edgemax < 4.0*edgemin) { // not really worth complicate aniso
       Isotrop=1; 
     } else {
      Isotrop=0;  // set this to 0 but note we switch to isotropic if needed    
//step 2: get and store all the edge vectors for the current vertex 
//     double coordvcur[3];
     double coordvother[3];
      pVertex vother;
      for (int i=0; i < numEdges; i++) {
        edge = V_edge(vertex,i);
        vother=E_otherVertex(edge, vertex);	
        V_coord(vother, coordvother );	
        for (int j=0 ;j<3; j++) {
               edgesIonV[i][j]=coordvother[j]-coordvcur[j];
        }
      }
//step 3 compute dot product table to find three most orthogonal edges
      double dotProdTable[numEdges][numEdges];
      for (int i=0; i < numEdges; i++) {
        for (int j=0; j < numEdges; j++) {
          dotProdTable[i][j]=dotProd(edgesIonV[i],edgesIonV[j]);
        }
      }
      double InvEdgeL[numEdges];
      int longestEdge;
      double InvEdgeMin=1e6;
      for (int i=0; i < numEdges; i++) {
         InvEdgeL[i]=1.0/sqrt(dotProdTable[i][i]);
         if (InvEdgeL[i]<InvEdgeMin) {
           InvEdgeMin=InvEdgeL[i];
           longestEdge=i;
         }
      }
//      double NormdotProdTable[numEdges][numEdges];
      for (int i=0; i < numEdges; i++) {
        for (int j=0; j < numEdges; j++) {
          NormdotProdTable[i][j]=dotProdTable[i][j]*
          InvEdgeL[i]*InvEdgeL[j];
        }
      }
      int edgeAlignCounts[numEdges][3];
      double high,med, low, temp;
      high=0.99;
      med=0.5;
      low=0.01;
      double tmp1;
      double projmin=1;
      if(1){
//  Simmetrix suggestion that longest edge should be one of vectors
      ifirst=longestEdge;
// for now, choose the second edge to be the one that is most orthogonal to first
      for (int i=0; i < numEdges; i++) {
           tmp1=fabs(NormdotProdTable[ifirst][i]);
           if(tmp1 < projmin) {
              projmin=tmp1;
              isecond=i;
           }
      }
// third vector is not an edge, rather it is orthogonal to the first 2 and then scaled by the longest "other" edge projection onto this third vector
      double size3[3];
      size3[0]= edgesIonV[ifirst][1]*edgesIonV[isecond][2]
              - edgesIonV[ifirst][2]*edgesIonV[isecond][1];
      size3[1]= edgesIonV[ifirst][2]*edgesIonV[isecond][0]
              - edgesIonV[ifirst][0]*edgesIonV[isecond][2];
      size3[2]= edgesIonV[ifirst][0]*edgesIonV[isecond][1]
              - edgesIonV[ifirst][1]*edgesIonV[isecond][0];
      tmp1 =1.0/sqrt(size3[0]*size3[0]+size3[1]*size3[1]+size3[2]*size3[2]);
      size3[0]*=tmp1;
      size3[1]*=tmp1;
      size3[2]*=tmp1;
      double projmax=0.0;
      double proj3;
      for (int i=0; i < numEdges; i++) {
         if((i!=ifirst)&&(i!=isecond)){ 
           tmp1=dotProd(edgesIonV[i],size3)*InvEdgeL[i];
           if(fabs(tmp1) > projmax) {
              projmax=abs(tmp1); // alignment (non-dim)
              proj3=tmp1/InvEdgeL[i]; // actual length projection onto size3
              ithird=i;
           }
         }
       }
// now scale the third vector
       size3[0]*=proj3;
       size3[1]*=proj3;
       size3[2]*=proj3;
       edgesIonV[ithird][0]=size3[0];
       edgesIonV[ithird][1]=size3[1];
       edgesIonV[ithird][2]=size3[2];
      dotProdTable[ithird][ithird]=proj3*proj3; // store size3's length here to be consistent with later use to get length of the third vector 
// now make second edge perpendicular to first and third
      size3[0]= edgesIonV[ithird][1]*edgesIonV[ifirst][2]
              - edgesIonV[ithird][2]*edgesIonV[ifirst][1];
      size3[1]= edgesIonV[ithird][2]*edgesIonV[ifirst][0]
              - edgesIonV[ithird][0]*edgesIonV[ifirst][2];
      size3[2]= edgesIonV[ithird][0]*edgesIonV[ifirst][1]
              - edgesIonV[ithird][1]*edgesIonV[ifirst][0];
// want square not this      tmp1 =1.0/sqrt(size3[0]*size3[0]+size3[1]*size3[1]+size3[2]*size3[2]);
      tmp1 =1.0/(size3[0]*size3[0]+size3[1]*size3[1]+size3[2]*size3[2]);
//project original v2 onto unit vector perp to V1 and v3.   
// to do this scale away its cross-product length to make it a unit vector then dot this vector with v2 and multiply by unit vector
// math v2<-(v2.size3/|size3|*size3/|size3| = (v2.size3)*size3/|size3|^2
// which is why above tmp1 has 1/|size3|^2 and next we *=tmp1 with the v2.size3
// to get the scaleing for size3
      tmp1*=dotProd(edgesIonV[isecond],size3);
       edgesIonV[isecond][0]=tmp1*size3[0];
       edgesIonV[isecond][1]=tmp1*size3[1];
       edgesIonV[isecond][2]=tmp1*size3[2];
       dotProdTable[isecond][isecond]=dotProd(edgesIonV[isecond],edgesIonV[isecond]);
// make second edge perpendicular to first and third

    } else {
// find the edge with the greatest orthonality to others
      for (int i=0; i < numEdges; i++) {
        for (int j=i+1; j < numEdges; j++) {
         tmp1=fabs(NormdotProdTable[i][j]);
            if(tmp1<projmin) {
               projmin=tmp1;
               ifirst=i;
               isecond=j;
            }
         }
      }
// find the third vector as the one with the weakest sum projection onto first and second
      projmin=2;
      ithird=-1;
      for (int i=0; i < numEdges; i++) {
         if((i!=ifirst)&&(i!=isecond)){ 
           tmp1=fabs(NormdotProdTable[ifirst][i]);
// average allows bad choice           tmp1+=fabs(NormdotProdTable[isecond][i]);
           tmp1=max(tmp1,fabs(NormdotProdTable[isecond][i]));
           if(tmp1 < projmin) {
              projmin=tmp1;
              ithird=i;
           }
         }
      }
      if(ithird==-1) {
         cout << "third edge failed" << endl;
      }

     } // should close not simmetrix suggested way
      eLength[0]=sqrt(dotProdTable[ifirst][ifirst]);    
      eLength[1]=sqrt(dotProdTable[isecond][isecond]);    
      eLength[2]=sqrt(dotProdTable[ithird][ithird]);    
      edgemin=1e6;
      edgemax=0.0;
      for (int k=0 ;k<3; k++) {
        if(eLength[k]>edgemax){
            edgemax=eLength[k];
            maxE=k;
        }
        if(eLength[k]<edgemin){
            edgemin=eLength[k];
            minE=k;
        }
      }
      for (int k=0 ;k<3; k++) if((k!=maxE) && (k!=minE)) midE=k;
      int d0,d1,d2;
      if(minE==0) d0=ifirst;
      if(minE==1) d0=isecond;
      if(minE==2) d0=ithird;
      if(midE==0) d1=ifirst;
      if(midE==1) d1=isecond;
      if(midE==2) d1=ithird;
      if(maxE==0) d2=ifirst;
      if(maxE==1) d2=isecond;
      if(maxE==2) d2=ithird;
      minE=d0;
      midE=d1;
      maxE=d2;
      double dots[3];
      dots[0]=fabs(dotProd(edgesIonV[minE],edgesIonV[midE]));
      dots[1]=fabs(dotProd(edgesIonV[maxE],edgesIonV[midE]));
      dots[2]=fabs(dotProd(edgesIonV[minE],edgesIonV[maxE]));
      double A20,A10,A21;
      eLength[0]=sqrt(dotProdTable[minE][minE]);    
      eLength[1]=sqrt(dotProdTable[midE][midE]);    
      eLength[2]=sqrt(dotProdTable[maxE][maxE]);    
      A20=eLength[2]/eLength[0];
      A10=eLength[1]/eLength[0];
      A21=eLength[2]/eLength[1];
//      if(coordvcur[0]<26) Isotrop=1;
//      if(A21<4 || A10 < 10) Isotrop=1; // 
      for (int k=0; k<3; k++) {
         if (dots[k]>0.2) {
           cout << "WARNING: two of the selected vectors are not orthogonal" << endl;
           Isotrop=1;
         }
      }
      if(minE==midE || minE==maxE || midE==maxE) {
         Isotrop=1;
         cout << "Setting Isotropic size on vertex " << icountVertsOnFaces << endl;
      }
     } // anisotrop      
     if(Isotrop==1) Isotrop=-1;
     if(Isotrop==1) {
         MSA_setVertexSize(simAdapter, 
                        vertex,
                        h[0]);
         icountIsotrop++;
     } else if (Isotrop==0){
          sizeRat=sizeRatio;
          double sizeRat1=1.0; 
          if(AnisoSimmetrix==-2) sizeRat1=sizeRatio;
//refine middle size if it will be larger than the largest after its split
          if(sizeRat1*eLength[1]>sizeRat*eLength[2]) 
             sizeRat1*=eLength[2]*sizeRat/eLength[1];
//   short edges cause trouble so rescale to length of middle edge
//   including  adustment down for sizeRat1
          double scl=sizeRat1*eLength[1]/eLength[0];
          OrgAnisoSize[0][0]= scl*edgesIonV[minE][0];
          OrgAnisoSize[0][1]= scl*edgesIonV[minE][1];
          OrgAnisoSize[0][2]= scl*edgesIonV[minE][2];
 
          OrgAnisoSize[1][0]= sizeRat1*edgesIonV[midE][0];
          OrgAnisoSize[1][1]= sizeRat1*edgesIonV[midE][1];
          OrgAnisoSize[1][2]= sizeRat1*edgesIonV[midE][2];
 
          OrgAnisoSize[2][0]= sizeRat*edgesIonV[maxE][0];
          OrgAnisoSize[2][1]= sizeRat*edgesIonV[maxE][1];
          OrgAnisoSize[2][2]= sizeRat*edgesIonV[maxE][2];
          int istop;
          if(coordvcur[1]<5e-6) {
             istop=1;
          }
          MSA_setAnisoVertexSize(simAdapter, 
                        vertex,
                       OrgAnisoSize);

          icountAnisotrop++;
     }

   } // the skip if not marked
  }  // iterator
  VIter_delete(vofIter);
  }  // iterator
  GFIter_delete(gfIter);
  } // region iterator...just a way to get region
  GRIter_delete(grIter);
//  5 lines above make a vertex iterator over model faces below is all
//  VIter_delete(vIter);
  fclose(itf);
  delete [] h;
//  if(PMU_rank()==0) {
    cout << "icountVertsNotBL " << PMU_rank() << " "  << icountVertsNotBL << endl;
    cout << "icountVertsOnFaces " << PMU_rank()  << " " << icountVertsOnFaces << endl;
    cout << "icountIsotrop " << PMU_rank()  << " " << icountIsotrop << endl;
    cout << "icountAnisotrop " << PMU_rank()  << " " << icountAnisotrop << endl;
 // }

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

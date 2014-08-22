#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "phastaIO.h"
#include <iostream>
#include "MeshSimInternal.h"
#include "mesh_interface.h"

#ifdef SUN
#include <cctype> // for tolower on linux
#endif

#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

extern char* iformat;
extern int adaptFlag;
extern int strategy;
using namespace std;
extern int multipleRestarts;
extern int numRestartsToRead;

extern int numVars;
extern pMeshDataId phasta_solution;
extern int numTotParts;
extern int numParts;

static int cscompare(const char* s1, const char* s2)
{
  while( *s1 == ' ') s1++;
  while( *s2 == ' ') s2++;
  while((*s1) && (*s2) &&(tolower(*(s1++))==tolower(*(s2++))))
    {
      while( *s1 == ' ') s1++;
      while( *s2 == ' ') s2++;
    }
  if (!(*s1) && !(*s2)) return 1;
  else return 0;
}

static void SwapArrayByteOrder(void* array, int nbytes, int nItems)
{
  int i,j;
  unsigned char  ucTmp;
  unsigned char* ucDst = (unsigned char*)array;

  for(i=0; i < nItems; i++) {
    for(j=0; j < (nbytes/2); j++)
      swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
    ucDst += nbytes;
  }
}
///////////////////////////////////////////////////////////////////////////////////////
//  restart only called if ICs are provided through restart files
//
//  We have zeroed the incoming qTot array because restart is going to fill
//  only the portion of this array that we request it to based on
//  our choice of nshgReqst and rRead (below called nvReqst for
//  consistency. Note that the array is shaped nv*nshgTot because
//  later, higher order modes will need to index into this TOTAL
//  solution array for all variables.  If they were not read it
//  will find zero's there which we assume will be ok or replaced
//  by attribute expressions. 
// 
// every proc carries the same copy of qTot
//
// called in writeEnsaFiles()
// 
// for multiple restart files only local solutions are read in (nshgTot --> nshg)
// 
// nvReqst <--> nvr
// nshgReqst <--> nshg
// passed out:
// *qTot
///////////////////////////////////////////////////////////////////////////////////////
void restart(double *qTot, int nshgTot, int nshgReqst, int nvReqst, int* lstep,
             pMesh mesh, int ipart)
{

    int nvr = nvReqst;
    int nshg = nshgReqst;
    int iarray[4];// contains: nshg,numVars,lstep
    // get the data either from file or from adaptor
    if(!adaptFlag  ){
        // convert from old version
        int restart;
        char rfname[40];
        sprintf(rfname,"restartc.inp.%d",PMU_rank()*numParts+ipart+1);
        
        openfile_( rfname, "read",  &restart );
        
        
        int isize = 3;
        
        readheader_( &restart, "solution", iarray,
                     &isize, "double", iformat );
        
        isize = iarray[0]*iarray[1];
        *lstep = iarray[2];
        double* qlocal = new double [ isize ];
        
        readdatablock_( &restart, "solution", qlocal, &isize,
                        "double" , iformat );
        
        // copy the needed part into the array passed in 
        if ( nshg > iarray[0] ) {
            cerr << "reading only " << iarray[0] << "modes from restart" << endl;
            cerr << nshg << " modes were requested " << endl;
            nshg = iarray[0];
        }
        if ( nvr  > iarray[1] ) {
            cerr << "reading only " << iarray[1] << "vars from restart" << endl;
            cerr << nvr << " modes were requested " << endl;
            nvr  = iarray[1];
        }
        // arrange qTot according to the nshg that are requested
        // (nshg may be smaller than the actual nshg(=iarray[0]) read in)
        // qTot carries "nvr" (=5) sequences of "nshg"-long sequences
        // loop num Variables
        for(int i=0; i < nvr; i++)
            // loop num shapefuns
            for(int j=0; j < nshg; j++)
                qTot[i*nshg+j] = qlocal[i*iarray[0]+j];
        
        delete [] qlocal;
    }//if(!adaptFlag)
    // mesh was previously adapted AND restart solution read in separately
    // NOT strategy 7: uniform splitting without restart files available
//    else if(adaptFlag &&  strategy!=7){
   else if(adaptFlag){      
        if(PMU_rank() == 0){
            cout<<"\n entering solution retriever from mesh: ADAPTED\n";
        }
        double* nodalSolution;
        pVertex vertex;
        VIter vIter = M_vertexIter( mesh );
        // variables written out
        int nshg_fine = M_numVertices(mesh);
        double* q = new double [ nshg_fine * nvr ];    
        int vcounter=0;
        
        while ( vertex = VIter_next( vIter ) ) {
            
	  if(EN_getDataPtr((pEntity)vertex, phasta_solution, 
			   (void**)&nodalSolution)){   
                
                int j=0;
                for( int i=0; i< nvr; i++ ) {
                    
                    //store the solution into q vector
                    //in the format [P(1),P(2).....P(nshg),u(1),....u(nshg),
                    //v(1),...v(nshg),w(1),......,T(1),.......]
                    q[vcounter+j] = nodalSolution[i]; 
                    j=j+nshg_fine;
//  #ifdef DEBUG
//                      if(PMU_rank() == 0){

//                          cout<<"\nin restart:retrieved a solution:"
//                              <<nodalSolution[i]<<"\n";
//                      }
//  #endif                
                }
                vcounter = vcounter+1;  
            }
            // there is no solution on this vertex
            else{
                cerr<<"\nerror in restart: encountered\n"
                    <<"a vertex that carries no solution (callback unsufficient)\n";
                
                V_info(vertex);
                SimPartitionedMesh_stop();
                exit(0);
            }    
        }// while ( vertex
        VIter_delete(vIter);    

        // arrange qTot according to the nshg that are requested
        // (nshg may be smaller than the actual nshg(=iarray[0]) read in)
        // qTot carries "nvr" (=5) sequences of "nshg"-long sequences
        // loop num Variables
        for(int i=0; i < nvr; i++){
            // loop num shapefuns
            for(int j=0; j < nshg; j++){
                qTot[i*nshg+j] = q[i*nshg_fine+j];

//  #ifdef DEBUG
//                  if(PMU_rank() == 0){
//                      cout<<"["<<PMU_rank()<<"]: restart - reading multiple restarts."
//                          <<" Local solu: q["<<j<<"]["<<i<<"]="<<qTot[i*nshg+j]<<"\n";
//                  }
//  #endif
            }
        }
        
        delete [] q;
    }//else if adaptFlag=1
}





#ifdef PARALLEL
/*******************************************************************
 * Compute the mesh partition information, ie which processor each
 * mesh region is on. And to attach the PID number to each region.
 *
 * There is also the constraint that the internal face and its SPEBC
 * slave lie on the same processor
 *
 * Chris Whiting (Fall 1997)
 *
 * Anilkumar Karanam (Fall 1998)
 *
 * J. Mueller, Spring 2007 - GEN-data functions replaced by STL maps
*******************************************************************
 */

#include <assert.h>
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"
#include "parallel.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#endif
/* global variables */
#include "phParAdapt.h"


extern int wGraph;
extern int justpart;
extern int parted;
extern map<pair<string,pGEntity>,int> GEntityDataContainerInt;
extern string intr;
extern string sper;

// name of fn is retained for historical reasons, but this is now just setting
// weights for the pmesh partitioning.
// called in  setup()
int partitionMesh(pParMesh pmesh)
{
    
  if(justpart == 1){
      if(PMU_rank()==0){

          pMesh mesh = PM_mesh(pmesh, 0);

          // UserWeight not used at the moment
          UserWeight(mesh);  // setting weights according to ndof
          
          pGModel model = M_model(pmesh);

          pGFace gface[2];
          
          GFIter gfiter = GM_faceIter(model);
          
          // intr set in setPeriodic
          // gface[0] is kept as a non-NULL only when it is an interior model face
          int tmp;

          while (gface[0] = GFIter_next(gfiter))
              if (GEN_dataI((pGEntity)gface[0],"intr",&tmp)) break;// interior model face

          GFIter_reset(gfiter);
          
          while (gface[1] = GFIter_next(gfiter))  {
              if (!GEN_dataI((pGEntity)gface[1],"intr",&tmp))// interior model face
                  if(GEN_dataI((pGEntity)gface[1],"sper",&tmp)) break;//SPEBC
          }
          
          GFIter_delete(gfiter);
          
          // I dont understand this.. cant there be 2 or more intr-sper pairs? why
          // should gface[0..1] be a pair? Also cant gface[1] exist and not gface[0]?

          // usef for weighting
          pMeshDataId id;
          
          if (gface[1])  {
              
              int i, n, Sregion = 0;
              
              RIter rIter = M_regionIter(mesh);
              
              pRegion region, savereg = 0;
              
              // loop mesh regions
              while (region = RIter_next(rIter))  {
                  
                  pPList verts = return_vertices(region);
                  
                  n = PList_size(verts);
                  
                  // create a weight set for each region

                  id = PM_newRegionWtId ("defaultWtSet");
                  
                  // loop vertices
                  for (i=0; i < n; i++)  {
                      
                      pVertex vertex = (pVertex) PList_item(verts, i);
                      
                      pGEntity gent = V_whatIn(vertex);
                      
                      // gface[0] could be NULL ???
                      if(gface[0]){
                          if (GF_inClosure(gface[0],gent) || GF_inClosure(gface[1],gent))
                              break;
                      }
                  }
                  
                  PList_delete(verts);
                  
                  if (i < n)  {  // vertex was on gface[0] or gface[1]
                      if (!savereg)
                          savereg = region;
                      else
                          
                          R_setWeight(region, id,-1); // negative wt => zero wt  
                      Sregion++;
                      
                  }// if (i < n) 

              }// loop regions
              
              assert(savereg != 0);
              R_setWeight(savereg, id, Sregion);
              
          }//if (gface[1]) 
      }
  }//justpart==1
  // on 0  procs if in  1stage - mesh is NOT parted
  else if (justpart ==0 && parted == 0){
      if(PMU_rank()==0){
          pMesh mesh = PM_mesh(pmesh, 0);
          
          // not used at the moment
          UserWeight(mesh);  // setting weights according to ndof
          
          pGModel model = M_model(pmesh);
          
          pGFace gface[2];
          
          GFIter gfiter = GM_faceIter(model);
          
          int tmp;

          while (gface[0] = GFIter_next(gfiter))
              if (GEN_dataI((pGEntity)gface[0],"intr",&tmp)) break;
          
          GFIter_reset(gfiter);
          
          while (gface[1] = GFIter_next(gfiter))  {
              if (!GEN_dataI((pGEntity)gface[1],"intr", &tmp))
                  if(GEN_dataI((pGEntity)gface[1],"sper",&tmp)) break;
          }
          
          GFIter_delete(gfiter);
      
          // I dont understand this.. cant there be 2 or more intr-sper pairs? why
          // should gface[0..1] be a pair? Also cant gface[1] exist and not gface[0]?
          
          // usef for weighting
          pMeshDataId id;
          
          if (gface[1])  {
              
              int i, n, Sregion = 0;
              
              RIter rIter = M_regionIter(mesh);
              
              pRegion region, savereg = 0;
              
              // loop mesh regions
              while (region = RIter_next(rIter))  {
                  
                  pPList verts = return_vertices(region);
                  
                  n = PList_size(verts);
              
                  // create a weight set for each region
                  
                  id = PM_newRegionWtId ("defaultWtSet");
                  
                  // loop vertices
                  for (i=0; i < n; i++)  {
                      
                      pVertex vertex = (pVertex) PList_item(verts, i);
                      
                      pGEntity gent = V_whatIn(vertex);
                      
                      // gface[0] could still be 0 ==> segfault
                      if(gface[0]){
                          if (GF_inClosure(gface[0],gent) || GF_inClosure(gface[1],gent))
                              break;
                      }
                  }
                  
                  PList_delete(verts);
                  
                  if (i < n)  {  // vertex was on gface[0] or gface[1]
                      if (!savereg)
                          savereg = region;
                      else
                          R_setWeight(region, id,-1); // negative wt => zero wt  

                      Sregion++;
                      
                  }// if (i < n) 
                  
              }// loop regions
              
              assert(savereg != 0);
              R_setWeight(savereg, id, Sregion);
              
          }//if (gface[1]) 
          
      }//if (PMU_rank() == 0)  
  }//if (justpart == 0)  
  
  // I hope above is right and works. Anyway, cannot make the spebc regions lie
  // on zeroth proc.. code should be made to work for spebc regions on any proc

  return 1; // due to userweight, always check weights for partitioning
}

#endif

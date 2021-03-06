#include <iostream>
#include <stdlib.h>
#include "InitialCondition.h"
#include "ccfunc.h"
#include <stdlib.h>
// #ifndef SIM
// #include <vector>
// #include "AttributeException.h"
// using namespace SCOREC_Util;
// #endif


using namespace std;

extern int ensa_dof, rRead;

// Right now user must be very careful to see that what he is reading
// out of the restart file complemets what is set through the GUI

// If anything is read twice , the GUI takes precedence

InitialCondition::InitialCondition(pGModel model)
{
  static int eflg[3] ={ 1, 1 ,1};
  static int stick=1;
  char* flag[7]={"initial pressure", "initial velocity", "initial temperature",
                 "initial scalar_1", "initial scalar_2", "initial scalar_3",
                 "initial scalar_4" };
//#ifdef SIM
  // put it in a loop... - SST
  attList[0] = GM_attrib(model, "initial pressure");
  attList[1] = GM_attrib(model, "initial velocity");
  attList[2] = GM_attrib(model, "initial temperature");
  attList[3] = GM_attrib(model, "initial scalar_1");
  attList[4] = GM_attrib(model, "initial scalar_2");
  attList[5] = GM_attrib(model, "initial scalar_3");
  attList[6] = GM_attrib(model, "initial scalar_4");
// #else
//   std::vector<Attribute*> atts;
//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial pressure");
//     attList[0] = atts[0];
//   } catch  (AttributeNotExistent) { attList[0] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial velocity");
//     attList[1] = atts[0];
//   } catch  (AttributeNotExistent) { attList[1] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial temperature");
//     attList[2] = atts[0];
//   } catch  (AttributeNotExistent) { attList[2] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial scalar_1");
//     attList[3] = atts[0];
//   } catch  (AttributeNotExistent) { attList[3] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial scalar_2");
//     attList[4] = atts[0];
//   } catch  (AttributeNotExistent) { attList[4] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial scalar_3");
//     attList[5] = atts[0];
//   } catch  (AttributeNotExistent) { attList[5] = 0 ; }

//   try{
//     atts = SCOREC_att::retrieveAttributePList((pGEntity)model,"initial scalar_4");
//     attList[6] = atts[0];
//   } catch  (AttributeNotExistent) { attList[6] = 0 ; }

// #endif

  for(int i =0; i< ensa_dof-2; i++) if (attList[i]){
    if(eflg[i]) {
       if (PMU_rank() == 0)
         cout<< flag[i] <<" has been set in the GUI "<<endl;
       eflg[i] = 0;
    }
  }

  if ((!attList[0] || !attList[1] || !attList[2]) && (rRead == 0)) {
    cerr << "\nNSpre Error: you must specify an initial condition or use "
         << "the restart option" << endl;
    exit(-1);
  }
  if (stick && PMU_rank()==0){
     cout <<" The Rest of the Initial Conditions (if any left) "<<endl
          <<" Should be coming from the restart file " << endl;
     stick = 0;
  }
}

void InitialCondition::eval(pVertex vertex, double *q)
{
  double x[3];
  int nsc;
  nsc=ensa_dof-5;
  V_coord(vertex,x);

//old
//    if(attList[0]) q[0] = Att_evalTensorOr0(attList[0], x);
//    if(attList[1]) Att_evalTensorOr1(attList[1], x, q+1);
//    if(attList[2]) q[4] = Att_evalTensorOr0(attList[2], x);
//    for (int i=0; i < nsc; i++)
//      if(attList[3+i]) q[5+i] = Att_evalTensorOr0(attList[3+i], x);
  if(attList[0]) q[0] = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x);
  if(attList[1]) {
    q[1] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 0, x);
    q[2] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 1, x);
    q[3] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 2, x);
  }
  if(attList[2]) q[4] = AttributeTensor0_evalDS((pAttributeTensor0) attList[2], x);
  for (int i=0; i < nsc; i++)
    if(attList[3+i]) q[5+i] = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x);
}

void InitialCondition::eval(pEdge edge, double *q, int p)
{
  int nsc = ensa_dof-5;
  double factor = 1.00000000000;

  if (p > 2){
    for (int i=0; i < 5+nsc; i++) q[i]=0.0;
    return;
  }

  pVertex v1 = E_vertex(edge,0);
  pVertex v2 = E_vertex(edge,1);

  double x1[3], x2[3];
  V_coord(v1, x1);
  V_coord(v2, x2);
  double x3[3] = {0.5*(x1[0]+x2[0]), 0.5*(x1[1]+x2[1]), 0.5*(x1[2]+x2[2])};
  double d1,d2,d3,vel1[3],vel2[3],vel3[3];

  if(attList[0]) {
//old
//      d1 = Att_evalTensorOr0(attList[0], x1);
//      d2 = Att_evalTensorOr0(attList[0], x2);
//      d3 = Att_evalTensorOr0(attList[0], x3);
//      q[0] = factor*(d1 + d2 - 2.0 * d3);
    d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x1);
    d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x2);
    d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x3);
    q[0] = (d1 + d2 - 2.0 * d3);
  }

  if(attList[1]){
//old
//      Att_evalTensorOr1(attList[1], x1, vel1);
//      Att_evalTensorOr1(attList[1], x2, vel2);
//      Att_evalTensorOr1(attList[1], x3, vel3);
    for(int i=0; i<3; i++) {
      vel1[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x1);
      vel2[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x2);
      vel3[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x3);
    }
    q[1] = factor*(vel1[0] + vel2[0] - 2.0 * vel3[0]);
    q[2] = factor*(vel1[1] + vel2[1] - 2.0 * vel3[1]);
    q[3] = factor*(vel1[2] + vel2[2] - 2.0 * vel3[2]);
  }

  if(attList[2]){
//old
//      d1 = Att_evalTensorOr0(attList[2], x1);
//      d2 = Att_evalTensorOr0(attList[2], x2);
//      d3 = Att_evalTensorOr0(attList[2], x3);
    d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x1);
    d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x2);
    d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x3);
    q[4] = factor*(d1 + d2 - 2.0 * d3);
  }

  for (int i=0; i < nsc && attList[3+i] ; i++) {
    if (attList[3+i]){
//old
//        d1 = Att_evalTensorOr0(attList[3+i], x1);
//        d2 = Att_evalTensorOr0(attList[3+i], x2);
//        d3 = Att_evalTensorOr0(attList[3+i], x3);
      d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x1);
      d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x2);
      d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x3);
      q[5+i] =factor*(d1 + d2 - 2.0 * d3);
    }
  }
}

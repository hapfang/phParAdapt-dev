#ifndef _H_shapeFuncInternals
#define _H_shapeFuncInternals

#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

int C_equal1(double real1,double real2);

double En(int ip, double r, double s);
int EnDrv(int ip, double r, double s, double *drv);
int EnDrv2(int ip, double r, double s, double drv[]);
double Fn(int i, int j, double r, double s);
int FnDrv(int i,int j, double r, double s, double drv[2]);
int FnDrv2(int i,int j, double r, double s, double drv[]);
double Bn(int i, int j, int k, double r, double s, double t);
int BnDrv(int i,int j, int k, double r, double s, double t, double drv[3]);
int BnDrv2(int i,int j, int k, double r, double s, double t, double drv[]);

double V_blendOnEntity(pVertex v, pEntity e, double *L);
double V_blendIndexed(int i, double *L);
double V_blendIndexedOnEdge(int i, double *L);
double E_blendOnFace(pEdge edge, pFace face, double *L);
double F_edgeBlendTri(int *index, double *L);
double F_edgeBlendQuad(int *index, double *L);
double E_blendOnRegion(pEdge edge, pRegion region, double *L);
double R_edgeBlendTet(int *index, double *rstw);
double F_blendOnRegion(pFace face, pRegion region, double *L);


double E_modeShape(int p, double *L);
double F_modeShapeTri(int p, int i, double *L);
double R_modeShapeTet(int p, int i, double *L);
int E_modeShapeDrv(int p, double *L, int,double drv[2]);
int F_modeShapeTriDrv(int p, int i, double *L, int,double mdrv[2]);
double F_modeShapeQuad(int p, int i, double *L);
double R_modeShapeHex(int p, int i, double *L);
int F_modeShapeQuadDrv(int p, int i, double *L, int,double mdrv[]);
int R_modeShapeTetDrv(int p, int i, double *L, int,double mdrv[3]);
int R_modeShapeHexDrv(int p, int i, double *L, int,double mdrv[]);

int V_index(pVertex v, pEntity ent, int *index);
int E_index(pEdge e, pEntity ent, int *index) ;
int F_index(pFace face, pEntity ent, int *index);
int V_blendOnEntityDrv(pVertex v, pEntity e, double *L, double mdrv[3]);
int E_blendOnFaceDrv(pEdge edge, pFace face, double *L, double bdrv[2]);
int F_edgeBlendTriDrv(int *, double *L, double drv[]);
int E_blendOnRegionDrv(pEdge edge, pRegion region, double *L, double bdrv[3]);
int F_blendOnRegionDrv(pFace face, pRegion region, double *L, double drv[3]);
int R_edgeBlendTetDrv(int *, double *L, double drv[]);
int R_faceBlendTetDrv(int *index, double *L, double drv[]);

int E_parDrv(int i, int j, pEntity elem, double drv[][2]);
int F_parDrv(int i, int j, int k, pEntity elem, int (**drv)[2]);

int V_blendOnEntityDrv2(pVertex v, pEntity e, double *L, double mdrv[]) ;
int E_blendOnFaceDrv2(pEdge edge, pFace face, double *L, double bdrv[2]) ;
int F_edgeBlendTriDrv2(int index[2], double *L, double drv[]) ;
int E_blendOnRegionDrv2(pEdge edge, pRegion region, double *L, double bdrv[3]) ;
int F_blendOnRegionDrv2(pFace face, pRegion region, double *L, double drv[3]) ;
int R_edgeBlendTetDrv2(int *, double *L, double drv[]);
int R_faceBlendTetDrv2(int *index, double *L, double drv[]);
#ifdef __cplusplus
}
#endif

#endif


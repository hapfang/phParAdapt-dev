#ifndef _H_shapeFuntion
#define _H_shapeFuntion

#ifdef __cplusplus
extern "C" {
#endif

int Entity_setP(pEntity entity, char *field, int p);
int Entity_queryP(pEntity entity, char *field);
int Entity_shapeMinPolyOrder(pEntity entity);
int Entity_numShapeFunction(pEntity entity, int porder);
double Entity_shapeFunction(pEntity entity, pEntity element, int p, int ith, double *L);
int Entity_shapeFuncDrv(pEntity entity, pEntity element, int p, int ith, 
                            double *L, double *shpDrv);
int Entity_shapeFuncDrv2(pEntity entity, pEntity element, int p, int ith, 
                            double *L, double *shpDrv);


/* Functions to generate the Legandre Polynomials and their derivatives */

 double LP(int j, double x);
 double LPdrv(int j, double x);
 double phi(int p, double x);
 double phiDrv(int p,double x);
 int HexShapeAndDrv(int p, double par[4], double N[], double
			      dN[][3]);

/* Blending functions */

 double Line_eB(double xi1);
 double dLEBdxi1(double xi1);
 double dLEBdxi2(double xi1);
 double dLEBdxi3(double xi1);

 double Quad_eB(double xi1, double xi2, int sign);
 double dQEBdxi1(double xi1, double xi2, int sign);
 double dQEBdxi2(double xi1, double xi2, int sign);
 double dQEBdxi3(double xi1, double xi2, int sign);

 double Quad_fB(double xi1, double xi2);
 double dQFBdxi1(double xi1, double xi2);
 double dQFBdxi2(double xi1, double xi2);
 double dQFBdxi3(double xi1, double xi2);

 double Hex_eB(double xi[3], int sign2, int sign3);
 double dHEBdxi1(double xi[3], int sign2, int sign3);
 double dHEBdxi2(double xi[3], int sign2, int sign3);
 double dHEBdxi3(double xi[3], int sign2, int sign3);

 double Hex_fB(double xi[3], int sign3);
 double dHFBdxi1(double xi[3], int sign3);
 double dHFBdxi2(double xi[3], int sign3);
 double dHFBdxi3(double xi[3], int sign3);

/* Entity Level functions */

 int mesh_edge(double xi1,int gOrd[3], int p, double* entfn,
			 double** edrv);
 int quad_face(double xi[3], int gOrd[3], int p, double*
			 entfn, double** edrv); 
 int hex_regn(double xi[3], int p, double*
			entfn, double** edrv);


#ifdef __cplusplus
}
#endif

#endif

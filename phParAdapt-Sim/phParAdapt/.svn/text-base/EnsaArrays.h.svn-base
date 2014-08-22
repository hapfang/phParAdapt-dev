
#ifndef H_EnsaArrays
#define H_EnsaArrays

class EnsaParameters;

class EnsaArrays {
public:
  EnsaArrays(globalInfo *info, const EnsaParameters &par);
  ~EnsaArrays();
  void write(pMesh mesh, globalInfo *info, const EnsaParameters &par,
             double zscale[3], int *ifath, int* nsons, int ipart);

  int ***getIEN() { return ien; }
  int ***getIENB() { return ienb; }
  int *getNBC() { return nBC; }
  int *getIBC() { return iBC; }
  int ***getIBCB() { return iBCB; }
  int *getIPER() { return iper; }
  int *getILWORK() { return ilwork; }
  int *getNCORP() { return ncorp; }
  double **getX() { return x; }
  double **getQ() { return q; }
  double ***getBCB() { return BCB; }
  double **getBC() { return BC; }

private:
  globalInfo *_info;
  void fixInd(globalInfo *info);
  int ***ien, ***ienb, *nBC, *iBC, ***iBCB, *iper, *ilwork, *ncorp;
  double **x, **q, ***BCB, **BC;
  int nshp;
};

#endif

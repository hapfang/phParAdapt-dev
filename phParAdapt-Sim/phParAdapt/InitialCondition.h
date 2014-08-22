//
// Initial condition class
//
#ifdef FMDB
#include "myAttribute.h"
#endif
#include "parallel.h"
class InitialCondition {
public:
  InitialCondition(pGModel );
  void eval(pVertex , double *);
  void eval(pEdge , double *, int);
private:
  pAttribute attList[7]; //  two thermo + velocity vector (3) + 4 scalars
};

/*----------------------------------------------------------------------------
   returns TRUE if two real numbers are equals within machine precision
----------------------------------------------------------------------------*/
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ZERO 1.0e-14

int C_equal1(double real1,double real2)
{
   double areal1,areal2,rldif,rltol ;

   areal1 = fabs(real1) ;
   areal2 = fabs(real2) ;
   rldif = fabs(real1-real2) ;
   if ( areal1 > 0. && areal2 > 0. )
   {
     if ( areal2 > areal1 )
       rltol = rldif / areal1 ;
     else
       rltol = rldif / areal2 ;
     if ( rltol > ZERO )
       return (0) ;
     else
       return (1) ;
   }
   else
   {
     if ( rldif > ZERO )
       return (0) ;
     else
       return (1) ;
   }
}   

#ifdef __cplusplus
}
#endif

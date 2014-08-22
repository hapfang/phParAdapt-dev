#include <string.h>
#include "parallel.h"

/* This function generates the descriptive keyphrase for each topology
   polynomial order block based on the block they are in
   the standard key phrase prefix is "connectivity interior" for the interior
   elements and "connectivity boundary" for the boundary elements.
   example: cubic trifaces pyramids on the boundary get the phrase
   "connectivity boundary cubic pyramid_triface"
*/

#ifdef __cplusplus
extern "C" {
#endif

//old
void generate_keyphrase(char* target, const char* prefix,
                        const blockKey* tpblock)
{
  strcpy(target,"connectivity ");
  strcat(target,prefix);  /* interior or boundary */

  switch(tpblock->maxpoly){
  case 1:
    strcat(target,"linear ");
    break;
  case 2:
    strcat(target,"quadratic ");
    break;
  case 3:
    strcat(target,"cubic ");
    break;
  case 4:
    strcat(target,"quartic ");
    break;
  }

  switch(tpblock->lcsyst){
  case 1:
    strcat(target,"tetrahedron ");
    break;
  case 2:
    strcat(target,"hexahedron ");
    break;
  case 3:
    if(!strcmp(prefix,"boundary "))
      strcat(target,"wedge triface ");
    else
      strcat(target,"wedge ");
    break;
  case 4:
    strcat(target,"wedge quadface ");
    break;
  case 5:
    if(!strcmp(prefix,"boundary "))
      strcat(target,"pyramid quadface ");
    else
      strcat(target,"pyramid ");
    break;
  case 6:
    strcat(target,"pyramid triface ");
    break;
  }
}
void generate_keyphraseNew(char* target, char* prefix, blockKey* tpblock) {

  strcpy(target,prefix);  /* interior or boundary */

  switch(tpblock->maxpoly){
  case 1:
    strcat(target,"linear ");
    break;
  case 2:
    strcat(target,"quadratic ");
    break;
  case 3:
    strcat(target,"cubic ");
    break;
  case 4:
    strcat(target,"quartic ");
    break;
  }

  switch(tpblock->lcsyst){
  case 1:
    strcat(target,"tetrahedron ");
    break;
  case 2:
    strcat(target,"hexahedron ");
    break;
  case 3:
    if(!strcmp(prefix,"boundary "))
      strcat(target,"wedge triface ");
    else 
      strcat(target,"wedge "); 
    break;
  case 4:
    strcat(target,"wedge quadface ");
    break;
  case 5:
    if(!strcmp(prefix,"boundary "))
      strcat(target,"pyramid quadface ");
    else 
      strcat(target,"pyramid ");
    break;
  case 6:
    strcat(target,"pyramid triface ");
    break;
  }

}

#ifdef __cplusplus
}
#endif

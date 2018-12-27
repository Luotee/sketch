/********************************************************************
Approximate frequent items in a data stream
G. Cormode 2002

Last modified: 2003-10

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prng.h"
#include "massdal.h"

/******************************************************************/

#include "countmin.h"

/******************************************************************/



CM_type * cm;

int main(int argc, char **argv) 
{
  int i, uptime, outtime; 
  int width=4, depth=5;
  int * stream;

  cm=CM_Init(width,depth);

  CM_Update(cm, 12, 5);
  CM_Update(cm, 19, 2);
  CM_Update(cm, 7, 7);

  printf("Estimate of 12 is %d\n",CM_PointEst(cm,12));
  printf("Estimate of 12 is %d\n",CM_PointMed(cm,12));
  printf("Estimate of 19 is %d\n",CM_PointEst(cm,19));
  printf("Estimate of 19 is %d\n",CM_PointMed(cm,19));
  printf("Estimate of 7 is %d\n",CM_PointEst(cm,7));
  printf("Estimate of 7 is %d\n",CM_PointMed(cm,7));

  printf("Space was %d\n", CM_Size(cm));
  CM_Destroy(cm);
  return 0;
}

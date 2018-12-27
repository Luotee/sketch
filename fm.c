/********************************************************************
Flajolet-Martin Distinct Counting
Count distinct elements, G. Cormode 2001-2004

The current version does inserts only. 
Inserts and deletes is also possible with a bit more work. 

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "prng.h"
#include "fm.h"


///////////////////////////////////////////////////////////////////////////
//////  Count Distinct Routines

#define distinct_width 64
int distinct_width_mod = distinct_width -1;

unsigned int approx_distinct_zeros_table[16] = 
  {4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
// used to count the number of trailing zeros in the binary representation
// of a hash value

FM_type * FM_Init(int fmsize, int seed) {
  // Seed is a seed to pass to the random number generator
  // fmsize is the number of repetitions to try


  // this creates a bunch of hash functions 
  int i;
  FM_type * result;
  prng_type * prng;

  prng=prng_Init(seed,2); // fix 
  result=(FM_type *) calloc(1,sizeof(FM_type));
  result->fmsize=fmsize;
  result->fm=calloc(fmsize,sizeof(unsigned int));
  result->hasha=calloc(fmsize,sizeof(unsigned int));
  result->hashb=calloc(fmsize,sizeof(unsigned int));

  for (i=0;i<fmsize; i++)
    {
      result->fm[i]=0;
      result->hasha[i]=prng_int(prng);
      result->hashb[i]=prng_int(prng);
    }
  prng_Destroy(prng);
  return (result);
}

int zeros_slow(unsigned int test) {

  /*****************************************************************/
  /* For use in the hash function, return the number of trailing   */
  /* zeros in the binary representation of test.  This is not a    */
  /* very efficient implementation, but it'll do for now           */
  /*****************************************************************/

  int i;
  int bits=32;

  if (test == 0) return bits;
  for (i=0; i<bits; i++) {
    if ((test & 1) == 1) break;  
    // stop when lowest bit is 1, else shift right & continue
    test >>=1;
  }
    // this doesn't have to be very fast unless you do this a lot
  return i;
}

int zeros(unsigned int test) {
  // count the number of zeros in the binary representation of a number
  int result=0;

  if (test == 0) return 0; // special case, unlikely to happen
  while ((test & 15) == 0)
    {
      result+=4; test>>=4;
    } // process the number a nibble at a time
  return result+approx_distinct_zeros_table[test&15];
} // look up the answer for the last part

void FM_Update(FM_type * fm, unsigned int item)
{
  int i;
  unsigned int hash;
  
  for (i=0;i<fm->fmsize;i++)
    {
      hash=hash31(fm->hasha[i],fm->hashb[i],item);
      fm->fm[i]|=(1<<zeros(hash));
    } // update each entry with the number of zeros in the hash of item
}

double FM_Distinct(FM_type * fm)
{

  int bits=distinct_width;
  int max, i, j;
  unsigned int bitmap;
  double result=0.0;

  for (i=0; i<fm->fmsize; i++) {
    max=bits;
    bitmap=fm->fm[i];
    for (j=bits-1; j>=0; j--) {
      if (((bitmap>>j) & 1)==0) max=j;
    }
    // for each of the vectors, find the position of the smallest
    // zero entry
    result=result+ ((float) max);
    // add this to the running total
  } 
  // find the average
  result = result /((double) fm->fmsize);
  result = 1.2928*pow(2.0,result);
  // scale with the factor given in FM85
  return (result);
}

void FM_Destroy(FM_type * fm){

  if (!fm) return;
  free(fm->hashb);
  free(fm->hasha);
  free(fm->fm);
  free(fm);
  fm=NULL;
} 

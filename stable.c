/****************************************************************
 * Sketches of vectors using stable distributions 
 * For Euclidean, Manhattan, L_p distances
 * File originally Gaussian Distribution generator, Nick Koudas
 * Now extended to pick from arbitrary stable distributions
 * with alpha in the range (0, 2], Graham Cormode 
 * File date: 2002-03-13 (ish)
 * Updated:  2003-12-20
 ****************************************************************/

// This is based on using double floating point types
// for extra speed, you can use floats instead

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prng.h"
#include "stable.h"
#include "massdal.h"

/**********************************************************************/
/* Finally, we get down to business: some exported routines.          */
/* makesk given a vector will produce the sketch of it                */
/* norm  given a sketch will give you the l_p norm of it              */
/* dist probably shouldn't be used at the moment, it's rather hacky   */
/* sumsketches adds sketches, and diffsketches substracts them...     */
/* So to find the L_2 difference of two vectors represented by        */
/*  sketches then call norm(diffsketches(sketch1, sketch2));          */
/**********************************************************************/

double Stable_norm(Stable_sk * sket) {

  /***************************************************************/
  /* Find the L_p norm of the vector that a sketch represents,   */
  /* and return this to the power p -- needed for distinct values*/
  /* and so on                                                   */
  /***************************************************************/

  double* holder;
  double sum=0.0;
  double est;
  int i;

  // use the j-l lemma approach for L_2 distance 
  if (sket->alpha==2.0) 
    {
      for (i=0; i<sket->sksize; i++) 
	{
	  est=sket->sk[i]*sket->sk[i];
	  sum+=est;
	  //printf("Estimator is %f \n",est);
	}
      sum = pow(sum/((double) sket->sksize),0.5);
      // Calculate the L_alpha norm 
    } 
  else 
    {
      holder=(double*) calloc(sket->sksize+1, sizeof(double));
      for (i=0; i<sket->sksize; i++) 
	holder[i+1]=fabs(sket->sk[i]); 
      // transfer the details into a suitable array 
      sum=DMedSelect(sket->sksize/2,sket->sksize,holder);
      //find the median of the arary, this is the estimator for the 
      // L_p norm of the vector
      free(holder);
      if (sket->alpha<0.01)
	sum=pow(sum,0.02);
      else
	sum=pow(sum,sket->alpha);
      //return not the L_p norm but the L_p norm ^p  
    }  
  return sum;
}

int Stable_comparable(Stable_sk * sk1, Stable_sk * sk2) {

  if (sk1->alpha!=sk2->alpha) // incomparable sketches
    return 0;
  if (sk1->sksize!=sk2->sksize) // incomparable sketches
    return 0;
  if (sk1->seed!=sk2->seed) 
    return 0;
  return 1;
}

double Stable_dist(Stable_sk * s1, Stable_sk * s2) {

  /**********************************************************************/
  /* Calculate the distance between two sketches s1 and s2              */
  /* Do this by treating them by first finding the absolute difference  */
  /* between each component.  The either use median (for alpha<2) or L2 */
  /* (for alpha=2) to approximate the L_alpha distance.                 */
  /**********************************************************************/

  int i;
  double sum = 0.0;
  float alpha;
  double* holder;

  if(Stable_comparable(s1,s2)!=1) return -1.0;

  alpha=s1->alpha;
  if (alpha==2.0) {
    for (i = 0; i < s1->sksize; i++) {
      sum+= pow(fabs((s1->sk[i]-s2->sk[i])), 2.0);    
    }
    sum= sum/(double) s1->sksize;
    sum=pow(sum,0.5);
  // Calculate the L_alpha distance
  } else {
    // alternate method:
    holder=(double *) calloc(s1->sksize+1,sizeof(double *));
    for (i=0; i<s1->sksize; i++) 
      holder[i+1]=(s1->sk[i]-s2->sk[i]);
    sum=DMedSelect(s1->sksize/2,s1->sksize,holder);
    free(holder);
  }
  return sum;
}

/**********************************************************************/

Stable_sk * Stable_Init(int sksize, double alpha, long masterseed) {

  Stable_sk * result;

  result=(Stable_sk *) calloc(sizeof(Stable_sk),1);

  result->sksize=sksize;
  result->sk=(double *) calloc(result->sksize,sizeof(double));
  result->alpha=alpha;

  result->seed=masterseed;
  // firstly, set up the parameters of the sketch, zero the entries
  result->prng=prng_Init(masterseed,2);

  return (result); 
}

void Stable_Update(Stable_sk * sk, int item, double val) {
  // update the sketch with an addition of val to item 

  int j; 
  long seed;
 
  seed = (item+sk->seed);
  prng_Reseed(sk->prng,seed);
  // for each entry in the vector, create a pseudo-random sequence 
  // based on its position
  for (j=0; j<sk->sksize; j++) {
      sk->sk[j]+=val*prng_stable(sk->prng, sk->alpha);
      // use this to make a sequence of stable variables and 
      // multiply these by the entry in the vector to maintain the sketch
    }
}  

Stable_sk * Stable_makesk(double * vect, int length, int sksize, 
			  double alpha, long masterseed) {

  /********************************************************************/
  /* make a sketch of a vector vect which has length length           */
  /* alpha is the l_p norm to use, masterseed should be some positive */
  /* integer to initiate the process.  Returns a sketch for vect      */
  /* different sketches must have the same alpha and masterseed if    */
  /* they are to be comparable                                        */
  /********************************************************************/

  int i;
  Stable_sk * result;

  result=Stable_Init(sksize,alpha,masterseed);
  for (i=0;i<length;i++) {
    Stable_Update(result,i,vect[i]);
  }    
  return result;
}

void Stable_AddSketch(Stable_sk * s1, Stable_sk * s2) {

  // if s1 represents a vector v1 and s2 represents a vector v2 then
  // sumsketches(s1,s2) returns a sketch that represents the vector
  // v1 + v2

  int i;

  if(Stable_comparable(s1,s2)==1) {
    for (i=0; i<s1->sksize; i++)
      s1->sk[i]+=s2->sk[i];
  }
}

void Stable_SubSketch(Stable_sk * s1, Stable_sk * s2) {

  // if s1 represents a vector v1 and s2 represents a vector v2 then
  // sumsketches(s1,s2) returns a sketch that represents the vector
  // v1 + v2

  int i;

  if(Stable_comparable(s1,s2)==1) {
    for (i=0; i<s1->sksize; i++)
      s1->sk[i]-=s2->sk[i];
  }
}

void Stable_Destroy(Stable_sk * sk)
{
  if (sk)
    {
      free(sk->prng);
      free(sk);
    }
}

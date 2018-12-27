/********************************************************************
AMS Sketches
G. Cormode 2003

This version: 2003-12

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA. 
*********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "ams.h"
#include "prng.h"
#include "massdal.h"

AMS_type * AMS_Init(int buckets, int depth)
{
  int i,j;
  AMS_type * result;
  prng_type * prng;

  prng=prng_Init(-6321371,2);

  result=calloc(1,sizeof(AMS_type));
  if (result==NULL) exit(1);
  result->depth=depth;
  result->buckets=buckets;
  result->count=0;
  for (i=0;i<6;i++)
    result->test[i]=calloc(depth,sizeof(long long));
  // create space for the hash functions

  result->counts=(int *) calloc(buckets*depth, sizeof(int));
  if (result->counts==NULL) exit(1); 

  for (i=0;i<depth;i++)
    {
      for (j=0;j<6;j++)
	{
	  result->test[j][i]=(long long) prng_int(prng);
	  if (result->test[j][i]<0) result->test[j][i]= -result->test[j][i];
	  // initialise the hash functions
	  // prng_int() should return a random integer
	  // uniformly distributed in the range 0..2^31
	}
    }
  prng_Destroy(prng);
  return (result);
}

void AMS_Update(AMS_type * ams, unsigned long item, int diff)
{
  // update the sketch
  // hash to one bucket in each row
  // then multiply by {+1, -1} chosen at random

  int j;
  unsigned int hash;
  int mult, offset;

  ams->count+=diff;
  offset=0;
  for (j=0;j<ams->depth;j++)
    {
      hash=hash31(ams->test[0][j],ams->test[1][j],item);
      hash=hash % (ams->buckets); 
      mult=fourwise(ams->test[2][j],ams->test[3][j],
		    ams->test[4][j],ams->test[5][j],item);
      if ((mult&1)==1)
	ams->counts[offset+hash]+=diff;
      else
	ams->counts[offset+hash]-=diff;
      offset+=ams->buckets;
    }
}


int AMS_Compatible(AMS_type * a, AMS_type * b){
  int i,j;
  // test whether two sketches have the same parameters
  // if so, then they can be added, subtracted, etc. 


  if (!a || !b) return 0;
  if (a->buckets!=b->buckets) return 0;
  if (a->depth!=b->depth) return 0;
  for (i=0;i<a->depth;i++)
    for (j=0;j<6;j++)
      if (a->test[j][i]!=b->test[j][i])
	return 0;
  return 1;
}

int AMS_Count(AMS_type * ams, int item)
{
  // compute the estimated count of item 

  int i;
  int offset;
  int * estimates;
  unsigned int hash;
  int mult;

  estimates=(int *) calloc(1+ams->depth, sizeof(int));
  offset=0;
  for (i=1;i<=ams->depth;i++)
    {
      hash=hash31(ams->test[0][i-1],ams->test[1][i-1],item);
      hash=hash % (ams->buckets); 
      mult=fourwise(ams->test[2][i-1],ams->test[3][i-1],
		    ams->test[4][i-1],ams->test[5][i-1],item);
      if ((mult&1)==1)
	estimates[i]=ams->counts[offset+hash];
      else
	estimates[i]=-ams->counts[offset+hash];
      offset+=ams->buckets;
    }
  if (ams->depth==1) i=estimates[1];
  else if (ams->depth==2) i=(estimates[1]+estimates[2])/2; 
  else
    i=MedSelect(1+ams->depth/2,ams->depth,estimates);
  free(estimates);
  return(i);
}

long long AMS_F2Est(AMS_type * ams)
{
  // estimate the F2 moment of the vector (sum of squares)

  int i,j, r;
  long long * estimates;
  long long result, z;

  estimates=(long long *) calloc(1+ams->depth, sizeof(long long));
  r=0;
  for (i=1;i<=ams->depth;i++)
    {
      z=0;
      for (j=0;j<ams->buckets;j++)
	{
	  z+=((long long) ams->counts[r]* (long long) ams->counts[r]);
	  r++;
	}
      estimates[i]=z;
    }
  if (ams->depth==1) result=estimates[1];
  else if (ams->depth==2) result=(estimates[1]+estimates[2])/2; 
  else
    result=LLMedSelect(1+ams->depth/2,ams->depth,estimates);
  free(estimates);
  return(result);
}

long long AMS_InnerProd(AMS_type * a, AMS_type * b){
  int i,j, r;
  long long * estimates;
  long long result, z;
  // estimate the innerproduct of two vectors using their sketches.

  if (AMS_Compatible(a,b)==0) return 0;
  estimates=(long long *) calloc(1+a->depth, sizeof(long long));
  r=0;
  for (i=1;i<=a->depth;i++)
    {
      z=0;
      for (j=0;j<a->buckets;j++)
	{
	  z+=((long long) a->counts[r]* (long long) b->counts[r]);
	  r++;
	}
      estimates[i]=z;
    }
  if (a->depth==1) result=estimates[1];
  else if (a->depth==2) result=(estimates[1]+estimates[2])/2; 
  else
    result=LLMedSelect(1+a->depth/2,a->depth,estimates);
  free(estimates);
  return(result);

}

int AMS_AddOn(AMS_type * dest, AMS_type * source){
  int i,j,r;

  // add one sketch to another

  r=0;
  if (AMS_Compatible(dest,source)==0) return 0;
  for (i=0;i<source->depth;i++)
    for (j=0;j<source->buckets;j++)
      {
	dest->counts[r]+=source->counts[r];
	r++;
      }
  return 1;
}

int AMS_Subtract(AMS_type * dest, AMS_type * source){
  int i,j,r;

  // subtract one sketch from another

  r=0;
  if (AMS_Compatible(dest,source)==0) return 0;
  for (i=0;i<source->depth;i++)
    for (j=0;j<source->buckets;j++)
      {
	dest->counts[r]-=source->counts[r];
	r++;
      }
  return 1;
}

int AMS_Size(AMS_type * ams){
  int size;
  
  // return the space used in bytes of the sketch

  size=(sizeof(int *))+(ams->buckets*ams->depth)*sizeof(int)+
    ams->depth*6*sizeof(long long)+sizeof(AMS_type);
  return size;
}

void AMS_Destroy(AMS_type * ams)
{
  // destroy the data structure

  int i;

  if (ams)
    {
      for (i=0;i<6;i++)
	free(ams->test[i]);
      free(ams->counts);
      free(ams);
    }
}

/********************************************************************
Combinatorial Group Testing to Find Frequent Items
G. Cormode 2002, 2003

Original Code: 2002-11
This version: 2003-10

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA. 
*********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "cgt.h"
#include "prng.h"

CGT_type * CGT_Init(int buckets, int tests, int lgn, int gran)
{
  // Create the data structure for Combinatorial Group Testing
  // Keep T tests.  Each test has buckets buckets
  // lgn is the bit depth of the items which will arrive
  // this code assumes lgn <= 32 since it manipulates unsigned ints
  // gran is the granularity at which to perform the testing
  // gran = 1 means to do one bit at a time,
  // gran = 4 means to do one nibble at time
  // gran = 8 means to do one octet at a time, etc. 

  int i;
  CGT_type * result;
  prng_type * prng;

  prng=prng_Init(-3254512,2);

  result=calloc(1,sizeof(CGT_type));
  if (result==NULL) exit(1);
  result->tests=tests;
  result->logn=lgn;
  result->gran=gran;
  result->buckets=buckets;
  result->subbuckets=1+(lgn/gran)*((1<<gran) - 1);
  result->count=0;
  result->testa=calloc(tests,sizeof(long long));
  result->testb=calloc(tests,sizeof(long long));
  // create space for the hash functions

  result->counts=calloc(buckets*tests,sizeof(int *));
  if (result->counts==NULL) exit(1); 
  // create space for the counts
  for (i=0;i<buckets*tests;i++)
    {
      result->counts[i]=calloc(result->subbuckets,sizeof(int));
      if (result->counts[i]==NULL) exit(1); 
    }
  for (i=0;i<tests;i++)
    {
      result->testa[i]=(long long) prng_int(prng);
      if (result->testa[i]<0) result->testa[i]= -result->testa[i];
      result->testb[i]=(long long) prng_int(prng);
      if (result->testb[i]<0) result->testb[i]= -result->testb[i];
      // initialise the hash functions
      // prng_int() should return a random integer
      // uniformly distributed in the range 0..2^31
    }
  prng_Destroy(prng);
  return (result);
}

void shell(unsigned long n, int a[])
{
  // A shell sort routine taken from the web
  // to sort the output of the Group Testing

  unsigned long i,j,inc;
  int v;
  inc=1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do { 
    inc /= 3;
    for (i=inc;i<=n;i++) { 
      v=a[i];
      j=i;
      while (a[j-inc] > v) {
	a[j]=a[j-inc];
	j -= inc;
	if (j < inc) break;
      }
      a[j]=v;
    }
  } while (inc > 1);
}


int cmp(const void *p, const void *q)
{
  int a = *((int *) p);  int b = *((int *) q);
  if (a<b) return -1; else if (b>a) return 1; else return 0;
}

unsigned int findone(int *count, int n, int gran, int thresh) 
{
  // find if there is a frequent item in a set of counts

  int i,k,l,offset;
  int countabove, sum, last;

  k=0;
  if (count[0]>=thresh) 
  // if the count is not above threshold, then reject
    {
      offset=1;
      for (i=n;i>0;i-=gran)
	{
	  k<<=gran;
	  countabove=0; sum=0; last=0;
	  for (l=1;l<(1<<gran);l++)
	    {
	      if (count[offset]>=thresh)
		{
		  countabove++;
		  last=l;
		}
	      sum+=count[offset++];
	    }
	  if (count[0]-sum>=thresh)
	    countabove++;
	  if (countabove!=1)
	    {  // check: if both halves of a group are above threshold,
	      // then reject the whole group
	      k=0;
	      break;
	    }	  
	  k+=last;
	  // Update the record of the identity of the frequent item
	}
    }
  return k; // return the identity of the frequent item if there was one
  // this will return zero if there was none. 
}

void loginsert(int *lists, int val, int length, int gran, int diff) 
{
  // add on a value of diff to the counts for item val
  int i;
  int bitmask, offset;

  bitmask=(1<<gran)-1; 
  lists[0]+=diff; // add onto the overall count for the group
  offset=((length/gran)*bitmask)-bitmask;
  for (i=length;i>0;i-=gran) 
    {
      if ((val&bitmask)!=0) // if the lsb = 1, then add on to that group
	lists[offset+(val&bitmask)]+=diff;
      val>>=gran; // look at the next set of bits
      offset-=bitmask;
    }
}

void CGT_Update(CGT_type *cgt, int newitem, int diff)
{
  // receive an update and process the groups accordingly

  int i;
  unsigned int hash;
  int offset=0;

  cgt->count+=diff;
  for (i=0;i<cgt->tests;i++) 
    {
      hash=hash31(cgt->testa[i],cgt->testb[i],newitem);
      hash=hash % (cgt->buckets); 
      loginsert(cgt->counts[offset+hash],newitem,cgt->logn,cgt->gran,diff);
      offset+=cgt->buckets;
    }
}

unsigned int * CGT_Output(CGT_type * cgt, int thresh)
{
  // Find the hot items by doing the group testing

  int i=0,j=0,k=0;
  unsigned int guess=0;
  unsigned int * results, *compresults;
  int hits =0;
  int last=-1;  
  int claimed=0;  
  int testval=0;
  int pass = 0;
  int hash=0;
  
  results=calloc(cgt->tests*cgt->buckets,sizeof(unsigned int));
  if (results==NULL) exit(1); 
  // make some space for the list of results
  
  for (i=0;i<cgt->tests;i++)
    {
      for (j=0; j<cgt->buckets; j++)      
	{      
	  guess=findone(cgt->counts[testval],cgt->logn,cgt->gran,thresh);
	  // go into the group, and see if there is a frequent item there
	  // then check item does hash into that group... 
	  if (guess>0) 
	    {
	      hash=hash31(cgt->testa[i],cgt->testb[i],guess);
	      hash=hash % cgt->buckets; 
	    }
	  if ((guess>0) && (hash==j))
	    {
	      pass=1;
	      for (k=0;k<cgt->tests;k++) 
		{
		  // check every hash of that item is above threshold... 
		  hash=hash31(cgt->testa[k],cgt->testb[k],guess);
		  hash=(cgt->buckets*k) + (hash % (cgt->buckets));
		  if (cgt->counts[hash][0]<thresh)
		    pass=0;
		}
	      if (pass==1)
		{ 
		  // if the item passes all the tests, then output it
		  results[hits]=guess;
		  hits++;
		}
	    }
	  testval++;
	}
    }
  if (hits>0)
    {
      // sort the output
      shell(hits-1,results);
      last=0; claimed=0;
      for (i=0;i<hits;i++)
	{ 
	  if (results[i]!=last)
	    {   // For each distinct item in the output...
	      claimed++;
	      last=results[i];
	    }
	}

      compresults=(unsigned int *) calloc(claimed+1,sizeof(unsigned int));
      compresults[0]=claimed;
      claimed=1; last=0;

      for (i=0;i<hits;i++)
	{ 
	  if (results[i]!=last)
	    {   // For each distinct item in the output...
	      compresults[claimed++]=results[i];
	      last=results[i];
	    }
	}
    }
  else
    {
      compresults=(unsigned int *) malloc(sizeof(unsigned int));
      compresults[0]=0;
    }
  free(results);
  return(compresults);
}  

int CGT_Size(CGT_type *cgt)
{
  int size;
  size=2*cgt->tests*sizeof (long long) + 
    cgt->buckets*cgt->tests*(cgt->subbuckets*sizeof(int))+sizeof(CGT_type);
  return(size);
}

void CGT_Destroy(CGT_type * cgt)
{
  // Free all the space used
  int i;

  free(cgt->testa);
  free(cgt->testb);

  for (i=0; i<cgt->buckets;i++)
    free(cgt->counts[i]);
  free(cgt->counts);
  free (cgt);
}

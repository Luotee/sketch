/********************************************************************
Find large absolute changes between two data streams
G. Cormode 2003,2004

Last modified: 2004-09-14

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA. 
*********************************************************************/

#include <stdlib.h>
#include "prng.h"
#include "massdal.h"
#include "change.h"

/******************************************************************/

#define min(x,y)	((x) < (y) ? (x) : (y))
#define max(x,y)	((x) > (y) ? (x) : (y))

void quitmemory(void * pointer)
     // quit if a memory allocation failed
{
  if (pointer==NULL) exit(1);
}

void loginsert(int *lists, int val, int length, int diff) 
{
  // internal routine used in update
  // lists is a list of 'length' counts
  // val is the item being added
  // diff is the amount (positive or negative)
  //    that its count changes by

  int i;

  // update the logn different tests for a particular item  
  lists[0]+=diff;
  for (i=length;i>0;i--) 
    {
      if ((val&1)==1)
	lists[i]+=diff;
      val>>=1;
    }
}

void floginsert(float *lists, int val, int length, float diff) 
{
  // internal routine used in update
  // lists is a list of 'length' counts
  // val is the item being added
  // diff is the amount (positive or negative)
  //    that its count changes by

  int i;

  // update the logn different tests for a particular item  
  lists[0]+=diff;
  for (i=length;i>0;i--) 
    {
      if ((val&1)==1)
	lists[i]+=diff;
      val>>=1;
    }
}

/******************************************************************/

AbsChange_type * AbsChange_Init(int width, int depth, int lgn)
{
  // initialize the data structure for finding absolute changes
  // width = 1/eps = width of hash functions
  // depth = number of independent repetitions to avoid misses
  // lgn = number of bits in representation of item indexes

  int i;
  prng_type * prng;
  AbsChange_type * absc;
  
  prng=prng_Init(3152131,2);
  // use the random number generator to choose the hash functions

  absc=(AbsChange_type *) calloc(1,sizeof(AbsChange_type));
  absc->depth=depth;
  absc->width=width;
  absc->size=width*depth;
  absc->lgn=lgn;
  // take the input parameters and put them into the struct


  absc->testa=(long long *) calloc(depth,sizeof(long long));
  absc->testb=(long long *) calloc(depth,sizeof(long long));
  absc->counts=(int **) calloc(absc->size,sizeof(int *));
  quitmemory(absc->testa);
  quitmemory(absc->testb);
  quitmemory(absc->counts);
  // make space for the hash functions

  for (i=0;i<absc->size;i++)
    {
      absc->counts[i]=(int *) calloc(1+lgn,sizeof(int));
      quitmemory(absc->counts[i]);
    }
  // make space for the counters

  for (i=0;i<depth;i++)
    {
      absc->testa[i]=(long long) prng_int(prng);
      if (absc->testa[i]<0) absc->testa[i]= -(absc->testa[i]);
      absc->testb[i]=(long long) prng_int(prng);
      if (absc->testb[i]<0) absc->testb[i]= -(absc->testb[i]);
    }
  // create the hash functions
  prng_Destroy(prng);

  return(absc);
}

int AbsChange_Size(AbsChange_type * absc)
{
  // output the size (in bytes) used by the data structure

  int size;

  size=absc->size*sizeof(int)*absc->lgn + absc->depth*2*sizeof(unsigned long)
    + sizeof(AbsChange_type);
  return (size);
}

void AbsChange_Destroy(AbsChange_type * absc)
{
  // free up the space that was allocated for the data structure

  int i;

  for (i=0;i<absc->size;i++)
    free(absc->counts[i]);
  free(absc->counts);
  free(absc->testb);
  free(absc->testa);
  free(absc);
}

unsigned long absfindone(int *count, int l, int thresh) 
{
  // internal routine used to find deltoids. 

  int i,c;
  unsigned long j,k;

  // search through a set of tests to detect whether there is a deltoid there
    // returns 0 if none found, returns id of found item otherwise. 

  j=1;
  k=0;

  if (abs(count[0])<thresh) 
    k=0;
  else
    {
      c = count[0];      
      for (i=l;i>0;i--)
	{
	  // main test: if one side is above threshold and the otherside is not
	  if (((abs(count[i])<thresh) && (abs(c-count[i])<thresh)) ||
	      ((abs(count[i])>=thresh) && (abs(c-count[i])>=thresh)))
	    {
	      k=0;
	      break; 
              // if test fails, bail out
	    }	  
	  if(abs(count[i])>=thresh) 
	    k+=j;
	  j=j<<1;
	  // build the binary representation of the item
	}
    }
  return k;
}

void AbsChange_Update(AbsChange_type * absc, unsigned long newitem, int diff) 
{
  // routine to update the count of an item
  // newitem is the idenfitier of the item being updated
  // diff is the change, positive or negative
  // absc is the data structure

  int i;
  unsigned long hash;

  // for each set of groups, find the group that the item belongs in, update it
  for (i=0;i<absc->depth;i++) 
    {
      hash=hash31(absc->testa[i],absc->testb[i],newitem);
      hash=hash % absc->width; 
      // use the hash function to find the place where the item belongs
      loginsert(absc->counts[i*absc->width+hash],newitem,absc->lgn,diff);
      // call external routine to update the counts
    }
}

unsigned long * AbsChange_Output(AbsChange_type * absc, int thresh)
{
  // take output from the data structure
  // thresh is the threshold for being a deltoid
  // absc is the data structure that holds the changes

  int i,j,k;
  unsigned long guess;
  unsigned long * results;
  int testval=0;
  unsigned long hash=0;
    
  // search for deltoids and put them in a list
  results=(unsigned long *) calloc(absc->size,sizeof(unsigned long));
  results[0]=0;
  // first, create a list
  
  for (i=0;i<absc->depth;i++)
    {
      for (j=0; j<absc->width; j++)      
	{      
	  // go over all the different tests and see if there is a 
	  // deltoid within each test
	  guess=absfindone(absc->counts[testval],absc->lgn,thresh);
	  if (guess>0) 
	    {
	      hash=hash31(absc->testa[i],absc->testb[i],guess);
	      hash=hash % (absc->width); 
	      // check item does hash into that bucket... 

	      if (hash==j)
		{
		  // supress duplicates in output
		    // may be a little slow for large lists
		      for (k=1;k<=results[0];k++)
			if (results[k]==guess) break;
		  if (results[k]!=guess)	
		    results[++results[0]]=guess;
		  // if everything worked out OK, add to the output
		}
	    }
	  // advance to next item
	  testval++;
	}
    }
  return (results);
}  

/******************************************************************/

VarChange_type * VarChange_Init(int width, int depth, int lgn, int streams)
{
  // initialize the data structure for finding variance changes
  // width = 1/eps = width of hash functions
  // depth = number of independent repetitions to avoid misses
  // lgn = number of bits in representation of item indexes
  // streams = number of streams to monitor

  int i, j;
  VarChange_type *  varc;
  prng_type * prng;
  
  varc=(VarChange_type *) malloc(sizeof(VarChange_type));
  // create the space for the structure
  
  prng=prng_Init(732572,2);

  varc->depth=depth;
  varc->width=width;
  varc->size=width*depth;
  varc->lgn=lgn;
  varc->streams=streams;
  // copy the parameters into the struct

  quitmemory(varc->testa=(long long *) calloc(depth,sizeof(long long)));
  quitmemory(varc->testb=(long long *) calloc(depth,sizeof(long long)));
  quitmemory(varc->counts=(int ***) calloc(varc->size,sizeof(int **)));
  // allocate memory for the hash functions

  for (j=0;j<varc->size;j++)
    {
      varc->counts[j]=(int **) calloc(streams,sizeof(int *));
      quitmemory(varc->counts[j]);
      for (i=0;i<streams;i++)
	varc->counts[j][i]= (int *) calloc(1+lgn,sizeof(int));
    }
  // allocate memory for the counters

  for (j=0; j<4;j++)
    {
      varc->wise[j]=(long long *) calloc(depth,sizeof(long long));
      quitmemory(varc->wise[j]);
    }
  // allocate memory for the 4wise independent hash functions

  for (i=0;i<depth;i++)
    {
      for (j=0; j<4; j++)
	{
	  varc->wise[j][i]=(long long) prng_int(prng);
	  if (varc->wise[j][i]<0) varc->wise[j][i]= -varc->wise[j][i];
	}
      varc->testa[i]=(long long) prng_int(prng);
      if (varc->testa[i]<0) varc->testa[i]= -varc->testa[i];
      varc->testb[i]=(long long) prng_int(prng);
      if (varc->testb[i]<0) varc->testb[i]= -varc->testb[i];
    }
  // create the hash functions

  prng_Destroy(prng);
  return(varc); 
}

void VarChange_Destroy(VarChange_type * varc)
{
  // remove the space allocated for the data structure

  int i, j;
  
  for (j=0; j<4;j++)
    free(varc->wise[j]);
  
  for (j=0;j<varc->size;j++)
    {
      for (i=0;i<varc->streams;i++)
	free(varc->counts[j][i]);
      free(varc->counts[j]);
    }
  free(varc->counts);
  free(varc->testb);
  free(varc->testa);
  free(varc);
}

void VarChange_Update(VarChange_type * varc, unsigned long newitem, 
		      int diff, int strm) 
{
  // update the data structure with a new item
  // newitem = the new item identifier
  // diff = the change in its count
  // strm = the stream it occurs in

  int i;
  unsigned long hash;

  for (i=0;i<varc->depth;i++) 
    {
      if ((fourwise(varc->wise[0][i],varc->wise[1][i],
		    varc->wise[2][i],varc->wise[3][i],newitem)&1)  == 1)
	diff=-diff;
      // compute the 4wise independent hash function onto +1/-1
      hash=hash31(varc->testa[i],varc->testb[i],newitem);
      hash=hash % (varc->width); 
      loginsert(varc->counts[varc->width*i + hash][strm],
		newitem,varc->lgn,diff);
      // insert into the count structure
      // can use the same function as in the absolute case
    }
}

int varfindone(int **count, int n, double thresh, int strms) 
{
  // internal routine to find variational changes 

  int i,j,k;
  double c, strs;
  double scores[2][64];
  double avgs[2][64];
  // the size of these arrays is hardcoded to 64
  // to avoid allocating and destroying memory the whole time
  // this limits things to 64 bit item identifiers
  // but my compiler doesn't have larger variables
  // so that shouldn't be problem for now
  
  j=1;
  k=0;
  strs= (double) strms;
  for (i=0;i<=n;i++)
    { avgs[0][i]=0.0; avgs[1][i]=0.0; }
  for (i=0;i<=n;i++) 
    {
      for (j=0;j<strms;j++) 
	{
	  avgs[0][i]+=count[j][i];
	  avgs[1][i]+=(count[j][0]-count[j][i]);
	}
      
      avgs[0][i]=avgs[0][i]/strs;
      avgs[1][i]=avgs[1][i]/strs;
    }

  // the first step is to compute the average value in each dimension

  for (i=0;i<=n;i++)
    {scores[0][i]=0.0; scores[1][i]=0.0;}

  for (j=0;j<strms;j++) 
    {
      scores[0][0]=scores[0][0]+
	(count[j][0]-avgs[0][0])*(count[j][0]-avgs[0][0]);

      for (i=1;i<=n;i++)
	{
	  scores[0][i]+=(count[j][i]-avgs[0][i])*(count[j][i]-avgs[0][i]);

	  scores[1][i]+=(count[j][0]-count[j][i]-avgs[1][i])*
	    (count[j][0]-count[j][i]-avgs[1][i]);
	}
    }
  // next, we compute the sum squared deviation from the mean
  // ie, the variance in each dimension
  // and store these variances in the scores array

  if (scores[0][0]<thresh) 
    {
      k=0;
    }
  else
    {
      c = scores[0][0];      
      j=1;
      for (i=n;i>0;i--)
	{
	  if (((scores[0][i]<thresh) && (scores[1][i]<thresh)) ||
	      ((scores[0][i]>=thresh) && (scores[1][i]>=thresh)))

	    // reject if both sides are above the threshold
	    // or both sides are below threshold
	    {
	      k=0;
	      break;
	    }	  
	  if(scores[0][i]>=thresh) 
	    k+=j;
	  j=j<<1;
	}
    }
  return k;
  // return the item that was found as a deltoid, if any
}

unsigned long * VarChange_Output(VarChange_type * varc, double thresh)
     // output the items with large variance
     // thresh = threshold
{
  int i,j, k, testval;
  unsigned long guess;
  unsigned long * results;
  unsigned long hash=0;
  
  results=(unsigned long *) calloc(varc->size,sizeof(unsigned long));
  // allocate space for the output

  testval=0;
  for (i=0;i<varc->depth;i++)
    {
      for (j=0; j<varc->width; j++)      
	{      
	  // walk over the the data structure

	  guess=varfindone(varc->counts[testval++],varc->lgn,thresh,
			   varc->streams);
	  // use the above procedure to test a set of counters

	  if (guess>0) 
	    {
	      hash=hash31(varc->testa[i],varc->testb[i],guess);
	      hash=hash % varc->width; 
	    }
	  // check item does hash into that bucket...
	  if ((guess>0) && (hash==j))
	    {
	      // supress duplicates in output
	      // may be a little slow for large lists
	      for (k=1;k<=results[0];k++)
		if (results[k]==guess) break;
	      if (results[k]!=guess)	
		results[++results[0]]=guess;
	      // if everything worked out OK, add to the output
	    }
	}
    }
  return (results);
}  

long long VarChange_EstimateVariance(VarChange_type * varc)
{
  // estimate the variance of the data that has been observed so far.

  int i,j,k, pointer;
  long long * results;
  long long mean, sqddev, l;

  results=(long long *) calloc(varc->streams+1,sizeof(long long));
  // allocate space for the intermediate results

  pointer=0;
  for (i=0;i<varc->depth;i++)
    {
      // we get one estimator for each row of the data structure
      sqddev=0;
      for (j=0;j<varc->width;j++)
	{
	  // compute the average of the values
	  mean=0;
	  for (k=0;k<varc->streams;k++)
	    mean+=varc->counts[pointer][k][0];
	  mean=mean/varc->streams;
	  for (k=0;k<varc->streams;k++)
	    {
	      l=varc->counts[pointer][k][0]-mean;
	      sqddev+=l*l; 
	      // compute the deviation from the mean
	      // and add the square of this to the variance
	    }
	  pointer++;
	}
      results[i]=sqddev;
      // note, we don't normalize this by the number of streams
    }
  if (varc->depth==1)
    sqddev=results[1];
  else 
    if (varc->depth==2)
      sqddev=(results[2]+results[1])/2;
    else 
      sqddev= LLMedSelect(1+varc->depth/2,varc->depth,results);
  return sqddev;
  // compute the median of the estimators and return that
}

int VarChange_Size(VarChange_type * varc)
{
  // output the size (in bytes) used by the data structure

  int size;

  size=
    varc->streams*(varc->size*sizeof(int)*varc->lgn + 
		 varc->depth*2*sizeof(unsigned long))
    + sizeof(VarChange_type);
  return (size);
}

/******************************************************************/

RelChange_type *  RelChange_Init(int width, int depth, int lgn)
{
  // initialize the data structure for finding variance changes
  // width = 1/eps = width of hash functions
  // depth = number of independent repetitions to avoid misses
  // lgn = number of bits in representation of item indexes

  int i;
  prng_type * prng;  
  RelChange_type * relc;

  prng=prng_Init(44545,2);

  relc=(RelChange_type *) calloc(1,sizeof(RelChange_type));
  relc->depth=depth;
  relc->width=width;
  relc->size=width*depth; 
  relc->lgn=lgn;
  // take the input parameters and put them into the struct

  relc->testa=(long long *) calloc(depth,sizeof(long long));
  relc->testb=(long long *) calloc(depth,sizeof(long long));
  relc->counts=(float **) calloc(2*relc->size,sizeof(float *));
  quitmemory(relc->testa);
  quitmemory(relc->testb);
  quitmemory(relc->counts);
  // make space for the hash functions

  for (i=0;i<2*relc->size;i++)
    {
      relc->counts[i]=(float *) calloc(1+lgn,sizeof(float));
      quitmemory(relc->counts[i]);
    }
  // make space for the counters

  for (i=0;i<depth;i++)
    {
      relc->testa[i]=(long long) prng_int(prng);
      if (relc->testa[i]<0) relc->testa[i]= -(relc->testa[i]);
      relc->testb[i]=(long long) prng_int(prng);
      if (relc->testb[i]<0) relc->testb[i]= -(relc->testb[i]);
    }
  // create the hash functions
  prng_Destroy(prng);

  return(relc);
}


int testzero(float x)
{
  // crude internal test for something to be too small to bother with
  if (x<=0.00001) return 1; else return 0;
}


unsigned long relfindone(float *count, float *icount, int n, float thresh) 
{
  // internal procedure to find the identity of a relative deltoid

  int i;
  unsigned long j,k;
  float c,e, f, g;
  
  j=1;
  k=0;

  c = count[0];      
  e = icount[0];
  if (testzero(e)) e=1.0;

  if ((c*e)<thresh) 
    k=0;
  else
    {
      for (i=n;i>0;i--)
	{
	  f=icount[i];
	  if (testzero(f)) f=1.0;
	  g=e-f;
	  if (testzero(g)) g=1.0;
	  if (( ((count[i]*f)<thresh) && ((c-count[i])*g<thresh) )  ||
	      ( ((count[i]*f)>=thresh) && ((c-count[i])*g>=thresh)))
	    // if both sides of a test are positive or negative, 
	    // reject the whole set of counters
	    {
	      k=0;
	      break;
	    }	  
	  if((count[i]*f)>=thresh) 
	    {
	      k+=j;
	    }
	  j=j<<1;
	  // build the binary representation of the deltoid
	}
    }
  return k;
  // return the identity of the found item, if any
}

void RelChange_Update(RelChange_type * relc, unsigned long newitem, 
		      float diff, int stream)
{
  // update the relative changes data structure
  // newitem = the item identifier
  // stream = 0 if it is the stream of unadjusted values
  // stream = 1 if it is the stream of inverted values
  // diff = the change for the item (it should be reciprocated
  //        if it is coming from the inverted stream)

  int i;
  unsigned long offset;
  unsigned long hash;

  if (stream==1) offset=relc->size; else offset=0;
  // for each set of groups, find the group that the item belongs in, update it
  for (i=0;i<relc->depth;i++) 
    {
      hash=hash31(relc->testa[i],relc->testb[i],newitem);
      hash=hash % relc->width; 
      // use the hash function to find the place where the item belongs
      floginsert(relc->counts[offset+hash],newitem,relc->lgn,diff);
      // call external routine to update the counts
      offset+=relc->width;
      //printf("Inserting type %d offset %d value %f\n",stream,offset,diff);
    }
}

unsigned long * RelChange_Output(RelChange_type * relc, float thresh)
{
  // output the relative deltoids
  // thresh = the threshold above which to declare a deltoid

  int i,j, k, testval=0;
  unsigned long guess=0,hash=0;
  unsigned long *results;
 
  results=(unsigned long *) calloc(relc->size,sizeof(unsigned long));
  quitmemory(results);

  for (i=0;i<relc->depth;i++)
    {
      for (j=0; j<relc->width; j++)      
	{      
	  guess=relfindone(relc->counts[testval],
			relc->counts[relc->size+testval],relc->lgn,thresh);
	  // for a given test, see if a deltoid can be recovered from it

	  if (guess!=0)
	    {
	      hash=hash31(relc->testa[i],relc->testb[i],(long long) guess);
	      hash=hash % relc->width; 
	      // check that it hashes into the right place

	      if (hash==j) 
		{
		  // supress duplicates in output
		  // may be a little slow for large lists
		  for (k=1;k<=results[0];k++)
		    if (results[k]==guess) break;
		  if (results[k]!=guess)	
		    results[++results[0]]=guess;
		  // if everything worked out OK, add to the output
		}
	    }
	  testval++;
	}
    }
  return (results);
}  

int RelChange_Size(RelChange_type * relc)
{
  // output the size (in bytes) used by the data structure

  int size;

  size=relc->size*sizeof(int)*relc->lgn*2 + relc->depth*2*sizeof(unsigned long)
    + sizeof(RelChange_type);
  return (size);
}

void RelChange_Destroy(RelChange_type * relc)
{
  // free up the space that was allocated for the data structure

  int i;

  for (i=0;i<relc->size;i++)
    free(relc->counts[i]);
  free(relc->counts);
  free(relc->testb);
  free(relc->testa);
  free(relc);
}

/******************************************************************/

/********************************************************************
Implementation of Lossy Counting algorithm to Find Frequent Items
Based on the paper of Manku and Motwani, 2002
Implementation by G. Cormode 2002, 2003

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
#include "lossycount.h"

LC_type * LC_Init(float phi)
{
  LC_type * result;

  result=(LC_type *) calloc(1,sizeof(LC_type));
  result->buckets=0;
  result->holdersize=0;
  result->epoch=0;

  result->window=(int) 1.0/phi;
  result->maxholder=result->window*4;  
  result->bucket=(Counter*) calloc(result->window+2,sizeof(Counter));
  result->holder=(Counter*) calloc(result->maxholder,sizeof(Counter));
  result->newcount=(Counter*) calloc(result->maxholder,sizeof(Counter));
  return(result);
}

void LC_Destroy(LC_type * lc)
{
  free(lc->bucket);
  free(lc->holder);
  free(lc->newcount);
  free(lc);
}

void countershell(unsigned long n, Counter a[])
{
  unsigned long i,j,inc;
  Counter v;
  inc=1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do { 
    inc /= 3;
    for (i=inc+1;i<=n;i++) { 
      v=a[i-1];
      j=i;
      while (a[j-inc-1].item > v.item) {
	a[j-1]=a[j-inc-1];
	j -= inc;
	if (j < inc) break;
      }
      a[j-1]=v;
    }
  } while (inc > 1);
}

int countermerge(Counter *newcount, Counter *left, Counter *right, 
		      int l, int r, int maxholder) 
{  // merge up two lists of counters. returns the size of the lists. 
  int i,j,m;

  if (l+r>maxholder)
    { // a more advanced implementation would do a realloc here...
      printf("Out of memory -- trying to allocate %d counters\n",l+r);
      exit(1);
    }      
  i=0;
  j=0;
  m=0;
  
  while (i<l && j<r)
    { // merge two lists
      if (left[i].item==right[j].item) 
	{ // sum the counts of identical items
	  newcount[m].item=left[i].item;
	  newcount[m].count=right[j].count;
	  while (left[i].item==right[j].item)
	    {
	      newcount[m].count+=left[i].count;
	      i++;
	    }
	  j++;
	}
      else if (left[i].item<right[j].item) 
	{ // else take the left item, creating counts appropriately
	  newcount[m].item=left[i].item;
	  newcount[m].count=0;
	  while (left[i].item==newcount[m].item)
	    {
	      newcount[m].count+=left[i].count;
	      i++;
	    }
	}
      else {
	newcount[m].item=right[j].item;
	newcount[m].count=right[j].count;
	j++;
      }
      newcount[m].count--;
      if (newcount[m].count>0) m++;
      else 
	{ // adjust for items which may have negative or zero counts
	  newcount[m].item=-1;
	  newcount[m].count=0;
	}
    }

  // now that the main part of the merging has been done
  // need to copy over what remains of whichever list is not used up

  if (j<r)
    { 
      while (j<r) 
	{
	  if (right[j].count > 1) 
	    {
	      newcount[m].item=right[j].item;
	      newcount[m].count=right[j].count-1;
	      m++;
	    }
	  j++;
	}
    }
  else 
    if (i<l)
      { 
	while(i<l)
	    {
	      newcount[m].item=left[i].item;
	      newcount[m].count=-1;
	      while ((newcount[m].item==left[i].item) && (i<l))
		{
		  newcount[m].count+=left[i].count;
		  i++;
		}
	      if (newcount[m].count>0) 
		m++;
	      else 
		{ 
		  newcount[m].item=-1;
		  newcount[m].count=0;
		}
	    }
      }
  
  return(m);
}


void LC_Update(LC_type * lc, int val)
{
  Counter *tmp;

  // interpret a negative item identifier as a removal
  if (val>0) 
    { 
      lc->bucket[lc->buckets].item=val;    
      lc->bucket[lc->buckets].count=1;
    }
   else
     {
       lc->bucket[lc->buckets].item=-val;
       lc->bucket[lc->buckets].count=-1;
     }
  
  if (lc->buckets==lc->window) 
    {
      countershell(lc->window,lc->bucket);
      lc->holdersize=countermerge(lc->newcount,lc->bucket,lc->holder,
				  lc->window,lc->holdersize,lc->maxholder);
      tmp=lc->newcount;
      lc->newcount=lc->holder;
      lc->holder=tmp;
      lc->buckets=0;
      lc->epoch++;
    }
  else 
    lc->buckets++;
}

int LC_Size(LC_type * lc)
{
  int size;
  size=(lc->maxholder+lc->window)*sizeof(Counter)+sizeof(LC_type);
  return size;
}

unsigned int *  LC_Output(LC_type * lc, int thresh)
{
  //int correct=0;
  //int claimed=0;
  int i,point;  

  unsigned int * results;

  results=(unsigned int *) calloc(lc->window, sizeof(unsigned int));
  point=1;
  // should do a countermerge here.

  for (i=0;i<lc->holdersize;i++) 
    {
      if (lc->holder[i].count+lc->epoch>=thresh)
	results[point++]=lc->holder[i].item;
    }
  results[0]=point-1;
  return results;
}

/********************************************************************
Find large changes between two data streams
G. Cormode 2003,2004

Last modified: 2004-09-19

This shows off the routines in the library change.c

It creates two streams synthetically from zipf distributions, and then
calls various routines on those streams in order to find the items
with large changes, or 'deltoids' in the language of "What's New:
Finding Significant Changes in Network Data".  The methods tested are
the Absolute Changes, Relative Changes and Variance Changes using
Group Testing, from the above paper.  We also compare to sampling, and
to finding the changes by building sketches.

The data is created as follows: the first stream is drawn from a zipf
distribution with parameter z, then values are hashed randomly into a
20-bit domain.  The second stream is also drawn from the same
distribution, but is offset by one.  The effect of this is that the
count of most items is almost zero in one of the streams.  In general,
such a stream is easier to deal with, since one just needs to find the
heavy hitters in each stream (real data is typically not so easy to
find deltoids in), but it gives a resonable demonstration of how to
use the methods in this library.

Note that we compute the exact threshold for finding the absolute
deltoids.  Based on knowledge of the formation of the streams, we are
able to use this as the basis of the threshold for variational and
relative changes as well.  In general, a different threshold is needed
for these change types.  See the paper for details.

The user can experiment with various parameters from the command line:
z -- the zipf parameter to create the streams.  0 = very uniform, 
     3.0 = very skewed.  Default is 0.8

n -- the length of streams to generate.  default is 50,000 (each)

phi -- the fraction for being a deltoid.  default = 0.005 (0.5%)

width -- the "width" of the data structure to use.  default is 200.
      this is small, values of 1000 or more are typical.

depth -- number of repetitions of the data structure to use. 
      default is 2, this is also small, but 4 or 5 is usually sufficient. 

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <strings.h>

#include "change.h"
#include "prng.h"
#include "massdal.h"
#include "countmin.h"

/******************************************************************/

#define min(x,y)	((x) < (y) ? (x) : (y))
#define max(x,y)	((x) > (y) ? (x) : (y))

float phi, zipfpar; // parameters for random streams 
int width, depth; // parameters of the data structure
unsigned long ** streams; //stream1, *stream2;
int slen;         // length of the stream to use

unsigned long * delts, * reldelts; // a list of deltoids found

typedef unsigned long duo[2];

void CheckArguments(int argc, char **argv) {

  /*****************************************************************/
  /* Examine the command line arguments.  Pick up the file to read */
  /* from the first argument, and also look for other parameters   */
  /*****************************************************************/
  int i;  

  if (argc < 1) {
    printf("%s n z phi width depth\n",argv[0]);
    exit(1);
  }
  printf("Called with: ");
  for (i=0; i<argc; i++)
    printf("%s ",argv[i]);
  printf("\n");

  // read in parameters from the command line, and set to 
  // some reasonable defaults if not present
  if (argc>1)
    slen=atoi(argv[1]);
  else 
    slen=50000;
  // first parameter is the length of streams to consider

  if (argc>2)
    zipfpar=atof(argv[2]);
  else 
    zipfpar=0.8;
  // second parameter gives the skewness of the zipf 
  // distribution to use to generate the streams

  if (argc>3)
    phi=atof(argv[3]);
  else phi=0.005;
  // phi is the fraction above which something is a deltoid
  // (significant change)
  // this is typically small, say 0.01 or 0.001

  if (argc>4)
    width=atoi(argv[4]);
  else width=200;
  // this is the width to use for the data structure
  // it should be proportional to 1/error desired.
  // in turn, error should be less than phi, but
  // in practice this isn't always necessary 
  // width should be a few hundred to a few thousand
  // maybe 10,000 or higher for best accuracy, depending on 
  // memory available

  if (argc>5)
    depth=atoi(argv[5]);
  else 
    depth=2;
  // depth is the number of repetitions of the testing with 
  // different hash functions. 
  // more repetitions drive down the probability of missing 
  // something exponentially. 
  // in practice, quite small values seem to suffice, 
  // but setting this to 5 or more gives very good accuracy
  // although increasing this slows the updates down, since
  // "depth" copies are run in parallel. 
}

unsigned long * CreateStream(int length, int offset, float z)
{
  // generate a stream based of values drawn from a zip distribution

  long a,b;
  float zet;
  int i; 
  unsigned long value, *stream;
  prng_type * prng;

  stream=(unsigned long *) calloc(length+1,sizeof(unsigned long));
      
  prng=prng_Init(44545,2);
  a = (long long) (prng_int(prng)% MOD);
  b = (long long) (prng_int(prng)% MOD);
  zet=zeta(length,z);
  for (i=1;i<=length;i++) 
    {
      value= 
	(offset+hash31(a,b,
		(int) floor(fastzipf(zipfpar,length,zet,prng))))&1048575;
      // get a value from the zipf dbn, and hash it to a new place
      // use offset to mix things up a bit 
      stream[i]=value;
    }
  stream[0]=length;
  prng_Destroy(prng);
  return(stream);
}

int dcmp(const void *x, const void *y)
{
  // used in quicksort to sort lists to get the exact results for comparison
  int *h1,*h2;
  
  h1=x; h2=y;

  if ((*h1)<(*h2)) 
    return -1; 
  else if ((*h1)>(*h2)) 
    return 1;
  else return 0;
}

void adddelt(unsigned long item, long count1, long count2, 
	     unsigned long * list, int thresh, 
	     unsigned long *list2, double dthresh)
{
  // add an item to the list of deltoids if its count is high enough

  if (abs(count1-count2)>=thresh)  
    list[++list[0]]=item;
  // this creates the list of absolute deltoids


  if (count2==0)
    {
      if ((double) count1> dthresh)
	list2[++list2[0]]=item;
    }
  else
    if (((double) count1/ (double) count2) > dthresh)
      list2[++list2[0]]=item;
  // the second part is compute a list of relative deltoids
}


void ShowResults(unsigned long * r) 
{
  // display the contents of a list of items whose length is stored 
  // in the 0th entry of the list 
  // this isn't called in this version of the code, but was for 
  // debugging purposes

  int i;

  for (i=1; i<=r[0]; i++)
    fprintf(stdout, "%ld, ",r[i]);
}

int Collect(unsigned long * str1, duo * c1)
{
  // go through a list of items, and merge together copies of the same
  // item, and give a count of how many times each item occurs.

  int i,prevptr, t, collect;

  qsort(&str1[1], str1[0], sizeof(unsigned long), dcmp);
  // sort the input list
  prevptr=0; t=0;
  collect=0;
  for (i=1;i<=str1[0];i++)
    {
      if (str1[i]!=str1[prevptr])
	{
	  c1[t][0]=str1[prevptr];
	  c1[t][1]=collect;
	  t++;
	  prevptr=i;
	  collect=0;
	  // record the number of times the previous item was seen, and 
	  // set up for the next one.
	}
      collect++;
    }
  c1[t][0]=str1[prevptr];
  c1[t][1]=collect;
  return(t);
}

int reformatstream(unsigned long *str1, unsigned long *str2, 
		    unsigned long * results, unsigned long * relresults)
{
  // a critical procedure in creating the exact solutions 
  // in order to evaluate the approximate solutions. 

  int i, j, t1, t2, thresh;
  duo *c1;
  duo *c2;
  long long netpos;
  double netchange, rc, relthresh;

  // allocate space for the sorted output  
  c1=(duo *) calloc(str1[0]+2,2*sizeof(unsigned long));
  c2=(duo *) calloc(str2[0]+2,2*sizeof(unsigned long));
  t1=Collect(str1,c1);
  t2=Collect(str2,c2);
  // first we sort the streams and compact duplicates

  netpos=0;i=1;j=1; netchange=0.0;
  // first pass to compute maximum difference, second to test threshhold
  while ((i<=t1) && (j<=t2))
    {
      if (c1[i][0]==c2[j][0])
	{
	  rc=((double) c1[i][1]) / ((double) c2[j][1]);
	  netchange+=rc;
	  netpos+=abs(c1[i++][1] - c2[j++][1]);
	}
      else 
	if (c1[i][0]<c2[j][0])
	  {
	    netchange+=(double) c1[i][1]; // normalize missing value to 1
	    netpos+=abs(c1[i++][1]);
	  }
	else if (c1[i][0]>c2[j][0])
	  // does not count towards netchange
	  netpos+=abs(c2[j++][1]);
    }
  while (j<=t2) 
    netpos+=abs(c2[j++][1]);
  while (i<=t1)
    {
      netchange+=(double) c1[i][1];
      netpos+=abs(c1[i++][1]);
    }
  // at this point, netpos records the total L1 difference 
  // between the two streams
  
  thresh=((float) netpos) * phi; 
  if (thresh==0) thresh=1;
  //derive the threshold for being a deltoid from the L1 difference

  relthresh= netchange * phi;
  // if desired, can also work out threshold for relative deltoids

  i=0; j=0;
  // a second pass over the two streams lets us compute the difference
  // in count for each item, and test whether this is greater than the
  // threshold for being a deltoid
  while ((i<=t1) && (j<=t2))
    {
      if (c1[i][0]==c2[j][0])
	{
	  adddelt(c1[i][0], c1[i][1] , c2[j][1], results, thresh, 
		  relresults, relthresh);
	  i++; j++;
	}
      else 
	if (c1[i][0]<c2[j][0])
	  {
	    adddelt(c1[i][0],c1[i][1],0, results, thresh, 
		    relresults, relthresh);
	    i++;
	  }
	else if (c1[i][0]>c2[j][0])
	  {
	    adddelt(c2[j][0],0,c2[j][1], results, thresh, 
		    relresults, relthresh);
	    j++;
	  }
    }
  while (j<=t2) 
    {
      adddelt(c2[j][0],0,c2[j][1], results, thresh, 
	      relresults, relthresh);
      j++;
    }
  while (i<=t1)
    {
      adddelt(c1[i][0],c1[i][1],0, results, thresh, 
	      relresults, relthresh);
      i++;
    }
  return thresh;
}

void ScoreResults(unsigned long * reslts, unsigned long * delts)
{
  // go through list, check if it coincides with the deltoids. 
  // this assumes the list is indexed from 0 and no repeats in list
  int i,j,hits;

  hits=0;
  for (i=1;i<=reslts[0];i++)
    for (j=1;j<=delts[0];j++)
      if (reslts[i]==delts[j]) {
	hits++;
	break;
      }
  printf("  Recall: %.2f Precision: %.2f\n",(float) hits/delts[0], 
	 (float) hits/reslts[0]);

  // the recall is the fraction of deltoids that were identified
  // the precision is the fraction of the guesses that were correct.
}


void sample(unsigned long **streams, float p)
{
  // a simple routine to simulate sampling a small fraction of 
  // the streams, and using these samples to search for deltoids
  // in general, sampling can be very bad at this task

  unsigned long *samp, *rsamp;
  unsigned long *ss[2];
  int i, strm;
  float r, stime;
  prng_type * prng;

  StartTheClock();
  samp=(unsigned long *) calloc((int) (p*2.0*streams[0][0]) + 2,
				sizeof(unsigned long));
  rsamp=(unsigned long *) calloc((int) (p*2.0*streams[0][0]) + 2,
				sizeof(unsigned long));
  // create space for the sample

  prng=prng_Init(34212,2);
  for (strm=0;strm<=1;strm++)
    {
      ss[strm]=calloc((int)(p*2.0*streams[strm][0])+2,sizeof(unsigned long));
      for (i=1; i<=streams[strm][0];i++)
	{
	  r=prng_float(prng); 
	  if (r<p) ss[strm][++ss[strm][0]]=streams[strm][i];
	  // add an item to the sample if the randomness says so
	}
    }
  stime=StopTheClock();
  reformatstream(ss[0],ss[1],samp, rsamp);
  // now compute the exact answer from the samples 
  // using above routines for the complete exact answer 

  fprintf(stdout,"Sampling rate is %f, sampled %ld items (%ld bytes) in %.1fms\n",
	  p, ss[0][0]+ss[1][0],(ss[0][0]+ss[1][0])*sizeof(long),stime);
  ScoreResults(samp,delts);

  /*  ScoreResults(samp,reldelts); */
  // can also find the quality of sampling for finding relative deltoids */

  free(ss[0]);
  free(ss[1]);
  free(samp);
}

void sketchbased(unsigned long ** streams, int w, int d, int thresh)
{
  // an approach to finding changes based on using sketches
  // in general, this approach doesn't work, since it needs
  // two passes: one to build sketches and find threshold
  // a second to query the streams.


  CM_type * cm; 
  unsigned long * results;
  int i,j,estdif;
  float sktime;

  results=(unsigned long *) calloc(streams[0][0],sizeof(unsigned long));
  // allocate space for the output
  // this is much bigger than is needed, but is done for simplicity

  cm=CM_Init(w,d,5722119);
  StartTheClock();
  for (i=1; i<=streams[0][0];i++)
    CM_Update(cm,streams[0][i],+1);
  // build a sketch of stream 1 width w depth d
  for (i=1; i<=streams[1][0];i++)
    CM_Update(cm,streams[1][i],-1);
  // build a sketch of stream 2 using the same parameters

  sktime=StopTheClock();
  for (i=1;i<=streams[1][0];i++)
    {
      estdif= abs(CM_PointMed(cm,streams[1][i]));
      // replay str2, and query difference of sketches 
      // to find the deltoids (use correct thresh as in group testing)
      if (estdif>=thresh) // add str2[i] to output
	{
	  for (j=1;j<=results[0];j++)
	    if (results[j]==streams[1][i]) break;
	  if (results[j]!=streams[1][i])
	    {
	      results[++results[0]]=streams[1][i];
	    }
	}
    }
  fprintf(stdout,"Sketch Based processed %ld items in %.1fms, using %d bytes\n",
	  streams[0][0]+streams[1][0],sktime,CM_Size(cm));
  ScoreResults(results,delts);
  free(results);
  CM_Destroy(cm);
}

int AbsWrap(unsigned long ** streams, int lgn, int thresh)
{
  // this routine calls the change detection methods based on group
  // testing described in Cormode-Muthukrishnan (Infocom 2004)
  // this implements the tests for large absolute changes

  AbsChange_type * absc;
  unsigned long * results;
  int sized, i;
  float gttime;

  absc=AbsChange_Init(width,depth,lgn);
  // initialize the data structure
  // with the width and depth to use, plus the number of bits 
  // in the items being input (eg 20 or 32)

  StartTheClock();
  for (i=1; i<=streams[0][0]; i++)
    AbsChange_Update(absc,streams[0][i],1);
  for (i=1; i<=streams[1][0]; i++)
    AbsChange_Update(absc,streams[1][i],-1);
  // update the data structure with the stream values
  // +1 for the first stream, -1 for the second stream

  gttime=StopTheClock();
  results=AbsChange_Output(absc,thresh);
  // get the output from the data structure

  sized=AbsChange_Size(absc);
  // compute the amount of space used 

  fprintf(stdout,"Group Testing processed %ld items in %.1fms, using %d bytes\n",
	 streams[0][0]+streams[1][0],gttime,sized);

  ScoreResults(results,delts);
  // compare the answers found to the exact answers

  free (results);
  AbsChange_Destroy(absc);
  // free up the space used

  return (sized);
}

void VarWrap(unsigned long ** streams, int lgn, int thresh, int nostreams)
{
  // this routine calls the change detection methods based on group
  // testing described in Cormode-Muthukrishnan (Infocom 2004)
  // this implements the tests for large variational changes

  VarChange_type * varc;
  unsigned long * results;
  int i, j, len;
  float gttime;

  varc=VarChange_Init(width,depth,lgn,2);
  // initialize the data structure for 2 streams

  StartTheClock();
  len=0;
  for (i=0;i<nostreams;i++)
    {
      len+=streams[i][0];
      for (j=1;j<=streams[i][0];j++)
	VarChange_Update(varc,streams[i][j],1,i);
    }
  gttime=StopTheClock();

  // The following is specific to the case of two streams
  // If there are only two streams, then the variance simplifies
  // to just sum (a[i] - b[i])^2/2
  // so we can use the same threshold from before, but square it and
  // halve the result and run the variance algorithm to find absolute
  // changes 
  results=VarChange_Output(varc,0.5*thresh*thresh);
  //  sized=AbsChange_Size(absc);
  fprintf(stdout, "Variance Test processed %d items in %.1fms, using %d bytes\n",
	 len,gttime,VarChange_Size(varc));
  ScoreResults(results,delts);
  free (results);
  VarChange_Destroy(varc);
}


void RelWrap(unsigned long ** streams, int lgn, int thresh)
{
  // this routine calls the change detection methods based on group
  // testing described in Cormode-Muthukrishnan (Infocom 2004)
  // this implements the tests for large relative changes

  RelChange_type * relc1;
  unsigned long * results;
  duo * inv1;
  int sized, i,t1;
  float gttime;

  inv1=(duo *) calloc(streams[1][0]+2,2*sizeof(unsigned long));
  t1=Collect(streams[1],inv1);

  relc1=RelChange_Init(10*width,2*depth,lgn);

  // initialize the data structure
  // with the width and depth to use, plus the number of bits 
  // in the items being input (eg 20 or 32)
  // relative change requires much more space to find deltoids

  StartTheClock();
  for (i=1; i<=streams[0][0]; i++)
    RelChange_Update(relc1,streams[0][i],1,0);

  for (i=1; i<=t1; i++)
    RelChange_Update(relc1,inv1[i][0],1.0/((float)inv1[i][1]),1);
  // update the data structure with the stream values

  gttime=StopTheClock();
  results=RelChange_Output(relc1,(double) thresh/2.0);
  // get the output from the data structure

  sized=RelChange_Size(relc1);
  // compute the amount of space used 

  fprintf(stdout,"Rel Change processed %ld items in %.1fms, using %d bytes\n",
	 streams[0][0]+streams[1][0],gttime,sized);

  ScoreResults(results,reldelts);
  // compare the answers found to the exact answers

  free (results);
  RelChange_Destroy(relc1);
  // free up the space used
}


/******************************************************************/

int main(int argc, char **argv) 
{
  int lgn, n, nostreams, i; 
  float extime, rate;  
  int asize;        // space used by the approximate methods
  int thresh;       // threshold for being a deltoid

  CheckArguments(argc,argv);
  fprintf(stderr,
	  "____________________________________________________________\n");
  fprintf(stderr,
	  "%s compiled at %s, %s\n", __FILE__, __TIME__, __DATE__);

  n=1048576; // this defines the universe size
  nostreams=2; // fix two streams

  streams=(unsigned long **) calloc(nostreams,sizeof(unsigned long *));
  StartTheClock();
  for (i=0;i<nostreams;i++)
    streams[i]=CreateStream(slen,52521*i,zipfpar);
  // create the streams based on the input parameters

  delts=(unsigned long *) calloc((int) 1/phi,sizeof(unsigned long));  
  reldelts=(unsigned long *) calloc((int) 1/phi,sizeof(unsigned long));  
  lgn=ceil(log((float) n)/log(2.0));
  thresh=reformatstream(streams[0],streams[1],delts,reldelts);
  extime=StopTheClock();
  fprintf(stdout, "Created streams in %.1f ms \n",extime);
  // compute the exact changes from the synthetic streams


  asize=AbsWrap(streams,lgn,thresh);
  // call the main group testing routines

  VarWrap(streams,lgn,thresh,nostreams);
  // call the group testing for variational changes

  RelWrap(streams,lgn,thresh);
  // call the group testing for relative changes

  rate=(float) asize/((slen*nostreams)*sizeof(unsigned long));
  sample(streams,rate);
  // compare against a sampling based solution

  sketchbased(streams,(lgn/4)*width,4*depth,thresh);
  // and compare against using sketches
  // with appropriate scaling of parameters to give equal space

  /* Done! */
  return 0;
}


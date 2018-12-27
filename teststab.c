/********************************************************************
Approximate frequent items in a data stream
G. Cormode 2002

Last modified: 2004-01

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prng.h"
#include "massdal.h"

/******************************************************************/

#include "stable.h"
#include "ams.h"
#include "ccfc.h"
#include "fm.h"

/******************************************************************/

float zipfpar, phi;
int range, distinct;
int *exact;
int d;
int n=1000; 
int netpos;  
int quartiles[4],max;
long long sumsq;

/******************************************************************/

int * CreateStream(int length)
{
  long a,b;
  float zet;
  int i; 
  long value;
  int * stream;
  prng_type * prng;

  n=1048575;
  exact=(int *) calloc(n+1,sizeof(int));
  stream=(int *) calloc(length+1,sizeof(int));
      
  prng=prng_Init(44545,2);
  a = (long long) (prng_int(prng)% MOD);
  b = (long long) (prng_int(prng)% MOD);

  netpos=0;

  zet=zeta(length,zipfpar);

  for (i=1;i<=length;i++) 
    {
      value= 
	(hash31(a,b,((int) floor(fastzipf(zipfpar,n,zet,prng)) ))&1048575);
      exact[value]++;
      netpos++;
      stream[i]=value;
      // printf("Stream %d is %d \n",i,value);
    }

  prng_Destroy(prng);

  return(stream);

}

/******************************************************************/

void RunExact()
{
  int i;

  sumsq=0; distinct=0;
  for (i=0;i<n;i++) 
    {
      sumsq+=   (long long) exact[i]*  (long long) exact[i];
      if (exact[i]>0) 
	{
	  distinct++;
	}
    }
}


/******************************************************************/

int main(int argc, char **argv) 
{
  int * stream, i;
  double z;

  AMS_type * ams;
  CCFC_type * ccfc;
  Stable_sk * sk1, * sk2, * sk0, *sk002;
  FM_type * fm;

  printf("____________________________________________________________\n");
  printf("%s compiled at %s, %s\n", __FILE__, __TIME__, __DATE__);
  
  if (argc>1)
    range=atoi(argv[1]);
  else
    range=12345;

  if (argc>2)
    zipfpar=atof(argv[2]);
  else zipfpar=0.8;

  if ((range<=0) || (zipfpar<0.0))
    {
      printf("Usage: %s range zipfpar\n",argv[0]);
      printf("range = number of values to generate\n");
      printf("zipfpar = parameter of zipfdistribution\n");
      exit(1);
    }

  stream=CreateStream(range);

  printf("\nTesting L2 Norm\n\n");

  RunExact();

  sk2=Stable_Init(128,2.0,54211);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      Stable_Update(sk2,stream[i],1.0);
    else
      Stable_Update(sk2,stream[i],-1.0);

  ams=AMS_Init(512,5);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      AMS_Update(ams,stream[i],1);      
    else
      AMS_Update(ams,-stream[i],-1);      

  ccfc=CCFC_Init(512,5,1,1);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      CCFC_Update(ccfc,stream[i],1);      
    else
      CCFC_Update(ccfc,-stream[i],-1);      

  printf("Correct value of L2 is   %lld \n",sumsq);
  printf("AMS  estimate of L2 is   %lld \n",AMS_F2Est(ams));
  printf("CCFC estimate of L2 is   %lld \n",CCFC_F2Est(ccfc));
  z=Stable_norm(sk2);
  printf("Stable estimate of L2 is %lf \n",z*z);

  printf("\nTesting L1 Norm\n\n");

  sk1=Stable_Init(128,1.0,13461);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      Stable_Update(sk1,stream[i],1.0);
    else
      Stable_Update(sk1,stream[i],-1.0);
  printf("Correct value of L1 is   %d \n",range);
  z=Stable_norm(sk1);
  printf("Stable estimate of L1 is %lf \n",z);

  printf("\nTesting L0 Norm\n\n");

  printf("Correct value of L0 is        %d \n",distinct);

  fm=FM_Init(128,112351);
  for (i=1;i<=range;i++)
    if (stream[i]>0) // current implementation of FM is inserts only
      FM_Update(fm,stream[i]);
    else
      FM_Update(fm,-stream[i]);
  z=FM_Distinct(fm);
  printf("FlajoletMartin estimate L0 is %lf \n",z);

  sk0=Stable_Init(128,0.0,13461);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      Stable_Update(sk0,stream[i],1.0);
    else
      Stable_Update(sk0,stream[i],-1.0);
  z=Stable_norm(sk0)/1.43;
  printf("Fast stable estimate of L0 is %lf \n",z);
  sk002=Stable_Init(128,0.02,13461);
  for (i=1;i<=range;i++) 
    if (stream[i]>0)
      Stable_Update(sk002,stream[i],1.0);
    else
      Stable_Update(sk002,stream[i],-1.0);
  z=Stable_norm(sk002)/1.43;
  printf("Slow stable estimate of L0 is %lf \n",z);

  AMS_Destroy(ams);
  Stable_Destroy(sk2);
  Stable_Destroy(sk1);
  FM_Destroy(fm);  

  printf("\n");
  /* Done! */
  return 0;
}


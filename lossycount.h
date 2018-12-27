// lossycount.h -- header file for Lossy Counting
// see Manku & Motwani, VLDB 2002 for details
// implementation by Graham Cormode, 2002,2003

typedef struct counter
{
  int item;
  int count;
} Counter;

typedef struct LC_type
{
  Counter *bucket;
  Counter *holder;
  Counter *newcount;
  int buckets;
  int holdersize;
  int maxholder;
  int window;
  int epoch;
} LC_type;

extern LC_type * LC_Init(float);
extern void LC_Destroy(LC_type *);
extern void LC_Update(LC_type *, int);
extern int LC_Size(LC_type *);
extern unsigned int * LC_Output(LC_type *,int);

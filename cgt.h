// cgt.h -- header file for Combinatorial Group Testing, Graham Cormode
// 2002,2003

typedef struct CGT_type{
  int tests;
  int logn;
  int gran;
  int buckets;
  int subbuckets;
  int count;
  int ** counts;
  int *testa, *testb;
} CGT_type;

extern CGT_type * CGT_Init(int, int, int, int);
extern void CGT_Update(CGT_type *, int, int); 
extern unsigned int * CGT_Output(CGT_type *, int);
extern void CGT_Destroy(CGT_type *);
extern int CGT_Size(CGT_type *);

// ccfc.h -- header file for Adaptive Group Testing, Graham Cormode, 2003
// using Count Sketches, proposed in CCFC 2002.

typedef struct CCFC_type{
  int tests;
  int logn;
  int gran;
  int buckets;
  int count;
  int ** counts;
  int *testa, *testb, *testc, *testd;
} CCFC_type;

extern CCFC_type * CCFC_Init(int, int, int, int);
extern void CCFC_Update(CCFC_type *, int, int); 
extern int CCFC_Count(CCFC_type *, int, int);
extern unsigned int * CCFC_Output(CCFC_type *, int);
extern long long CCFC_F2Est(CCFC_type *);
extern void CCFC_Destroy(CCFC_type *);
extern int CCFC_Size(CCFC_type *);


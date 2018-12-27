// ams.h -- header file for Alon-Matias-Szegedy sketches
// using pairwise hash functions to speed up updates, Graham Cormode 2003

typedef struct AMS_type{
  int depth;
  int buckets;
  int count;
  int * counts;
  int *test[6];
} AMS_type;

extern AMS_type * AMS_Init(int, int);
extern void AMS_Update(AMS_type *, unsigned long, int); 
extern long long AMS_F2Est(AMS_type *);
extern long long AMS_InnerProd(AMS_type *, AMS_type *); 
extern int AMS_Subtract(AMS_type *, AMS_type *);
extern int AMS_AddOn(AMS_type *, AMS_type *);
extern void AMS_Destroy(AMS_type *);
extern int AMS_Size(AMS_type *);


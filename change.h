typedef struct AbsChange_type{
  int depth;
  int width;
  int lgn; 
  int size;
  long long *testa, *testb;
  int ** counts;
} AbsChange_type;

extern AbsChange_type * AbsChange_Init(int, int, int);
extern void AbsChange_Update(AbsChange_type *, unsigned long, int); 
extern unsigned long * AbsChange_Output(AbsChange_type *, int); 
extern void AbsChange_Destroy(AbsChange_type *);
extern int AbsChange_Size(AbsChange_type *);

typedef struct VarChange_type{
  int depth;
  int width;
  int lgn; 
  int size;
  long long *testa, *testb;
  int *** counts;
  long long *wise[4];
  int streams;
} VarChange_type;

extern VarChange_type * VarChange_Init(int, int, int, int);
extern void VarChange_Update(VarChange_type *, unsigned long,int,int); 
extern unsigned long * VarChange_Output(VarChange_type *, double);
extern void VarChange_Destroy(VarChange_type *);
extern long long VarChange_EstimateVariance(VarChange_type *);
extern int VarChange_Size(VarChange_type *);


typedef struct RelChange_type{
  int depth;
  int width;
  int lgn; 
  int size;
  long long *testa, *testb;
  float ** counts;
} RelChange_type;

extern RelChange_type * RelChange_Init(int, int, int);
extern void RelChange_Update(RelChange_type *, unsigned long, float, int); 
extern unsigned long * RelChange_Output(RelChange_type *, float); 
extern void RelChange_Destroy(RelChange_type *);
extern int RelChange_Size(RelChange_type *);



typedef struct FM_type {
  int fmsize; // length of the sketch
  unsigned int * fm;
  unsigned int * hasha;
  unsigned int * hashb; 
} FM_type;


extern FM_type * FM_Init(int, int);
extern void FM_Update(FM_type *, unsigned int);
extern double FM_Distinct(FM_type * fm);
extern void FM_Destroy(FM_type *);

 

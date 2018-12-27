//frequent.h -- simple frequent items routine
// see Misra&Gries 1982, Demaine et al 2002, Karp et al 2003
// implemented by Graham Cormode, 2002,2003

typedef struct itemlist ITEMLIST;
typedef struct group GROUP;

struct group 
{
  int diff;
  ITEMLIST *items;
  GROUP *previousg, *nextg;
};

struct itemlist 
{
  int item;
  GROUP *parentg;
  ITEMLIST *previousi, *nexti;
  ITEMLIST *nexting, *previousing;
  
};

typedef struct freq_type{

  ITEMLIST **hashtable;
  GROUP *groups;
  int k;
  int tblsz;
  long long a,b;
} freq_type;


extern freq_type * Freq_Init(float);
extern void Freq_Destroy(freq_type *);
extern void Freq_Update(freq_type *, int);
extern int Freq_Size(freq_type *);
extern unsigned int * Freq_Output(freq_type *,int);

#include "prng.h"

typedef struct Stable_sk {
  double alpha; // The norm we are working in 
  int sksize; // length of the sketch
  double *sk; // the sketch
  long seed;  // seed for prng
  prng_type * prng; // the prng
} Stable_sk;

//#define PI 3.141592653589793
// It is very important to know the ratio between the radius of a circle 
// and its circumference

extern Stable_sk * Stable_Init(int, double, long);
extern void Stable_Update(Stable_sk *, int , double);
extern double Stable_norm(Stable_sk *);
extern void Stable_Destroy(Stable_sk *);

 

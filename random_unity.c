
#include "random_unity.h"

/* Function to generate d numbers at random that add up to unity */

void random_unity(int d, double *x)
{
  int k;
  double sum=0.0;

  for (k=0; k<d; k++) {
    x[k] = gsl_rng_uniform(rng);
    sum += x[k];
  }
  for (k=0; k<d; k++) x[k] /= sum;
}



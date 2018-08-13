/* Program to perform K-group EM/BP community detection using the
 * degree-corrected SBM on an arbitrary network read from a GML file, with
 * discrete (categorical) metadata stored in the "label" field
 *
 * Written by Mark Newman  28 NOV 2014
 */

/* Program control */

#define VERBOSE        // Set to print progress updates to stderr

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "readgml.h"

/* Constants */

#define K 2            // Number of groups

#define BP_ACC 1e-4    // Required accuracy for BP to terminate
#define EM_ACC 1e-4    // Required accuracy for EM to terminate

#define NOCONVERGE     // Set to abort if the solution doesn't converge
#define BP_MAXSTEP 20  // Maximum number of BP steps before aborting
#define EM_MAXSTEP 100 // Maximum number of EM steps before aborting

#define SMALL 1.0e-100

/* Globals */

NETWORK G;             // Struct storing the network
int twom;              // Twice the number of edges

int *x;                // Metadata
int *nx;               // Number of nodes with each distinct metadata value
double **nrx;          // Expected number in each group with value
char **mlabel;         // Metadata strings
int nmlabels;          // Number of distinct metadata strings

double **gmma;         // Prior parameters (spelled "gmma" because "gamma"
                       //   is a reserved word in C math.h)
double omega[K][K];    // Mixing parameters

double ***eta;         // Messages
double **q;            // One-point marginals

gsl_rng *rng;          // Random number generator


/* Get metadata from the labels */

#define MAXMETA 300

void get_metadata()
{
  int u,i;

  /* Make space for the metadata numbers and labels */

  x = malloc(G.nvertices*sizeof(int));
  mlabel = malloc(MAXMETA*sizeof(char*));
  nmlabels = 0;

  /* Go through the vertices */

  for (u=0; u<G.nvertices; u++) {

    /* Check to see if this label is already in the list of labels */

    for (i=0; i<nmlabels; i++) {
      if (strcmp(G.vertex[u].label,mlabel[i])==0) break;
    }

    /* If not, add it */

    if (i==nmlabels) {
      mlabel[nmlabels++] = G.vertex[u].label;  // Just set pointers equal
    }

    /* Record this as the metadata type for this vertex */

    x[u] = i;
  }

  /* Count how many nodes there are in each metadata group */

  nx = calloc(nmlabels,sizeof(int));
  for (u=0; u<G.nvertices; u++) nx[x[u]]++;

#ifdef VERBOSE
  fprintf(stderr,"Found %i distinct metadata values:\n",nmlabels);
  for (i=0; i<nmlabels; i++) fprintf(stderr," %i %s\n",i,mlabel[i]);
#endif
}


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


/* Do BP */

int bp()
{
  int i,j;
  int u,v;
  int r,s;
  int steps;
  double deltaeta,maxdelta;
  double logeta,neweta,sum,norm,largest;
  double d[K];
  double logpre[K];
  double logqun[K];
  double ***logetaun;

  // Make space for the new log-etas, which are called "logetaun" because
  // the are initially calculated in unnormalized form

  logetaun = malloc(G.nvertices*sizeof(double**));
  for (u=0; u<G.nvertices; u++) {
    logetaun[u] = malloc(G.vertex[u].degree*sizeof(double*));
    for (i=0; i<G.vertex[u].degree; i++) {
      logetaun[u][i] = malloc(K*sizeof(double));
    }
  }

  // Main BP loop

  steps = 0;
  do {

    /* Calculate the expected group degrees */
    
    for (r=0; r<K; r++) {
      d[r] = 0.0;
      for (u=0; u<G.nvertices; u++) d[r] += q[u][r]*G.vertex[u].degree;
    }
  
    /* Calculate the log-prefactors (without the leading factor of d_i or
     * the prior) */

    for (r=0; r<K; r++) {
      logpre[r] = 0.0;
      for (s=0; s<K; s++) logpre[r] -= omega[r][s]*d[s];
    }

    /* Calculate new values for the one-vertex marginals */

#ifdef VERBOSE
    fprintf(stderr,"Calculating one-vertex marginals...    \r");
#endif
    for (u=0; u<G.nvertices; u++) {
      for (r=0; r<K; r++) {
	logqun[r] = log(gmma[r][x[u]]) + G.vertex[u].degree*logpre[r];
	for (i=0; i<G.vertex[u].degree; i++) {
	  sum = 0.0;
	  for (s=0; s<K; s++) sum += eta[u][i][s]*omega[r][s];
	  if (sum<SMALL) sum = SMALL;
	  logqun[r] += log(sum);
	}
	if (r==0) largest = logqun[r];
	else if (logqun[r]>largest) largest = logqun[r];
      }

      /* Normalize */

      norm = 0.0;
      for (r=0; r<K; r++) {
	logqun[r] -= largest;
	norm += exp(logqun[r]);
      }
      for (r=0; r<K; r++) q[u][r] = exp(logqun[r])/norm;
    }

    /* Calculate (unnormalized) new values for the (log) messages */

#ifdef VERBOSE
    fprintf(stderr,"Calculating messages...              \r");
#endif
    for (u=0; u<G.nvertices; u++) {
      for (i=0; i<G.vertex[u].degree; i++) {
	v = G.vertex[u].edge[i].target;
	for (r=0; r<K; r++) {
	  logeta = log(gmma[r][x[v]]) + G.vertex[v].degree*logpre[r];
	  for (j=0; j<G.vertex[v].degree; j++) {
	    if (G.vertex[v].edge[j].target!=u) {
	      sum = 0.0;
	      for (s=0; s<K; s++) sum += eta[v][j][s]*omega[r][s];
	      if (sum<SMALL) sum = SMALL;   // Prevent -Inf
	      logeta += log(sum);
	    }
	  }
	  logetaun[u][i][r] = logeta;
	}
      }
    }

    /* Normalize and calculate largest change */

#ifdef VERBOSE
    fprintf(stderr,"Normalizing messages...\r");
#endif
    maxdelta = 0.0;
    for (u=0; u<G.nvertices; u++) {
      for (i=0; i<G.vertex[u].degree; i++) {
	norm = 0.0;
	largest = logetaun[u][i][0];
	for (r=1; r<K; r++) {
	  if (logetaun[u][i][r]>largest) largest = logetaun[u][i][r];
	}
	for (r=0; r<K; r++) {
	  logetaun[u][i][r] -= largest;
	  norm += exp(logetaun[u][i][r]);
	}	  
	for (r=0; r<K; r++) {
	  neweta = exp(logetaun[u][i][r])/norm;
	  deltaeta = fabs(neweta-eta[u][i][r]);
	  if (deltaeta>maxdelta) maxdelta = deltaeta;
	  eta[u][i][r] = neweta;
	}
      }
    }

#ifdef VERBOSE
    fprintf(stderr,"BP steps %i, max change = %g                   \r",
	    steps,maxdelta);
#endif

  } while ((maxdelta>BP_ACC)&&(++steps<=BP_MAXSTEP));

#ifdef VERBOSE
  fprintf(stderr,"\n");
#endif

  // Free space

  for (u=0; u<G.nvertices; u++) {
   for (i=0; i<G.vertex[u].degree; i++) free(logetaun[u][i]);
    free(logetaun[u]);
  }
  free(logetaun);

  return steps;
}


// Function to calculate new values of the parameters

double params()
{
  int u,v;
  int i,j;
  int r,s;
  double norm,esum,quvrs;
  double d[K];
  double term[K][K];
  double sum[K][K];
  double L;

  // Calculate some basics

  for (r=0; r<K; r++) {
    d[r] = 0.0;
    for (i=0; i<nmlabels; i++) nrx[r][i] = 0.0;
  }

  for (u=0; u<G.nvertices; u++) {
    for (r=0; r<K; r++) {
      d[r] += q[u][r]*G.vertex[u].degree;
      nrx[r][x[u]] += q[u][r];
    }
  }

  // Calculate new values of the gammas

  for (r=0; r<K; r++) {
    for (i=0; i<nmlabels; i++) gmma[r][i] = nrx[r][i]/nx[i];
  }

  // Calculate the new values of the omegas

  // Zero out the sum variables

  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) sum[r][s] = 0.0;
  }
  esum = 0.0;

  // Perform the sums

  for (u=0; u<G.nvertices; u++) {
    for (i=0; i<G.vertex[u].degree; i++) {
      v = G.vertex[u].edge[i].target;

      // Find which edge leads back from v to u

      for (j=0; j<G.vertex[v].degree; j++) {
	if (G.vertex[v].edge[j].target==u) break;
      }
      if (j==G.vertex[v].degree) {
	fprintf(stderr,"Error!\n");
	exit(23);
      }

      // Calculate the terms and the normalization factor

      norm = 0.0;
      for (r=0; r<K; r++) {
	for (s=0; s<K; s++) {
	  term[r][s] = omega[r][s]*eta[u][i][r]*eta[v][j][s];
	  norm += term[r][s];
	}
      }

      // Add to the running sums

      for (r=0; r<K; r++) {
	for (s=0; s<K; s++) {
	  quvrs = term[r][s]/norm;
	  sum[r][s] += quvrs;
	  esum += quvrs*log(quvrs);
	}
      }
    }
  }

  // Calculate the new values of the omega variables (after calculating
  // the likelihood using the old omegas)

  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) omega[r][s] = sum[r][s]/(d[r]*d[s]);
  }

  // Calculate the expected log-likelihood

  // Internal energy first

  L = 0.0;
  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) L += 0.5*sum[r][s]*log(omega[r][s]);
    for (i=0; i<nmlabels; i++) {
      if (gmma[r][i]>0.0) L += nx[i]*gmma[r][i]*log(gmma[r][i]);
    }
  }

  // Now the entropy

  L -= 0.5*esum;
  for (u=0; u<G.nvertices; u++) {
    for (r=0; r<K; r++) {
      if ((q[u][r]>0.0)&&(G.vertex[u].degree>0)) {
	L += (G.vertex[u].degree-1)*q[u][r]*log(q[u][r]);
      }
    }
  }

  return L;
}


main(int argc, char *argv[])
{
  int u,v,i,r,s;
  int step;
  int bpsteps;
  double norm,deltac,maxdelta;
  double c[K][K];
  double oldc[K][K];
  double ru[K];
  double L;

  // Initialize random number generator

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,time(NULL));

  // Read the network and the metadata from stdin

#ifdef VERBOSE
  fprintf(stderr,"Reading network...\n");
#endif
  read_network(&G,stdin);
  for (u=twom=0; u<G.nvertices; u++) twom += G.vertex[u].degree;
  get_metadata();

  // Make space for the marginals and initialize to random initial values

  eta = malloc(G.nvertices*sizeof(double**));
  q = malloc(G.nvertices*sizeof(double*));

  for (u=0; u<G.nvertices; u++) {
    q[u] = malloc(K*sizeof(double));
    random_unity(K,q[u]);
  }

  // Make space for the messages and initialize to the same values as the
  // marginals

  for (u=0; u<G.nvertices; u++) {
    eta[u] = malloc(G.vertex[u].degree*sizeof(double*));
    for (i=0; i<G.vertex[u].degree; i++) {
      eta[u][i] = malloc(K*sizeof(double));
      v = G.vertex[u].edge[i].target;
      for (r=0; r<K; r++) eta[u][i][r] = q[v][r];
    }
  }

  nrx = malloc(K*sizeof(double*));
  for (r=0; r<K; r++) nrx[r] = malloc(nmlabels*sizeof(double));

  // Malloc space for parameters gmma and choose random initial values

  gmma = malloc(K*sizeof(double*));
  for (r=0; r<K; r++) gmma[r] = malloc(nmlabels*sizeof(double));
  for (i=0; i<nmlabels; i++) {
    random_unity(K,ru);
    for (r=0; r<K; r++) gmma[r][i] = ru[r];
  }

  // Choose random values for the omegas, but with a bias toward
  // assortative choices (change if necessary for other networks)

  for (r=0; r<K; r++) {
    for (s=0; s<K; s++) {
      if (r==s) c[r][s] = 1 + gsl_rng_uniform(rng);
      else if (r<s) c[r][s] = gsl_rng_uniform(rng);
      else c[r][s] = c[s][r];
      omega[r][s] = c[r][s]/twom;
    }
  }

  // EM loop

#ifdef VERBOSE
  fprintf(stderr,"Starting EM algorithm...\n");
#endif
  step = 0;
  do {

    // Run BP to calculate the messages and one-vertex marginals

    bpsteps = bp();

    // Calculate the new values of the parameters

    L = params();

    // Calculate the new values of the c variables

    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) {
	oldc[r][s] = c[r][s];
	c[r][s] = omega[r][s]*twom;
      }
    }

    // Find the largest change in any of the c's

    maxdelta = 0.0;
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) {
	deltac = fabs(c[r][s]-oldc[r][s]);
        if (deltac>maxdelta) maxdelta = deltac;
      }
    }

    // Print out new values of the parameters

#ifdef VERBOSE
    fprintf(stderr,"EM step %i, max change = %g\n",step,maxdelta);
    fprintf(stderr,"gamma =\n");
    for (r=0; r<K; r++) {
      for (i=0; i<nmlabels; i++) fprintf(stderr," %.6f",gmma[r][i]);
      fprintf(stderr,"\n");
    }

    fprintf(stderr,"c =\n");
    for (r=0; r<K; r++) {
      for (s=0; s<K; s++) fprintf(stderr," %.6f",c[r][s]);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
#endif

    if (++step>EM_MAXSTEP) {
#ifdef NOCONVERGE
      fprintf(stderr,"Solution failed to converge in %i EM steps\n",
	      EM_MAXSTEP);
      break;
#endif
    }

  } while (maxdelta>EM_ACC);

#ifdef NOCONVERGE
  if (bpsteps>BP_MAXSTEP) {
    fprintf(stderr,"BP failed converge on final EM step\n");
  }
#endif

#ifdef VERBOSE
  fprintf(stderr,"Log-likelihood = %g\n\n",L);
#endif

  // Output the results

  for (u=0; u<G.nvertices; u++) {
    printf("%i %s",u,mlabel[x[u]]);
    for (r=0; r<K; r++) printf(" %.6f",q[u][r]);
    printf("\n");
  }
}

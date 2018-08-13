
#include "main.h"

int main(int argc, char *argv[])
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

  return 1;
}

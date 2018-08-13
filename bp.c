#include "bp.h"

/* Do BP */

int bp(void)
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


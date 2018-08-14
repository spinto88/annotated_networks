#include "params.h"

// Function to calculate new values of the parameters

double params(NETWORK G, int k_comm, int *x, int *nx, int nmlabels, double **nrx, double **gmma, double **omega, double ***eta, double **q)
{
  int u,v;
  int i,j;
  int r,s;
  double norm, esum, quvrs;
  double d[k_comm];
  double term[k_comm][k_comm];
  double sum[k_comm][k_comm];
  double L;

  // Calculate some basics

  for (r=0; r<k_comm; r++)
  {
    d[r] = 0.0;
    for (i=0; i<nmlabels; i++) 
      nrx[r][i] = 0.0;
  }

  for (u=0; u<G.nvertices; u++)
  {
    for (r=0; r<k_comm; r++) 
    {
      d[r] += q[u][r]*G.vertex[u].degree;
      nrx[r][x[u]] += q[u][r];
    }
  }

  // Calculate new values of the gammas

  for (r=0; r<k_comm; r++) 
  {
    for (i=0; i<nmlabels; i++) 
      gmma[r][i] = nrx[r][i]/nx[i];
  }

  // Calculate the new values of the omegas

  // Zero out the sum variables

  for (r=0; r<k_comm; r++) 
  {
    for (s=0; s<k_comm; s++) 
      sum[r][s] = 0.0;
  }
  esum = 0.0;

  // Perform the sums

  for (u=0; u<G.nvertices; u++) 
  {

    for (i=0; i<G.vertex[u].degree; i++) 
    {
      v = G.vertex[u].edge[i].target;
      // Find which edge leads back from v to u

      for (j=0; j<G.vertex[v].degree; j++) 
      {
	if (G.vertex[v].edge[j].target==u) 
          break;
      }
      if (j==G.vertex[v].degree) 
      {
	fprintf(stderr,"Error!\n");
	exit(23);
      }

      // Calculate the terms and the normalization factor

      norm = 0.0;
      for (r=0; r<k_comm; r++) 
      {
	for (s=0; s<k_comm; s++) 
        {
	  term[r][s] = omega[r][s]*eta[u][i][r]*eta[v][j][s];
	  norm += term[r][s];
	}
      }

      // Add to the running sums

      for (r=0; r<k_comm; r++) 
      {
	for (s=0; s<k_comm; s++) 
        {
	  quvrs = term[r][s]/norm;
	  sum[r][s] += quvrs;
	  esum += quvrs*log(quvrs);
	}
      }
    }
  }

  // Calculate the new values of the omega variables (after calculating
  // the likelihood using the old omegas)

  for (r=0; r<k_comm; r++) 
  {
    for (s=0; s<k_comm; s++) 
      omega[r][s] = sum[r][s]/(d[r]*d[s]);
  }

  // Calculate the expected log-likelihood

  // Internal energy first

  L = 0.0;
  for (r=0; r<k_comm; r++) 
  {
    for (s=0; s<k_comm; s++) 
      L += 0.5*sum[r][s]*log(omega[r][s]);
    for (i=0; i<nmlabels; i++) 
    {
      if (gmma[r][i]>0.0) 
        L += nx[i]*gmma[r][i]*log(gmma[r][i]);
    }
  }

  // Now the entropy

  L -= 0.5*esum;
  for (u=0; u<G.nvertices; u++) 
  {
    for (r=0; r<k_comm; r++) 
    {
      if ((q[u][r]>0.0)&&(G.vertex[u].degree>0)) 
      {
	L += (G.vertex[u].degree-1)*q[u][r]*log(q[u][r]);
      }
    }
  }

  return L;
}

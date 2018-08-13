#include "get_metadata.h"

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

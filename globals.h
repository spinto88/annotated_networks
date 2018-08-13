
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
#include "network.h"

/* Constants */

#define K 3            // Number of groups

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


/* Program to perform K-group EM/BP community detection using the
 * degree-corrected SBM on an arbitrary network read from a GML file, with
 * discrete (categorical) metadata stored in the "label" field
 *
 * Written by Mark Newman  28 NOV 2014
 */

/* Modified to do the connection with Python by Seba Pinto in Aug 2018 */

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Constants */

#define BP_ACC 1e-4    // Required accuracy for BP to terminate
#define EM_ACC 1e-4    // Required accuracy for EM to terminate

#define BP_MAXSTEP 20  // Maximum number of BP steps before aborting
#define EM_MAXSTEP 100 // Maximum number of EM steps before aborting

#define SMALL 1.0e-100


/* Network definition */

// Header file for VERTEX, EDGE, and NETWORK data structures 
//
// Mark Newman  11 AUG 06

#ifndef _NETWORK_H
#define _NETWORK_H

typedef struct {
  int target;        // Index in the vertex[] array of neighboring vertex.
                     // (Note that this is not necessarily equal to the GML
                     // ID of the neighbor if IDs are nonconsecutive or do
                     // not start at zero.)
  double weight;     // Weight of edge.  1 if no weight is specified.
} EDGE;

typedef struct {
  int id;            // GML ID number of vertex
  int degree;        // Degree of vertex (out-degree for directed nets)
  char *label;       // GML label of vertex.  NULL if no label specified
  EDGE *edge;        // Array of EDGE structs, one for each neighbor
} VERTEX;

typedef struct {
  int nvertices;     // Number of vertices in network
  int directed;      // 1 = directed network, 0 = undirected
  VERTEX *vertex;    // Array of VERTEX structs, one for each vertex
} NETWORK;

#endif

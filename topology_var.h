/*
 * This file defines the common constants and extern variables 
 * for all topologies
 */

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _SMALL_TEST
#define MAX_NODE 513           // Max Number of Nodes in the network
#define MAX_DEGREE 149             // Max Degree Supported..Atiq 6/15/16
#define MAX_NPROCS 128          // Max Number of Processing nodes in the n/w
#else 
#define MAX_NODE 65000
#define MAX_DEGREE 64             // Max Degree Supported
#define MAX_NPROCS 30000          // Max Number of Processing nodes in the n/w
#endif

#define TRUE 1                    // Define TRUE
#define FALSE 0                   // Define FALSE

#define MAX_LONG_PATH_LEN    64

#define MAX_ALLPATH 20000         // largest number of paths between each SD 
                                  // pair

#ifdef _TORUS                     // Define the maximum path length for the 
#define MAX_PATH_LEN 64           // topology
#else                             // The length is 64 for torus
#define MAX_PATH_LEN 15           // for others, may change 
#endif                            // depends on the topology

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif

#ifndef _REAL_PATH_NUM
#define _REAL_PATH_NUM 1
#endif

#define min(a, b) (((a) < (b)) ? (a):(b))

// Define the enumerator for the topologies
typedef enum Topology{
  XGFT, 
  DRAGONFLY, 
  HYPERCUBE, 
  TORUS, 
  JELLYFISH, 
  JELLYFLY,
  TIANHE,
  EDISON,
  NPJELLYFISH,
  KAUTZ,
  SLIMFLY
} Topology;

// Define the enumerator for the routing schemes
typedef enum Routing{
  SINGLEPATH,
  MULTIPATH
} Routing;

typedef struct HopType{
   int node;
   int index;
} Hop;

typedef int Path[MAX_PATH_LEN];
typedef long long int BWMBPS;

static const int numPath = _REAL_PATH_NUM; // This is used for the 
                                           // number of paths

// Define the generic function pointers

void (* routing_algorithm)(int src, int dst, Path path);
void (* model_routing_algorithm)(int src, int dst, int *len, int *rsrc, int *rdst);
void (* model_multipathrouting_algorithm)(int src, int dst, Path kpath[]);

#include "xgft.h"

/* these are the extern variables that should be declared in either the
   driver or simulation engine file
*/

int totNode; /**< Total number of nodes, including processing and switch elements. */
int totSE; /**< Total number of switching elements. */
int totPE; /**< Total number of processing elements */

int nprocs;                       // Actual Number of Processors
Topology topology;                // The topology used

/**
 * Indicates the type of routing that is being used.
 * This is defined by the user at the command line.
 * 
 * @see topology.h
 * 
 */
int routing;

Routing routingType;

/** The per hop latency is defined by the user at commandline. */
int per_hop_latency;
/** The software overhead is defined by the user at commandline. */
int software_overhead;
/** The host bandwidth is different from the switch bandwidth. 
 It is defined by the user at commandline. */

/** 
 * This stores the connection matrix for the topology.
 *
 * @see MAX_NODE
 * @see MAX_DEGREE
 */
int graph[MAX_NODE][2*MAX_DEGREE];
int links[MAX_NODE][MAX_DEGREE];

int graph_m[MAX_NODE][MAX_DEGREE]; 
int graph_m1[MAX_NODE][MAX_DEGREE]; 

/**
 * Stores the bandwidth usage for each SD pair.
 *  
 * @see MAX_NODE
 * @see MAX_DEGREE
 */
long long int bandwidth[MAX_NODE][MAX_DEGREE];
long long int llvar[MAX_NODE][MAX_DEGREE];

/**
 * Stores the load on each edge.
 *  
 * @see MAX_NODE
 * @see MAX_DEGREE
 */
#ifdef _VLB_ROUTING
double load_graph[MAX_NODE][MAX_DEGREE];
#else
int load_graph[MAX_NODE][MAX_DEGREE];
#endif

int per_hop_latency;
int software_overhead;
long long int host_bw;

#endif

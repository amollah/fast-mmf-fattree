#ifndef XGFT_H
#define XGFT_H
#include "topology.h"

#define MAX_H 4

// Routing schemes for the XGFT
// Use block 1 through 20

#define XGFT_DMODK_ROUTING     1
#define XGFT_RANDOM_ROUTING    2
#define XGFT_ADAPTIVE_ROUTING  3
#define XGFT_RR_ROUTING        4
#define XGFT_RR1_ROUTING       5
#define XGFT_RR2_ROUTING       6
#define XGFT_KPATH_ROUTING     7
#define XGFT_ECMP_ROUTING      8
#define XGFT_LADAPTIVE_ROUTING 9
#define XGFT_VLB_ROUTING       10

#define XGFT_MAX_NUM_PATH      1800       // (m/2)*(m/2)....
#define XGFT_MAX_DISJOINT_PATH 50    // W_1

/* xgft related routines */

void (* xgft_routing_algorithm)(int src, int dst, int *path);
void (* xgft_model_routing_algorithm)(int src, int dst, int *len, int *rsrc, int *rdst);

void xgft_build_topology(int h, int *M, int *W, long long int *BW);

void xgft_print_topology();

void xgft_dmodk_routing(int src, int dst, int *path);
void xgft_random_routing(int src, int dst, int *path);
void xgft_adaptive_routing(int src, int dst, int *path);
void xgft_roundrobin_routing(int src, int dst, int *path);
void xgft_roundrobin1_routing(int src, int dst, int *path);
void xgft_roundrobin2_routing(int src, int dst, int *path);
void xgft_kpath_routing(int src, int dst, Path kpath[numPath]);
void xgft_ecmp_routing(int src, int dst, Path kpath[numPath]);
int xgft_allpath_routing(int src, int dst, 
			  int allpath[MAX_ALLPATH][MAX_PATH_LEN]);

void xgft_model_dmodk_routing(int src, int dst, int *len, int *rsrc, int *rdst);
int  xgft_check_path(int *path);
void xgft_print_path(int *path);
void xgft_topology_init(int h, int *m, int *w, long long int *bw, int r);

/*xgft topology helper functions*/
void print_label(int label[MAX_H+1]);
void compute_label(int node, int label[MAX_H+1]);
void copy_label(int *dst, int *src);
void compute_nodeid(int * nodeid, int label[MAX_H+1]);
void check_if_node(int *label, int id);
int get_graph_index(int dst_switch, int sw);
int xgft_find_nca_height(int s, int d);
int xgft_find_npath(int s, int d, int *nca);
int xgft_find_path(int s, int d, int index, Hop *pathHop);
/* external variables */

extern int xgft_h;
extern int xgft_M[MAX_H];
extern int xgft_W[MAX_H];
extern long long int xgft_BW[MAX_H];  /* bandwidth in each layer */

#endif

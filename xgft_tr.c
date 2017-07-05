#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "topology.h"

int xgft_h;
int xgft_M[MAX_H];
int xgft_W[MAX_H];
long long int xgft_BW[MAX_H];  /* bandwidth in each layer */
                               /* in the unit of bps */
Hop pathHop[MAX_PATH_LEN];

int baseL[MAX_H+1];
int sizeL[MAX_H+1];

/* some local function */
static void compute_baseL(int h, int *M, int *W);


/**
 * simple clone function to copy XGFT node labels
 * @param dst destination label array
 * @param src source label array
 */
void copy_label(int *dst, int *src)
{
  int i;
  int h = xgft_h;
  for (i=0; i<=h; i++) {
    dst[i] = src[i];
  }
}

/** compute_baseL computes the base index for each layer of nodes
 * on-time routine to populate baseL(global) array
 * baseL[i] is the starting ID of the first node at level i
 * starting ID increments with i, since in XGFT, IDs are assigned in bottom-up
 * // L[0] are processing nodes while L[i] are switches
*/
static void compute_baseL(int h, int *M, int *W){
  int i, j;
  int tmp;

  for (i=0; i<=h; i++) {
    tmp = 1;
    for(j=0; j<i; j++) {
      tmp *= W[j];
    }
    for (j=i; j<h; j++)
      tmp *=M[j];
    sizeL[i] = tmp;
  }

  tmp = 0;
  for (i=0; i<=h; i++) {
    baseL[i] = tmp;
    tmp += sizeL[i];
  }
  totPE = sizeL[0];
  totNode = baseL[h] + sizeL[h];
  nprocs = totPE;
}

/**
 * translates a given node ID into the corresponding label
 * each lavel is of the format<h,M_h,M_h-1,...,M1>
 */
void compute_label(int nodeid, int label[MAX_H+1])
{
  int i;
  int localid;  
  int K[MAX_H];
  int h;

  h = xgft_h;

  if ((nodeid < 0) || (nodeid > baseL[h] + sizeL[h])) {
    printf("Nodeid %d out of bound.\n", nodeid);
    exit(0);
  }

  for (i=0; i<=h; i++) {
    if (baseL[i] > nodeid) break;
  }
  label[0] = i-1;  /* determine the level */
  
  localid = nodeid - baseL[label[0]];

  for (i=0; i<h; i++) {
    if (i<label[0]) K[i] = xgft_W[i];
    else K[i] = xgft_M[i]; 
  }

  /* compute the encoding of localid */
  for (i=0; i<h; i++) {
    label[i+1] = localid % K[i];
    localid = localid / K[i];
  }
}

void compute_nodeid(int * nodeid, int label[MAX_H+1])
{
  int i, j;
  int localid;  
  int K[MAX_H];
  int h;

  h = xgft_h;

  if ((label[0] < 0)  || (label[0] > h)) {
    printf("level %d does not exist.\n", label[0]);
    exit(0);
  }

  for (i=0; i<h; i++) {
    if (i<label[0]) K[i] = xgft_W[i];
    else K[i] = xgft_M[i]; 
  }

  for (i=0; i<h; i++) {
    if (label[i+1] > K[i]) {
      printf("label[%d] = %d out of bound %d\n", i+1, label[i+1], K[i]);
      for (j=0; j<=h; j++) {
	printf("label[%d] = %d\n", j, label[j]);
      }
      for (j=0; j<h; j++) {
	printf("K[%d] = %d\n", j, K[j]);
      }
      exit(0);
    }
  }

  localid = label[h];

  for (i=h-1; i>=1; i--) {
    localid = localid * K[i-1] + label[i];
  }

  *nodeid = localid + baseL[label[0]];
}

/**
 * Given a source and destination for an XGFT this function will compute
 * the height that the nearest common ancestor will be found for the two nodes.
 *
 * @param s - the source node
 * @param d - the destination node.
 * @return The height of the NCA.
 */
int xgft_find_nca_height(int s, int d){
   int src_label[MAX_H+1];
   int dst_label[MAX_H+1];
   int nca = MAX_H;
   extern int xgft_h;
   
   compute_label(s, src_label);
   compute_label(d, dst_label);
   
   nca = xgft_h;
   while(src_label[nca] == dst_label[nca]) --nca;
   
   return nca;
   
}

/**
 * Calculate the number of paths available between an Source-destination pair.
 *
 * @param s - the source node
 * @param d - the destination node.
 * @param nca - the nearest common ancestor will be stored in here if not null
 * @return The number of paths between the SD-pair.
 */
int xgft_find_npath(int s, int d, int *nca){
	int nca_tmp =0;
	int npath = 1;
	int i = 0;

	nca_tmp = xgft_find_nca_height(s, d);
	for(i = 0; i < nca_tmp; ++i){
		npath *= xgft_W[i];
	}
	if(nca != NULL)
		*nca = nca_tmp;
	return npath;
}

/**
 * Check it the label with id is a node.
 *
 * @param label - the label to be checked.
 * @param id - the ID of the node.
 */
void check_if_node(int *label, int id) {
   if (label[0] != 0) {
      printf("ERROR check_if_node: Node %d, ", id);
      print_label(label);
      printf(" is not a processing node.\n");
      exit(0);
   }
}

/**
 * This method find the index of src_switch for the switch sw.
 * If the source switch can not be found the program will be terminated.
 *
 * @param dst_switch - the switch that we are trying to connect to.
 * @param sw - the switch we are search from.
 * @return src_switch index in graph
 *
 */
int get_graph_index(int dst_switch, int sw) {
   int i = 0;
   while (graph[sw][i] != dst_switch && graph[sw][i] != -1) {
      ++i;
   }

   if (graph[sw][i] == dst_switch)
      return i;
   else {
      printf("ERROR get_graph_index: no edge found (%d, %d)\n", sw,
             dst_switch);
      exit(0);
   }
}
#if 0
/**
 * This function prints out the contents of the pathHop array to see the path
 * that was calculated. By default this function is left out of compilation as
 * it is purely a debug function. Remote the preprocessor #if to use this function.
 */
static void dump_path_hop(){
   int node, index, h;
   for(h = 0; h < MAX_PATH_LEN && pathHop[h].node != -1; ++h){
      node = pathHop[h].node;
      index = pathHop[h].index;
      printf(": %d (%d -> %d) :", node, index, graph[node][index]);
   }
   printf("\n");
}
#endif

static int add_hop(int node, int index, int hop, Hop pathHop[MAX_PATH_LEN]){
   pathHop[hop].node = node;
   pathHop[hop].index = index;
   return (hop+1);
}

/**
 * This function find the path between a source destination pair. Of all shortest paths available the index value indicates which 
 * of those paths should be calculated. The path is stored in the pathHop array.
 *
 * @param s - the source node ID in the network.
 * @param d - the destination node ID in the network.
 * @param index - the index of the path you want to create.
 * @param pathHop - the array to store the calculated path in.
 */
int xgft_find_path(int s, int d, int index, Hop *pathHop){
   int nca = xgft_find_nca_height(s, d);
   int src_label[MAX_H + 1];
   int dst_label[MAX_H + 1];
   int r_label[MAX_H +1];
   int t_label[MAX_H +1];
   int degree = 0;
   int curr_node = 0, next_node=0;
   int hop = 0;
   int l, w, k, j;

   compute_label(s, src_label);
   compute_label(d, dst_label);
   check_if_node(src_label, s);
   check_if_node(dst_label, d);


   copy_label(r_label, src_label);
   r_label[0] = 1;
   r_label[1] = 0;
   curr_node = s;
   
   compute_nodeid(&next_node, r_label);
   degree = get_graph_index(next_node, curr_node);
 
   hop = add_hop(curr_node, degree, hop, pathHop);
   copy_label(t_label, r_label);
   for(l = 1; l < nca; ++l){
      w = 1;
      /*for(k = 0; k < (nca - l); ++k)*/
      for(k = nca-1; k > l; --k)
         w *= xgft_W[k];
      
      j = ( index / w ) % xgft_W[l];
      t_label[ l + 1 ] = j;
      ++t_label[0];
      curr_node = next_node; 
      compute_nodeid(&next_node, t_label);
      degree = get_graph_index(next_node, curr_node);
      hop = add_hop(curr_node, degree, hop, pathHop);
      if(hop >= MAX_PATH_LEN){
         printf("Path length greater than can be stored: %d, max %d\n", hop, MAX_PATH_LEN);
         exit(0);
      }
   }

   for(l = nca - 1; l >=0; --l){
      t_label[ l + 1 ] = dst_label[ l + 1 ];
      --t_label[0];
      curr_node = next_node;
      compute_nodeid(&next_node, t_label);
      degree = get_graph_index(next_node, curr_node);
      hop = add_hop(curr_node, degree, hop, pathHop);
      if(hop >= MAX_PATH_LEN){
         printf("Path length greater than can be stored: %d, max %d\n", hop, MAX_PATH_LEN);
         exit(0);
      }
   }
   hop = add_hop(-1, -1, hop, pathHop);
   //  printf("Map : %d -> %d:\n", s, d);
   //  dump_path_hop();
   return hop;
}


/*
// build_xgft builds the extended generalized fat-tree
*/
void xgft_build_topology(int h, int *M, int *W, long long int *BW)
{
  int i, j, k;
  int nodeid;
  int label[MAX_H+1];
  int code[MAX_H+1];


  for (i=0; i<totNode; i++) 
    for (j=0; j<MAX_DEGREE; j++) {
      graph[i][j] = -1;
    }

  for (i=0; i<totNode; i++) {
    /* compute connection for each node */
    compute_label(i, label);
    /* up first then down */
    j = 0;
    if (label[0] != h) { /* can go up */
      for (; j<W[label[0]]; j++) {
        copy_label(code, label);
	/*        print_label(label); printf("xxxxx\n"); */
        code[code[0]+1] = j;
	code[0] ++;
	/*        print_label(code); printf("yyyy\n"); */
	compute_nodeid(&nodeid, code);
	graph[i][j] = nodeid;
        bandwidth[i][j] = BW[label[0]];
      }
    }

    k=0;
    if (label[0] != 0) {
      for (; k<M[label[0]-1]; k++) {
        copy_label(code, label);
        code[0] --;
        code[code[0]+1] = k;
	/*        print_label(code); printf("zzzzz\n"); */
	compute_nodeid(&nodeid, code);
	graph[i][j+k] = nodeid;
        bandwidth[i][j+k] = BW[code[0]];
      }
    }
  }
}

void print_label(int label[MAX_H+1])
{
  int i;
  int h;
  
  h = xgft_h;
  printf("(");
  for (i=0; i<h; i++) printf("%d, ", label[i]);
  printf("%d)", label[i]);
}

void xgft_print_topology()
{
  int i, j;
  int label[MAX_H+1];

  for (i=0; i<totNode; i++) {
    compute_label(i, label);
    printf("Node %d, ", i);
    print_label(label);
    printf(", connects to\n");
    for (j=0; graph[i][j] != -1; j++) {
      printf("  node %d (bw=%lld), ", graph[i][j], bandwidth[i][j]);
      compute_label(graph[i][j], label);
      print_label(label);
      printf("\n");
    }
  }
}

void xgft_print_path(int *path)
{
  int i = 0;
  int label[MAX_H+1];
  while (path[i] != -1) {
    compute_label(path[i], label);
    printf("node %d, ", path[i]);
    print_label(label);
    printf("\n");
    i++;
  }
}

int xgft_allpath_routing(int src, int dst, 
                                 int allpath[MAX_ALLPATH][MAX_PATH_LEN]){
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[MAX_H+1];
  int t_label[MAX_H+1];
  int index[MAX_H];
  int i,j,k,l,w = 1,npath = 1;
  int high;  
  int tdst, ttdst;
  int h;  

  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    allpath[0][0] = src;
    allpath[0][1] = -1;
    return 0;
  }

  high = h;
  tdst = dst;
  while (src_label[high] == dst_label[high]) high --;


  for(i=0; i<high; i++) npath *= xgft_W[i];

  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } 
  else {
    //printf("There are %d path between %d %d high = %d\n",npath,src,dst,high);
    for(i=0;i<MAX_ALLPATH;i++)for(j=0;j<MAX_PATH_LEN;j++)allpath[i][j]=-1;
    for(i=0; i<MAX_H; i++)index[i] = 0; index[0] = index[0];

    for(i=0; i< npath; i++){
      allpath[i][0] = src;
      copy_label(r_label, src_label);
      r_label[0] = 1;
      r_label[1] = 0;
      compute_nodeid(&allpath[i][1], r_label);
      copy_label(t_label, r_label);
      ttdst = tdst; ttdst = ttdst;
      for (l=1; l < high; l++) {
        //w = 1;for(k=0; k < (h - l); k++) w*= xgft_W[k];
	//        w = 1;for(k=0; k < (high - l); k++) w*= xgft_W[k];
        w = 1;for(k=high-1; k>l; k--) w*= xgft_W[k];
        j = (i/w) % xgft_W[l];
        t_label[l+1] = j;
        t_label[0] ++;
        compute_nodeid(&allpath[i][l+1], t_label);
      }
      
      /* down path */
      for (l = high -1; l>=0; l--) {
        t_label[l+1] = dst_label[l+1];
        t_label[0] --;
        j = 2*high - l;
        compute_nodeid(&allpath[i][j], t_label);
      }
      allpath[i][2*high+1] = -1;
    }
  }
  return npath;
  //xgft_print_topology();
  //for(i=0; i<npath; i++) xgft_print_path(allpath[i]);
}

void xgft_alldisjointpath_routing(int src, int dst, 
                             int djpath[XGFT_MAX_DISJOINT_PATH][MAX_PATH_LEN]){
  int i,j,k;
  int found = 1, npath =1 , high, m = 1, n;
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int allpath[XGFT_MAX_NUM_PATH][MAX_PATH_LEN];
  int path[MAX_PATH_LEN];
  //int djpath[XGFT_MAX_DISJOINT_PATH][MAX_PATH_LEN];
  
  //h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);
  
  for(i=0; i<MAX_PATH_LEN; i++)path[i] = -1;
  
  xgft_allpath_routing(src, dst, allpath);
  
  xgft_dmodk_routing(src, dst, path);
  
  //printf("%D-Mod-K path:\n");
  //xgft_print_path( path );
  
  high = xgft_h;
  
  for(i=0;i<XGFT_MAX_DISJOINT_PATH;i++)
    for(j=0;j<MAX_PATH_LEN;j++)
      djpath[i][j]=-1;
  
  while (src_label[high] == dst_label[high]) high --;
  for(i=0; i<high; i++) npath *= xgft_W[i];
  
  for(i=0; i<npath; i++){
    found = 1;
    for(j=0;j<2*xgft_h;j++){
      if(allpath[i][j] != path[j]){
        found = 0;
        break;
      }
    }
    if(found)break;
  }
  
  n = (npath < xgft_W[1])?npath:xgft_W[1];
  for(k=0; k < high-1; k++) m*= xgft_W[k];
  
  for(k=0; k < n; k++){
    memcpy(djpath[k],allpath[i],sizeof(Path));
    i = (i+m)%npath;
  }
  
  //printf("K = %d Disjoint paths:\n",k);
  //for(k=0; k < xgft_W[1]; k++)xgft_print_path( djpath[k] );
}

/*
void numPathpath_routing( int src, int dst, Path kpath[numPath]){
  int i;
  int ndjpath = xgft_W[2];
  int djpath[ XGFT_MAX_DISJOINT_PATH ] [ MAX_PATH_LEN ];
  
  if( ndjpath > XGFT_MAX_DISJOINT_PATH ){
    printf("Number of disjoint paths greater than XGFT_MAX_DISJOINT_PATH.\n");
    exit(0);
  }
  
  xgft_alldisjointpath_routing(src, dst, djpath);
  
  memcpy(kpath, djpath, numPath * sizeof(Path));
  
  //for(i=0; i < numPath; i++)xgft_print_path( kpath[i] );
}
*/

void xgft_kpath_routing( int src, int dst, Path kpath[numPath]){
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[MAX_H+1];
  int t_label[MAX_H+1];
  //int index[MAX_H];
  //int pathindex = 0;
  int i,j,l, npath = 1;
  int high;  
  int tdst,ttdst;
  int h; 
  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    kpath[0][0] = src;
    kpath[0][1] = -1;
    return;
  }

  high = h;
  tdst = dst;
  while (src_label[high] == dst_label[high]) high --;
  for(i=0; i<high; i++) npath *= xgft_W[i];

  // npath should be min ( npath, numPath, xgft_W[1]
  if( npath > numPath ) npath = numPath;
  if( npath > xgft_W[1] ) npath = xgft_W[1];
  
  //npath = (xgft_W[1] > numPath)?numPath:xgft_W[1];
  
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } else {
    for(i=0;i< numPath;i++)for(j=0;j<MAX_PATH_LEN;j++)kpath[i][j]=-1;
    
    for(i=0; i< npath; i++){
      kpath[i][0] = src;
      copy_label(r_label, src_label);
      r_label[0] = 1;
      r_label[1] = 0;
      compute_nodeid(&kpath[i][1], r_label);
      copy_label(t_label, r_label);
      ttdst = tdst;
      
      for (l=1; l<high; l++) {
        if( l == 1)
          t_label[l+1] = (ttdst + i) % xgft_W[l];
        else
          t_label[l+1] = ttdst % xgft_W[l]; 
        ttdst = ttdst / xgft_W[l];
        t_label[0] ++;
        compute_nodeid(&kpath[i][l+1], t_label);
        //printf("%d=>",kpath[i][l+1]);
      }
      
      //printf("\n");
      /* down path */
      for (l = high -1; l>=0; l--) {
        t_label[l+1] = dst_label[l+1];
        t_label[0] --;
        j = 2*high - l;
        compute_nodeid(&kpath[i][j], t_label);
      }
      kpath[i][2*high+1] = -1;
    }
  }
  
  //for(i=0; i<npath; i++) xgft_print_path(kpath[i]);
}

void xgft_dmodk_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[MAX_H+1];
  int i;
  int high;  
  int tdst;
  int h;

  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  high = h;
  tdst = dst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } else {
    /* printf("high = %d\n", high); */
    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
    /* up */
    for (i=1; i<high; i++) {
      /*      r_label[i+1] = dst_label[i] % W[i];  */
      r_label[i+1] = tdst % xgft_W[i]; 
      tdst = tdst / xgft_W[i];
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    /* down path */
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  }
  /*  print_path(path);
  //  if (check_path(path)) {
    //  printf("Path exists.\n");
  //  }
  */
}

void xgft_ecmp_routing(int src, int dst, Path kpath[numPath]){
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[_REAL_PATH_NUM][MAX_H+1];
  //int flag[600];
  int i,j,p,q;
  int high;  
  int tdst;
  int h;
  int npath,dport;

  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  
  for(i=0; i< numPath; i++)for(j=0; j < MAX_PATH_LEN; j++)kpath[i][j] = -1;
  
  if (src == dst){
    kpath[0][0] = src;
    kpath[0][1] = -1;
    return;
  }

  high = h;
  //tdst = dst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } else {
    npath = numPath;
    /* printf("high = %d\n", high); */
    for(p=0;p<npath;p++){
      tdst = dst;
      kpath[p][0] = src;
      copy_label(r_label[p], src_label);
      r_label[p][0] = 1;
      r_label[p][1] = 0;
      compute_nodeid(&kpath[p][1], r_label[p]);
      
      for (i=1; i<high; i++) {
        dport = tdst % xgft_W[i];
        q = rand() % numPath;
        //printf("kpath[%d][%d]:rand() = %d tdst = %d sport = %d\n",p,i,q,tdst,dport);
        r_label[p][i+1] = (dport + q)%xgft_W[i];
        tdst = tdst / xgft_W[i];
        r_label[p][0] ++;
        compute_nodeid(&kpath[p][i+1], r_label[p]);
      }
      /* down path */
      for (i = high -1; i>=0; i--) {
        r_label[p][i+1] = dst_label[i+1];
        r_label[p][0] --;
        compute_nodeid(&kpath[p][high-1 - i+high+1], r_label[p]);
      }
      kpath[p][high-1-i+high+1] = -1;
      //xgft_print_path( kpath[p] );
    }
  }
}

void xgft_model_dmodk_routing(int src, int dst, int *len, int *rsrc, int *rdst)
{
  int path[MAX_PATH_LEN];
  int i;

  xgft_dmodk_routing(src, dst, path);
  for (i=0; path[i+1] != -1; i++) {
    rsrc[i] = path[i];
    rdst[i] = path[i+1];
  }
  *len = i;
}

void xgft_random_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[MAX_H+1];
  int i;
  int high;  
  int tdst;
  int h;

  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  high = h;
  tdst = dst; tdst = tdst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } else {
    /* printf("high = %d\n", high); */
    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
    /* up */
    for (i=1; i<high; i++) {
      /*      r_label[i+1] = dst_label[i] % W[i];  */
      r_label[i+1] = rand() % xgft_W[i]; 
      /*      tdst = tdst / W[i]; */
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    /* down path */
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  }
  /*  print_path(path);
  //  if (check_path(path)) {
    //  printf("Path exists.\n");
  //  }
  */
}

/*
 // adaptive_routing computes the route based on the metrices stored
 // on load_graph
 */

void xgft_adaptive_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int r_label[MAX_H+1];
  int smallest_a[MAX_DEGREE+1];
  int smallest;
  int s_c;
  int i, j;
  int high;  
  int tdst;
  int h;

  h = xgft_h;

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  high = h;
  tdst = dst; tdst = tdst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);
  } else {
    /*    printf("high = %d\n", high); */
    path[0] = src;

    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;

    compute_nodeid(&path[1], r_label);
    /* up */

    for (i=1; i<high; i++) {

      smallest = 1000000;
      for (j=0; j<xgft_W[i]; j++) 
	if (load_graph[path[i]][j] < smallest) smallest = load_graph[path[i]][j];


      s_c = 0; 
      for (j=0; j<xgft_W[i]; j++) 
        if (load_graph[path[i]][j] == smallest) {
	  smallest_a[s_c] = j;
	  s_c++;
	}

      if (s_c == 0) { 
	printf("This cannot be right.\n");
	exit(0);
      }
      
      /*      r_label[i+1] = dst_label[i] % W[i];  */
      r_label[i+1] = smallest_a[rand() % s_c]; 
      /* r_label[i+1] = smallest_a[0]; 
      //      graph_r[path[i]][r_label[i+1]] ++;   */
      /* can't update the load graph, should change the load when others
         are set-up
      */
      /*      tdst = tdst / W[i]; */
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    /* down path */
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  }
  /*
  //  print_path(path);
  if (check_path(path)) {
    //  printf("Path exists.\n");
  }     */

}

void compute_pattern(int M, int W, int psrc[100][65], int pdst[100][65])
{
  int i, j;
  int mygcd;
  int ii = 0, jj = 0, kk = 0;

  if (M>W) mygcd = W;
  else mygcd = M;

  while ((W/mygcd*mygcd != W) || (M /mygcd * mygcd != M)) mygcd--;
 
  for (i=0; i<W; i++) {
    for (j=0; j<M; j++) {
     pdst[i][j] = ii; ii = (ii+1) % W;
      psrc[i][j] = jj;
      jj = (jj+1) % M;
      kk++;
      if (kk== M*W/mygcd) {
	jj = (jj+1) % M;
	kk = 0;
      }
    }
  }
  /*
  for (i=0; i<W; i++) {
    for (j=0; j<M-1; j++) {
      printf("(%d, %d), ", psrc[i][j], pdst[i][j]);
    }
    printf("(%d, %d)\n", psrc[i][j], pdst[i][j]);
  }
  exit(0);
  */
}
    


void xgft_roundrobin_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int top_label[MAX_H+1];
  int top_num;
  int r_label[MAX_H+1];
  int i;
  int high;  
  int tdst;
  int srcsw, dstsw;
  int h;
  int M[MAX_H];
  int W[MAX_H];
 
  h = xgft_h;
  for (i=0; i<MAX_H; i++) {M[i] = xgft_M[i]; W[i] = xgft_W[i];}

  if ((h > 3)) {
    printf("RR only supports 1-, 2-, and 3-level trees at this time. h=%d\n",
	   h);
    exit(0);
  }

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  srcsw = src / M[0];
  dstsw = dst / M[0];

  high = h;
  tdst = dst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);  
  } else if (high <=2) {
    // printf("high = %d\n", high);
    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
    // up
    for (i=1; i<high; i++) {
      //      r_label[i+1] = dst_label[i] % W[i]; 
      if (src/M[0] < dst/M[0])
	r_label[i+1] = (((dstsw - srcsw-1+2*srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      else 
	r_label[i+1] = ((((M[1]-1) - dstsw + srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      tdst = tdst / W[i];
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  } else if (high == 3) {
    
    int psrc[100][65];
    int pdst[100][65];
    int g_size = M[0]*M[1];
    int g_src = src_label[1] + src_label[2] * M[0];
    /* int g_dst = dst_label[1] * M[1] + dst_label[2]; */

    int mygcd, whichrow, ii;

    if (M[1] > W[1]) mygcd = W[1];
    else mygcd = M[1];
    while ((mygcd >0) && ((M[1]/mygcd*mygcd != M[1]) ||
			  (W[1]/mygcd*mygcd != W[1]))) mygcd--;
 
    if (mygcd <= 0) {
      printf("This is impossible. mygcd = %d\n", mygcd);
      exit(0);
    }

    compute_pattern(M[1], W[1], psrc, pdst);

    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);

    if (src_label[3] < dst_label[3]) {
      top_num = (((dst_label[3] + src_label[3]-1) * g_size + 
		  g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
      whichrow = top_num % W[1];
      for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
      if (ii == M[1]) {
	printf("This is impossible. 2\n");
	exit(0);
      }
      top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);
      /*
      top_num = ((dst_label[3] + src_label[3]-1) * g_size +
		 g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
      */
    } else {
      top_num = ((((M[2]-1) - dst_label[3] + src_label[3]) * g_size + 
		  g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
      whichrow = top_num % W[1];
      for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
      if (ii == M[1]) {
	printf("This is impossible. 2\n");
	exit(0);
      }
      top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);

    }

    /*
      top_num = (((M[2] -1) - dst_label[3] + src_label[3]-1) * g_size +
		 g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
    */

    top_label[0] = 3;
    top_label[1] = 0;
    top_label[2] = top_num % W[1];
    top_label[3] = top_num / W[1];
    top_label[4] = -1;
    /*
    printf("top num = %d, ", top_num);
    print_label(top_label); printf("\n");
    */

    r_label[0] = 2;
    r_label[2] = top_label[2];
    r_label[3] = src_label[3];
    compute_nodeid(&path[2], r_label);
    r_label[0] = 3;
    r_label[3] = top_label[3];
    compute_nodeid(&path[3], r_label);

    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
    
  }
  //  print_path(path);
}

void xgft_roundrobin1_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int top_label[MAX_H+1];
  int top_num;
  int r_label[MAX_H+1];
  int i;
  int high;  
  int tdst;
  int srcsw, dstsw;
  int h;
  int M[MAX_H];
  int W[MAX_H];
 
  h = xgft_h;
  for (i=0; i<MAX_H; i++) {M[i] = xgft_M[i]; W[i] = xgft_W[i];}

  if ((h > 3)) {
    printf("RR only supports 1-, 2-, and 3-level trees at this time. h=%d\n",
	   h);
    exit(0);
  }

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  srcsw = src / M[0];
  dstsw = dst / M[0];

  high = h;
  tdst = dst; tdst = tdst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);  
  } else if (high <=2) {
    // printf("high = %d\n", high);
    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
    // up
    /*
    for (i=1; i<high; i++) {
      //      r_label[i+1] = dst_label[i] % W[i]; 
      if (src/M[0] < dst/M[0])
	r_label[i+1] = (((dstsw - srcsw-1+2*srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      else 
	r_label[i+1] = ((((M[1]-1) - dstsw + srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      tdst = tdst / W[i];
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    */
  
    if (high == 2) {
      r_label[0] =2;
      {
	int src0, dst0, even;
	int nsrc0, ndst0, m0;
	
	src0 = src % M[0];
	dst0 = dst % M[0];
	even = M[0] / W[1] * W[1];
	m0 = M[0] - even;
	nsrc0 = src0 - even;
	ndst0 = dst0 - even;
	
	if ((dst0 < even)) { // d mod k
	  r_label[2] = dst0 % W[1];
	} else if ((src0 < even)) { // s mod k
	  r_label[2] = src0 % W[1];
	} else { // RR for left overs
	  if (src/M[0] < dst/M[0])
	    r_label[2] = (((dstsw - srcsw-1+2*srcsw) * m0 * m0) +
			  nsrc0 * m0 + ndst0) % W[1];
	  else 
	    r_label[2] = ((((M[1]-1) - dstsw + srcsw) * m0 * m0) +
			  nsrc0 * m0 + ndst0) % W[1];
	}
      }
      compute_nodeid(&path[2], r_label);
    }

    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  } else if (high == 3) {
    
    int psrc[100][65];
    int pdst[100][65];
    int g_size = M[0]*M[1];
    int g_src = src_label[1] + src_label[2] * M[0];
    /* int g_dst = dst_label[1] * M[1] + dst_label[2]; */

    int mygcd, whichrow, ii;

    if (M[1] > W[1]) mygcd = W[1];
    else mygcd = M[1];
    while ((mygcd >0) && ((M[1]/mygcd*mygcd != M[1]) ||
			  (W[1]/mygcd*mygcd != W[1]))) mygcd--;
 
    if (mygcd <= 0) {
      printf("This is impossible. mygcd = %d\n", mygcd);
      exit(0);
    }

    compute_pattern(M[1], W[1], psrc, pdst);

    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
   
    if (src_label[3] < dst_label[3]) {
      top_num = (((dst_label[3] + src_label[3]-1) * g_size + 
		  g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
      whichrow = top_num % W[1];
      for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
      if (ii == M[1]) {
	printf("This is impossible. 2\n");
	exit(0);
      }
      top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);
      /*
      top_num = ((dst_label[3] + src_label[3]-1) * g_size +
		 g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
      */
    } else {
      top_num = ((((M[2]-1) - dst_label[3] + src_label[3]) * g_size + 
		  g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
      whichrow = top_num % W[1];
      for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
      if (ii == M[1]) {
	printf("This is impossible. 2\n");
	exit(0);
      }
      top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);

    }

    /*
      top_num = (((M[2] -1) - dst_label[3] + src_label[3]-1) * g_size +
		 g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
    */

    top_label[0] = 3;
    top_label[1] = 0;
    top_label[2] = top_num % W[1];
    top_label[3] = top_num / W[1];
    top_label[4] = -1;
    /*
    printf("top num = %d, ", top_num);
    print_label(top_label); printf("\n");
    */

    r_label[0] = 2;
    r_label[2] = top_label[2];
    r_label[3] = src_label[3];
    compute_nodeid(&path[2], r_label);
    r_label[0] = 3;
    r_label[3] = top_label[3];
    compute_nodeid(&path[3], r_label);

    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
    
  }
  //  print_path(path);
}



void xgft_roundrobin2_routing(int src, int dst, int *path)
{
  int src_label[MAX_H+1];
  int dst_label[MAX_H+1];
  int top_label[MAX_H+1];
  int top_num;
  int r_label[MAX_H+1];
  int i;
  int high;  
  int tdst;
  int srcsw, dstsw;
  int h;
  int M[MAX_H];
  int W[MAX_H];
 
  h = xgft_h;
  for (i=0; i<MAX_H; i++) {M[i] = xgft_M[i]; W[i] = xgft_W[i];}

  if ((h > 3)) {
    printf("RR only supports 1-, 2-, and 3-level trees at this time. h=%d\n",
	   h);
    exit(0);
  }

  compute_label(src, src_label);
  compute_label(dst, dst_label);

  if (src_label[0] != 0) {
    printf("Node %d, ", src);
    print_label(src_label);
    printf(" is not a processing node.\n");
    exit(0);
  }

  if (dst_label[0] != 0) {
    printf("Node %d, ", dst);
    print_label(dst_label);
    printf(" is not a processing node.\n");
    exit(0);
  }
  if (src == dst){
    path[0] = src;
    path[1] = -1;
    return;
  }

  srcsw = src / M[0];
  dstsw = dst / M[0];

  high = h;
  tdst = dst; tdst = tdst;
  while (src_label[high] == dst_label[high]) high --;
  if (high <= 0) { 
    printf("This can't be right.\n");
    exit(0);  
  } else if (high <=2) {
    // printf("high = %d\n", high);
    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);
    // up
    /*
    for (i=1; i<high; i++) {
      //      r_label[i+1] = dst_label[i] % W[i]; 
      if (src/M[0] < dst/M[0])
	r_label[i+1] = (((dstsw - srcsw-1+2*srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      else 
	r_label[i+1] = ((((M[1]-1) - dstsw + srcsw) * M[0] * M[0]) +
			src % M[0]*M[0] + dst % M[0]) % W[i];
      tdst = tdst / W[i];
      r_label[0] ++;
      compute_nodeid(&path[i+1], r_label);
    }
    */
  
    if (high == 2) {
      r_label[0] =2;
      {
	int src0, dst0, even;
	int nsrc0, ndst0, m0;
	
	src0 = src % M[0];
	dst0 = dst % M[0];
	even = M[0] / W[1] * W[1];
	m0 = M[0] - even;
	nsrc0 = src0;
	ndst0 = dst0 - even;
	
	if ((dst0 < even)) { // d mod k
	  r_label[2] = dst0 % W[1];
	} else if ((src0 < even)) { // s mod k
		r_label[2] = src0 % W[1]; 
	} else { // RR for left overs
	  if (src/M[0] < dst/M[0])
	    r_label[2] = (((dstsw - srcsw-1+2*srcsw) * m0 * m0) +
			  nsrc0 * m0 + ndst0) % W[1];
	  else 
	    r_label[2] = ((((M[1]-1) - dstsw + srcsw) * m0 * m0) +
			  nsrc0 * m0 + ndst0) % W[1];
	}
      }
      compute_nodeid(&path[2], r_label);
    }

    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
  } else if (high == 3) {
    
    int psrc[100][65];
    int pdst[100][65];
    int g_size = M[0]*M[1];
    int g_src = src_label[1] + src_label[2] * M[0];
    /*    int g_dst = dst_label[1] * M[1] + dst_label[2]; */
    int even;

    int mylcm, whichrow, ii;

    if (M[1] > W[2]) mylcm = M[1];
    else mylcm = W[2];
    while ((mylcm/M[1]*M[1] != mylcm) ||
	   (mylcm/W[2]*W[2] != mylcm)) mylcm++;
 
    even  = mylcm/M[1]*W[1];

    compute_pattern(M[1], W[1], psrc, pdst);

    path[0] = src;
    copy_label(r_label, src_label);
    r_label[0] = 1;
    r_label[1] = 0;
    compute_nodeid(&path[1], r_label);


    if (M[0] < even) { // not enough PEs, RRR for all
      if (src_label[3] < dst_label[3]) {
	top_num = (((dst_label[3] + src_label[3]-1) * g_size + 
		    g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
	whichrow = top_num % W[1];
	for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
	if (ii == M[1]) {
	  printf("This is impossible. 2\n");
	  exit(0);
	}
	top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);
	/*
	  top_num = ((dst_label[3] + src_label[3]-1) * g_size +
	  g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
	*/
      } else {
	top_num = ((((M[2]-1) - dst_label[3] + src_label[3]) * g_size + 
		    g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
	whichrow = top_num % W[1];
	for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
	if (ii == M[1]) {
	  printf("This is impossible. 2\n");
	  exit(0);
	}
	top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);	
      }
      top_label[0] = 3;
      top_label[1] = 0;
      top_label[2] = top_num % W[1];
      top_label[3] = top_num / W[1];
      top_label[4] = -1;
    } else {
      if (dst_label[1] < even) {
	/* d_mod_k */
        int ndst = dst_label[3]*M[1]*even + dst_label[2]*even + dst_label[1];
	top_label[0] = 3;
	top_label[1] = 0;
        top_label[2] = ndst % W[1];
        top_label[3] = (ndst / W[1]) % W[2];
	top_label[4] = -1;        
      } else if (src_label[1] < even) {
	/* s-mod-k */
        int ndst = src_label[3]*M[1]*even + src_label[2]*even + src_label[1];
	top_label[0] = 3;
	top_label[1] = 0;
        top_label[2] = ndst % W[1];
        top_label[3] = (ndst / W[1]) % W[2];
	top_label[4] = -1;        
      } else {
        /* RRR */
	int g_size = (M[0]-even)*M[1];
	int g_src = (src_label[1] - even) + src_label[2] * (M[0]-even);
	/* int g_dst = (dst_label[1] -even) * M[1] + dst_label[2]; */

	if (src_label[3] < dst_label[3]) {
	  top_num = (((dst_label[3] + src_label[3]-1) * g_size + 
		      g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
	  whichrow = top_num % W[1];
	  for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
	  if (ii == M[1]) {
	    printf("This is impossible. 2\n");
	    exit(0);
	  }
	  top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);
	  /*
	    top_num = ((dst_label[3] + src_label[3]-1) * g_size +
	    g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
	  */
	} else {
	  top_num = ((((M[2]-1) - dst_label[3] + src_label[3]) * g_size + 
		      g_src % g_size * g_size + dst_label[1]*M[1])) / M[1];
	  whichrow = top_num % W[1];
	  for(ii=0; (psrc[whichrow][ii] != dst_label[2]) && (ii < M[1]); ii++);
	  if (ii == M[1]) {
	    printf("This is impossible. 2\n");
	    exit(0);
	  }
	  top_num = (top_num*M[1] + ii) % (W[0]*W[1]*W[2]);	
	}
	top_label[0] = 3;
	top_label[1] = 0;
	top_label[2] = top_num % W[1];
	top_label[3] = top_num / W[1];
	top_label[4] = -1;
      }
    }
    /*
      top_num = (((M[2] -1) - dst_label[3] + src_label[3]-1) * g_size +
		 g_src % g_size * g_size + g_dst % g_size) % (W[0]*W[1]*W[2]);
    */

    /*
    printf("top num = %d, ", top_num);
    print_label(top_label); printf("\n");
    */

    r_label[0] = 2;
    r_label[2] = top_label[2];
    r_label[3] = src_label[3];
    compute_nodeid(&path[2], r_label);
    r_label[0] = 3;
    r_label[3] = top_label[3];
    compute_nodeid(&path[3], r_label);

    // down path
    for (i = high -1; i>=0; i--) {
      r_label[i+1] = dst_label[i+1];
      r_label[0] --;
      compute_nodeid(&path[high-1 - i+high+1], r_label);
    }
    path[high-1-i+high+1] = -1;
    
  }
  //  print_path(path);
}



int xgft_check_path(int *path)
{
  int i, j;
  for (i=0; path[i+1] != -1; i++) {
    for (j=0; graph[path[i]][j] != -1; j++)
      if (graph[path[i]][j] == path[i+1]) break;
    if (graph[path[i]][j] == -1) {
      xgft_print_path(path);
      printf(" does not exist!!");
      exit(0);
      return FALSE;
    }
  }
  return TRUE;
}

/* static int histogram[10000000]; */

void xgft_topology_init(int h, int *m, int *w, long long int *bw, int r) 
{
  int i, j;



  xgft_h = h;
  if (MAX_H < xgft_h) {
    printf("MAX_H %d is too small, less than xgft_h %d\n", MAX_H, xgft_h);
    exit(0);
  }

  for(i=0; i<xgft_h; i++) {
    xgft_M[i] = m[i];
    xgft_W[i] = w[i];
    xgft_BW[i] = bw[i];
  }
  compute_baseL(h, m, w);
  xgft_build_topology(h, m, w, bw);

  topology = XGFT;

  /* reset the load graph when topology is initialized */
  for (i=0; i<totNode; i++) {
    for (j=0; graph[i][j] != -1; j++)
      load_graph[i][j] = 0;
  }

  printf("Simulating XGFT with parameters: height = %d\n  m = [", h);
  for (i=0; i<h-1; i++) printf("%d, ", m[i]);
  printf("%d]\n  w = [", m[i]);

  for (i=0; i<h-1; i++) printf("%d, ", w[i]);
  printf("%d]\n  bw = [", w[i]);

  for (i=0; i<h-1; i++) printf("%lld, ", bw[i]);
  printf("%lld]\n", bw[i]);

  for (i=0; i<=h; i++) {
    printf("sizeL[%d] = %d\n", i, sizeL[i]);
  }
  for (i=0; i<=h; i++) {
    printf("baseL[%d] = %d\n", i, baseL[i]);
  }

  printf("%d PEs, %d switches, %d nodes (switch+PE)\n", 
	 totPE, totNode-totPE, totNode);
  
  routing = r;
  if (r == XGFT_DMODK_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_dmodk_routing;
    model_routing_algorithm = xgft_model_routing_algorithm = xgft_model_dmodk_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: Destination mod K routing\n");
  } else if (r == XGFT_RANDOM_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_random_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: Random routing\n");
  } else if (r == XGFT_ADAPTIVE_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_adaptive_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: Adaptive routing\n");
  } else if (r == XGFT_RR_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_roundrobin_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: RR routing\n");
  } else if (r == XGFT_RR1_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_roundrobin1_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: RR1 routing\n");
  } else if (r == XGFT_RR2_ROUTING) {
    routing_algorithm = xgft_routing_algorithm = xgft_roundrobin2_routing;
    routingType = SINGLEPATH;
    printf("Routing algorithm: RR2 routing\n");
  } else if (r == XGFT_KPATH_ROUTING) {
    routingType = MULTIPATH;
    printf("Routing algorithm: KPATH\n");
  } else if (r == XGFT_ECMP_ROUTING) {
    routingType = MULTIPATH;
    printf("Routing algorithm: ECMP\n");
  }else if(r == XGFT_VLB_ROUTING){
     printf("Routing algorithm: VLB\n");
  }else {
    printf("Routing scheme %d not supported on xgft.\n", routing);
    exit(0);
  }

  /*
  for (i=0; i<MAX_NODE; i++)
    for (j=0; j<MAX_DEGREE; j++)
      graph_c[i][j] = 0;

  
  xgft_print_topology();
   while (1) {
    printf("src dst = "); fflush(0);
    k = scanf("%d %d", &i, &j); getchar();
    routing_algorithm(i, j, path);
    xgft_print_path(path);
  }
  */

  /*
  for (i=0; i<totPE; i++)
    for (j=0; j<totPE; j++) {
      int path[MAX_PATH_LEN];
      count ++;
      if (count % 1000000 == 0) {
	printf("count = %d/%d\n", count, totPE*(totPE-1));
	fflush(0);
      } 
      if (i != j) {
	routing_algorithm(i, j, path);

	for (k=1; path[k] != -1; k++) {
	  for (l=0; (graph[path[k-1]][l] != path[k]) && (l<MAX_DEGREE); l++);
	  if (l >= MAX_DEGREE) {
	    printf("(%d %d) k-1= %d, path[k-1] = %d, path[k] = %d Impossible.\n", i, j, k-1,
		   path[k-1], path[k]);
            xgft_print_path(path);
            for (l=0; graph[path[k-1]][l] != -1; l++)
	      printf("graph[%d][%d] = %d\n", path[k-1],l, graph[path[k-1]][l]);
	    exit(0);
	  }
	  graph_c[path[k-1]][l] ++;
	}
      }
    }

  for (i=0; i<10000000; i++) histogram[i] = 0;

  for (i=0; i<totNode; i++)
    for (j=0; (j<MAX_DEGREE) && (graph[i][j] != -1); j++) {
      if (graph_c[i][j] >= 10000000) {
	printf("histogram table too small. %d\n", graph_c[i][j]);
	exit(0);
      }
      
      histogram[graph_c[i][j]] ++;
    }
  
  count = 0;
  count1 = 0;
  for (i=1; i<10000000; i++) {
    if (histogram[i] != 0) {
      printf("%d: %d\n", i, histogram[i]);
      count += histogram[i];
      count1 += ((long long)i) * histogram[i];
    }
  }	
  
  printf("total = %d, %d linked used by %lld sd pairs\n", count, count1, 
	 ((long long)totPE) *(totPE-1));
  */
  if (DEBUG_LEVEL > 0) {
    xgft_print_topology();
  }  
}
  



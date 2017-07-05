
#include "topology.h"


#ifdef _SMALL
#define MAX_ITEM 40000LL
#else
#define MAX_ITEM 400000000LL
#endif


int node_level[MAX_NODE];                   // populated in calculate_node_levels()
int label_array[MAX_NODE][MAX_H+1];// precomputed labels for each nodeid
double rate_allocation_vector[MAX_ITEM];  //only used in  NLP2 algorithm

//calculates 'labels' for each node once and stores for future lookup
void  calculate_node_levels(){
  long long int i; 
  int label[MAX_H+1];                   
  for (i=0;i<totNode;i++){
    compute_label(i, label);
    node_level[i] = label[0];
  }
}            

int fan_out(int node_id){
          int fanout=1;     
          switch(node_level[node_id]){
            case 0: fanout = 1; 
                    break;
            case 1: fanout = xgft_W[1];
                    break;
            case 2: fanout = xgft_W[1] * xgft_W[2];              
                    break;
            case 3: fanout = xgft_W[1]*xgft_W[2]*xgft_W[3];
                    break;
            default: printf("Incorrect Node level detected");
          }
  return fanout;

}


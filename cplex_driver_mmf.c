// Seperate diver to simulate theoretical optimal(non-MMF) throughput from models
//use 'make gen' and 'cplex_gen.x' executable

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "topology_var.h"

extern void calculate_node_levels();

extern long long int generic_cplex_mmf_from_trace_file(char *filename, int* iteration_count, double * exec_time);

extern long long int xgft_mmf_cplex_from_trace_file(char *filename, int* iteration_count, double *exec_time);
extern long long int xgft_mmf_nonlp_from_trace_file (char* filename, int* iteration_count, double *exec_time);
extern long long int xgft_mmf_dmodk_from_trace_file (char* filename, int* iteration_count, double *exec_time);
extern long long int xgft_mmf_ppf_from_trace_file (char* filename, int* iteration_count, double *exec_time);

//in a separate file
extern void xgft_mmf_OPT (char* filename, double* throughput, int* iteration_count, double *exec_time);
extern double xbar_mmf_nlp_var2(char* filename, int* flow_count );


extern double xbar_from_trace_file(char* filename, int *flow_count);
extern double xbar_mmf_from_trace_file(char* filename, int* flow_count );



struct timeval t;
int all2allpattern_input=0;
int MCF_FLAG=0;				// 0 = MMF(Max-min fair), 1 = MCF(max concurrent flow)

int main(int argc, char *argv[]) 
{
  int i;
  FILE *ofd;
  int traffic_count;
  int iteration_count=0;
  double exec_time=0.0;
  long long  aggr_throughput=0;
  double  xbar_throughput=0;
  double  opt_patam_throughput=0;
  if (argc < 3) {
    printf("Usage: \n");
      printf("Usage: ./fast_mmf.out [xgft/pgft] h m1 m2 ... mh w1 w2 ...wh p1 p2 ... ph [GEN/LP/PF/DMODK/OPT] trafficfile \n");
    exit(0);
  }

  if ( (strcmp(argv[1], "xgft") == 0 )|| (strcmp(argv[1],"pgft")==0)) {
    int h;
    int M[MAX_H];
    int W[MAX_H];
    long long int BW[MAX_H];

    h = atoi(argv[2]);

    if (argc < 3*h+5 || argc > 3*h+7) {
      printf("Usage: ./fast_mmf.out xgft h m1 m2 ... mh w1 w2 ...wh p1 p2 ... ph [GEN/LP/PF/DMODK/OPT] trafficfile \n");
      exit(0);
    }
    if (argc >=3*h+6 && argv[3*h+5]!=0){
      printf("CALCULATING CONCURRENT FLOW INSTEAD OF MMF FLOW \n");
      MCF_FLAG=1;
    }
    if (argc==3*h+7 && argv[3*h+6]!=0){
      printf("Processing All-to-all Traffic pattern, not to be invoked from utilities/get_index.x \n");
      all2allpattern_input=1;
    }

    for (i=0; i<h; i++) M[i]  = atoi(argv[2+i+1]);
    for (i=0; i<h; i++) W[i]  = atoi(argv[2+h+i+1]);
    for (i=0; i<h; i++) BW[i] = atoi(argv[2+2*h+i+1])*((long long int)1000000);

    xgft_topology_init(h, M, W, BW, XGFT_KPATH_ROUTING);

    xbar_throughput =  xbar_mmf_from_trace_file(argv[3*h+4], &traffic_count);

    if(strcmp(argv[3*h+3], "GEN")==0)
      aggr_throughput = generic_cplex_mmf_from_trace_file(argv[2+3*h+2], &iteration_count, &exec_time);    
    else if(strcmp(argv[3*h+3], "LP")==0)
      aggr_throughput = xgft_mmf_cplex_from_trace_file(argv[2+3*h+2], &iteration_count, &exec_time);
    else if(strcmp(argv[3*h+3], "NON_LP")==0){ // how is it different from PF?
      calculate_node_levels();				//calculates 'labels' for each node once and stores for future lookup
      aggr_throughput = xgft_mmf_nonlp_from_trace_file(argv[2+3*h+2], &iteration_count, &exec_time);
    }
    else if( (strcmp(argv[3*h+3], "PF")==0) || (strcmp(argv[3*h+3], "NLP")==0) ){
      calculate_node_levels();				//calculates 'labels' for each node once and stores for future lookup
      aggr_throughput = xgft_mmf_ppf_from_trace_file(argv[2+3*h+2], &iteration_count, &exec_time);
    }
    else if( (strcmp(argv[3*h+3], "DMODK")==0) || (strcmp(argv[3*h+3], "DMK")==0) )
      aggr_throughput = xgft_mmf_dmodk_from_trace_file(argv[2+3*h+2], &iteration_count, &exec_time);
    else if(strcmp(argv[3*h+3], "OPT")==0){
      xgft_mmf_OPT(argv[3*h+4],&opt_patam_throughput, &iteration_count, &exec_time);
      xbar_throughput = xbar_mmf_nlp_var2(argv[3*h+4], &traffic_count);
      aggr_throughput = (long long int) opt_patam_throughput;
    }
    else {
      printf("Unknown formulation method %s \n", argv[3*h+3] ); 
      exit(0);
    }


    printf("\nAggregate throughput for traffic %s,  method %s is %lld\n", argv[3*h+4], argv[3*h+3], aggr_throughput);
    printf("Crossbar network throughput is %lf\n", xbar_throughput);

    printf("%s xput\t\t Crossbar \t LFTI\t Iteration#\t Time\n", argv[3*h+3]);
    printf( "%8.6lf\t %8.6lf\t %8.6lf\t %d\t %8.6lf (s)\n",
     aggr_throughput*1.0 / traffic_count,
     xbar_throughput*1.0/traffic_count, 
     aggr_throughput*1.0/xbar_throughput,
     iteration_count, 
     exec_time);

    /*char outfilename[100];
    strcpy(outfilename, argv[3*h+4]);
    strcat(outfilename, argv[3*h+3]);
    strcat(outfilename, ".out");
    if ((ofd = fopen(outfilename, "w")) == NULL) {
      printf("Can't open model_result for write.\n");
      exit(0);
    }
    fprintf(ofd,"%s throughput\t\t crossbar throughput\t LFTI\t Iteration#\t Time\n", argv[3*h+3]);
    fprintf(ofd, "%8.6lf\t %8.6lf\t %8.6lf\t %d\t %8.6lf (s)\n",
     xbar_throughput*1.0/traffic_count, aggr_throughput*1.0 / traffic_count,
     aggr_throughput*1.0/xbar_throughput,iteration_count, time_in_sec);
    if( ofd )fclose(ofd); 
    */

 }
 else {
    printf("ERROR: This program only supports xgft .\n");
    exit(0);
  }


  //Export throughput to file model_result in a format consistent with model engine of ../model/model_engine.c
  // consistency required to calculate LFTI using unaltered ../utilities/get_index.x
  //NOTE: the order of xbar and aggr throughput in model_result is opposite to the one in other models
  //should i change the model, or update get_index.c?
  if ((ofd = fopen("model_result", "w")) == NULL) {
    printf("Can't open model_result for write.\n");
    exit(0);
  }
  fprintf(ofd, "%8.6lf %8.6lf %8.6lf %8.6lf %8.6lf %8.6lf %d %8.6lf %8.6lf %8.6lf %8.6lf\n",
  0.0,
  0.0,
  0.0,
  0.0,
  xbar_throughput*1.0/traffic_count,
  aggr_throughput*1.0 / traffic_count,
  iteration_count,
  exec_time, 0.0, 0.0, 0.0);
  
  if( ofd )fclose(ofd); 
  
  return 0;
}



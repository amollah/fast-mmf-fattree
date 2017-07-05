// Seperate diver to simulate theoretical optimal(non-MMF) throughput from models
//use 'make gen' and 'cplex_gen.x' executable

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "topology_var.h"

extern void xgft_mmf_NLP2 (char* filename, double* throughput, int* iteration_count, double *exec_time);
extern double xbar_mmf_nlp(char* filename, int* flow_count );
extern double xbar_mmf_nlp_var2(char* filename, int* flow_count );


struct timeval t;
int all2allpattern_input=0;
int MCF_FLAG=0;				// 0 = MMF(Max-min fair), 1 = MCF(max concurrent flow)

int main(int argc, char *argv[])
{
  struct timeval start, end;
  int i;
  FILE *ofd;
  int traffic_count;
  int iteration_count=0;
  double exec_time=0.0;
  double aggr_throughput=0;
  double  xbar_throughput=0;
  if (argc < 2) {
    printf("Usage: \n");
    printf("  ./a.out xgft h m1 m2 ... mh w1 w2 ...wh bw1 bw2 ... bwh traffic_file NLP2\n");

    exit(0);
  }

 if (strcmp(argv[1], "xgft") == 0) {

    int h;
    int M[MAX_H];
    int W[MAX_H];
    long long int BW[MAX_H];

    h = atoi(argv[2]);
    for (i=0; i<h; i++) M[i]  = atoi(argv[2+i+1]);
    for (i=0; i<h; i++) W[i]  = atoi(argv[2+h+i+1]);
    for (i=0; i<h; i++) BW[i] = atoi(argv[2+2*h+i+1])*((long long int)1000000);


    xgft_topology_init(h, M, W, BW, XGFT_KPATH_ROUTING);
    if(strcmp(argv[3*h+3], "NLP2")==0){
      printf("Entered NLP2 routine!!!!\n\n\n");
      xgft_mmf_NLP2(argv[3*h+4],&aggr_throughput, &iteration_count, &exec_time);
      xbar_throughput = xbar_mmf_nlp_var2(argv[3*h+4], &traffic_count);
    }
    else{printf("XBARv1 ONLY"); aggr_throughput=-1;
      xbar_throughput = xbar_mmf_nlp(argv[3*h+4], &traffic_count);

}
    printf("\nAggregate throughput for traffic %s,  method %s is %lf\n", argv[3*h+4], argv[3*h+3], aggr_throughput);
    printf("Crossbar network throughput is %lf\n", xbar_throughput);

    double time_in_sec = ((end.tv_sec - start.tv_sec)+ (end.tv_usec - start.tv_usec)/1000000.0);
    printf("%s throughput\t crossbar throughput\t LFTI\t Iteration#\t Time\n", argv[3*h+3]);
    printf( "%8.6lf\t %8.6lf\t %8.6lf\t %d\t %8.6lf (s)\n",
     aggr_throughput*1.0/traffic_count, xbar_throughput*1.0 / traffic_count,
     aggr_throughput*1.0/xbar_throughput,iteration_count, exec_time);

    char outfilename[100];
    strcpy(outfilename, argv[3*h+4]);
    strcat(outfilename, argv[3*h+3]);
    strcat(outfilename, ".out");
    if ((ofd = fopen(outfilename, "w")) == NULL) {
      printf("Can't open model_result for write.\n");
      exit(0);
    }
    fprintf(ofd,"%s throughput\t crossbar throughput\t LFTI\t Iteration#\t Time\n", argv[3*h+3]);
    fprintf(ofd, "%8.6lf\t %8.6lf\t %8.6lf\t %d\t %8.6lf (s)\n",
     aggr_throughput*1.0/traffic_count, xbar_throughput*1.0 / traffic_count,
     aggr_throughput*1.0/xbar_throughput,iteration_count, time_in_sec);
    if( ofd )fclose(ofd);

  }
  else {
    printf("ERROR: This program only supports xgft\n");
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



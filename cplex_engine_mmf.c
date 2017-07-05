#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "topology.h"
#include <unistd.h>
#include <sys/time.h>


#define BW_SATURATION_GAP 1

#ifdef _SMALL
#define MAX_ITEM 40000LL
#else
#define MAX_ITEM 400000000LL
#endif


char outputfilename[100];			//to store the output filename
static int sdpaths[MAX_ALLPATH][MAX_PATH_LEN];  //used for VLB
int path[MAX_PATH_LEN];				//used for  dmodk

extern int all2allpattern_input;	//for special case handling of large/all2all pattern file
extern int MCF_FLAG;	//flag to choose between concurrent flow metric and max-min fair metric


//Prototypes
void all2all_update_link(char *filename, int sw,int port, int update_load);


//External function prototypes
double timediff(struct timeval start, struct timeval end);
int fan_out(int node_id);
void  calculate_node_levels();
void init_list();
long long int listmalloc();
void runcmd(char *);

//returns elapsed time of the given interval in seconds
/*double timediff(struct timeval start, struct timeval end)
{
 return ((end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec))/1000000.0;
}
*/

//! This function initializes the graph_m array.
void dfly_cplex_init() {
}

//! This function initializes the graph_m array.
//called by jfish_cplex formulation only
void model_init() {
  int i, j;

  for (i=0; i<totNode; i++) {
    for (j=0; graph[i][j] != -1; j++) {
      graph_m[i][j] = 0;
    }
  }
}


//! This function prints the path for a source and a destination
static void print_allpath(int s, int d,
			  int allpath[MAX_ALLPATH][MAX_PATH_LEN],
			  int num_path){
  int i, j;
  printf("\n%d paths between %d and %d:\n",num_path, s,d);
  for (i=0; i<num_path; i++) {
    printf("Path[%d]: ", i); fflush(0);
    printf("%d", allpath[i][0]);
    for(j=0; allpath[i][j+1] != -1; j++)
      printf("->%d", allpath[i][j+1]);
    printf("\n");
  }
}

//! This function initializes the paths to -1
void init_allpath(int allpath[MAX_ALLPATH][MAX_PATH_LEN]){
  int i,j;
  for(i=0; i < MAX_ALLPATH; i++){
    for(j=0; j< MAX_PATH_LEN; j++){
      allpath[i][j] = -1;
    }
  }
}


// a large linkedlist


extern long long int value[MAX_ITEM];
extern long long int next[MAX_ITEM];
extern long long int head;


static long long int  bandwidth_mmf[MAX_NODE][MAX_DEGREE];
static long long int load[MAX_NODE][MAX_DEGREE];
static long long int flow_vector[MAX_ITEM];
static long long int dirty_flow[MAX_ITEM];


static long long int path_per_flow[MAX_ITEM];

extern int node_level[MAX_NODE];		   // populated in calculate_node_levels()

static double temp_flow_vector[MAX_ITEM];          //only used in NON_LP algorithm..to prevent precision error
static long long int final_flow_vector[MAX_ITEM];  //only used in  NON_LP algorithm


//calculates 'labels' for each node once and stores for future lookup
/*
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
*/

/**
  * Utility function
  * @param cmd any shell command string
  * passes unix command to system() for execution
  */
/*void runcmd(char *cmd)
{
  int ret;
  fflush(0);
  ret = system(cmd);
  fflush(0);
  if ((ret == -1)) {
    printf("Cmd: '%s' failed. %d, %d\n", cmd, ret, ret);
    exit(0);
  }
}*/

/**
  * MMF output function
  * @param n number of flows
  * During/After calculating the max-min-fair rates, dumps the rate allocation status
  * reads from global array flow_vector (used in most/all formulations)
  * reads from global array  dirty_flow (used in XGFT_LP formulation)
 */
void dumpvector(long long int n){
  long long int i;
  printf("FLOW VECTOR: \n");
  for(i=0;i<n;i++)
  printf("[%lld]: %lld\n", i, flow_vector[i]);
  printf("DIRTY: \n");
  for(i=0;i<n;i++)
  printf("[%lld]: %lld\n", i, dirty_flow[i]);

}


/**
  * Linked list routines
  * init_list(): initialized the global arrays: value[] and next[]
  * Array implementation of link-list, same arrays(value and next) used to store multiple list
  */
/*void init_list(){
  long long int i;
  head = 0;
  for (i=0; i<MAX_ITEM-1; i++) {
    next[i] = i+1;
    value[i] = MAX_ITEM;
  }
  next[MAX_ITEM-1] = MAX_ITEM;
  value[MAX_ITEM-1] = MAX_ITEM;
}
*/

/**
  * Linked list routines
  * LEGACY..NEVER REFERENCED
  */

/*
void listfree(long long int index)
{
  next[index] = head;
  value[index] = MAX_ITEM;
  head = index;
}
*/


/**
  * Linked list routines
  * reads up the global variable head (used like a 'free' pointer) and returns its value
  * @return index to the next available/usable data slot(lowest index from available slots)
  */
/*
long long int listmalloc(){
  long long int res;

  if (head >= MAX_ITEM) {
    printf("List out of memory.\n");
    exit(0);
  }
  res = head;
  head = next[head];
  return res;
}
*/
/**
 * node-to-flow map structures(source-to-flow and destination-to-flow)
 * indexes(nodeid) used as keys, chain of associated flows returned as values
 * srchead: returns the head of the list of flows/SDpairs originating from nodeid
 * srctail: returns the tail of the list of flows/SDpairs originating from nodeid
 * dsthead: returns the head of the list of flows/SDpairs destined to nodeid
 * dsttail: returns the tail of the list of flows/SDpairs destined to nodeid
*/
long long int srchead[MAX_NODE];
long long int srctail[MAX_NODE];
long long int dsthead[MAX_NODE];
long long int dsttail[MAX_NODE];

/**
  * Initialization routine
  *  for source-to-flow and destination-to-flow maps
  */
void init_src_dst_list()
{
  int i;
  for (i=0; i<totNode; i++) {
    srchead[i] = srctail[i] = dsthead[i] = dsthead[i] = MAX_ITEM;
  }
}

/**
  * insert item at the tail of source-flow map
  * @param value the (flow id) to be mapped/inserted
  * @param ii the node id/src to be mapped with
  */
void insertsrc(long long int val, long long int ii)
{
  long long int i;
  i = listmalloc();
  value[i] = val;
  next[i] = MAX_ITEM;

  if (srchead[ii] == MAX_ITEM) {
    srchead[ii] = srctail[ii] = i;
  } else {
    next[srctail[ii]] = i;
    srctail[ii] = i;
  }
}

void insertdst(long long int val, long long int ii)
{
  long long int i;
  i = listmalloc();
  value[i] = val;
  next[i] = MAX_ITEM;
  if (dsthead[ii] == MAX_ITEM) {
    dsthead[ii] = dsttail[ii] = i;
  } else {
    next[dsttail[ii]] = i;
    dsttail[ii] = i;
  }
}

#define MAX_SWITCH 4000

long long int swhead[MAX_SWITCH][MAX_DEGREE];
long long int swtail[MAX_SWITCH][MAX_DEGREE];

void init_sw_list()
{
  int i, j;
  for (i=0; i<totNode - totPE; i++) 
    for (j=0; j < MAX_DEGREE; j++) {
      swhead[i][j] = swtail[i][j] = MAX_ITEM;
  }
}

void insertsw(long long int val, int ii, int port)
{
  long long int i;
  i = listmalloc();
  value[i] = val;
  next[i] = MAX_ITEM;
  if (swhead[ii][port] == MAX_ITEM) {
    swhead[ii][port] = swtail[ii][port] = i;
  } else {
    next[swtail[ii][port]] = i;
    swtail[ii][port] = i;
    //printf("## %lld\n", i);
  }
}

//! This function computes the model index using the traffic 
//! from the trace_file



//iteration value returned
long long int generic_cplex_mmf_from_trace_file(char *filename,  int* iteration_count, double *exec_time){

  struct timeval t1,t2;
  FILE *fd, *ofd;
  char buf[1000], *ch;
  long long int i, j, k;
  long long int tot_path;
  long long int ii, jj, kk;
  long long int counter;
  int next_s, next_d;
  long long int cc=0, dd=0;

  /*Variables/storage specific to the iterative algorithm*/
  FILE *cplex_fd;			//file descriptor to read in LP solution
  char cmd[1000];			//command string to invoke CPLEX from shell
  double lp_solution=0;			// solution objective value
  int iterator=0;
  int all_saturated=0;			//boundary condition flag

  //initilizing bandwidth input to LP formulation
  gettimeofday(&t1, NULL);

  for(i=0;i<MAX_NODE;i++){
      for(j=0;j<MAX_DEGREE;j++){
            bandwidth_mmf[i][j] = bandwidth[i][j];
            load[i][j] = MAX_ITEM;
    }
  }

  // printf("The BW here for switch %lld and port %d is %lld\n",i,j,bandwidth_mmf[i][j] );
  for(i=0;i<MAX_ITEM;i++){
    flow_vector[i]=-1;
    //dirty_flow[i]=0;
    path_per_flow[i]=0;
  }

  static int flag = 0;
  if (flag == 0) {
    printf("This is generic cplex_engine.\n");
    flag = 1;
  }

  //------------------------------------------
  // reading from the file and saving informtaio to data structures

  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }


  tot_path = 0;
  init_list();
  init_src_dst_list();


  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {

    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    // no need to route here, one variable for each SD pair
    insertsrc(tot_path, next_s);
    insertdst(tot_path, next_d);

    tot_path += 1;


    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
    //printf("############## %d %d\n", next_s, next_d); 
 }



if(fd) fclose(fd);
printf("\n Done reading trafficfile\n");

//------------------------------------------


while(!all_saturated){
  iterator++;
  all_saturated=1;
  lp_solution=0;
  cc=0;

  // need to either open file traffic again or not close the fd in the bottom!!
  if ((ofd = fopen("xgft_cplex_generic.lp", "w")) == NULL) {
    printf("Can't open xgft_cplex_generic.lp for write.\n");
    exit(0);
  }

  /*STEP 1:  generate objectives and constraints at nodes */

  // STEP 2: generate optimization objective
  fprintf(ofd, "Maximize\nobj: ");
  fprintf(ofd, "x");
  fprintf(ofd, "\nSubject To\n");


//CONSTRAINT STARTS

// Find out all the edges and figure out 

// For all demand d, range is totpath!!

  for (i=0;i<totNode;i++){
    for (j=0;j<MAX_DEGREE;j++){
          if (graph[i][j] != -1)
          {
             fprintf(ofd, "c%lld: ", cc++);
             int first=1;

             for (k=0;k<tot_path;k++)
                 {
                  if(first) {
                    fprintf(ofd, "x%lld_%lld_%d ",k,i,graph[i][j]);
                    first=0;
                  }
                  else fprintf(ofd, "+ x%lld_%lld_%d ",k,i,graph[i][j]);
                  if(k%3==0) fprintf(ofd, "\n");
                 }
                fprintf(ofd, " <= %lld\n",bandwidth_mmf[i][j]);

             }  

      }

  }

  // for each demand the sum of all flows entering the edge and leaving the edge are equal

  //implements both path conservation constraints and first set of traffic constraints
  for (k=0; k<tot_path; k++){
    for (i=0;i<totNode;i++)
     {
  
        int first=1;
        int flaff=0;
        long long int lambda=0;
        // if (i != s/d)

        // for the value of i find the corresponding s and d and match them against k

         long long int ptr=0, ptr2=0;  
         ptr = srchead[i];
         ptr2 = dsthead[i];

         while (ptr != MAX_ITEM){
               if (value[ptr] == k)flaff=1;
               ptr= next[ptr];
         }

         while (ptr2 != MAX_ITEM){

               if (value[ptr2] == k) flaff=2;
               ptr2 = next[ptr2];
	       }
         if (flaff !=0 ){
           if(flow_vector[k]==-1) continue;
           if(flaff==1)  lambda = flow_vector[k];
           else if(flaff==2) lambda = -1*flow_vector[k];
           else printf("Something must be wrong!");
          }


      fprintf(ofd, "c%lld: ", cc++);
      for (j=0;j<MAX_DEGREE;j++)
      {
              if (graph[i][j] != -1)
              {
                  if(first)
                   {fprintf(ofd, "x%lld_%lld_%d ",k,i,graph[i][j]);
                    first=0;}
                  else fprintf(ofd, "+ x%lld_%lld_%d ",k,i,graph[i][j]);
                  if(j%3==0) fprintf(ofd, "\n");

                  fprintf(ofd, " - x%lld_%d_%lld ",k,graph[i][j], i);
               }
        }

        fprintf(ofd, " = %lld\n", lambda); 

    }


  }

  for(i=0;i<totNode;i++){
    if(srchead[i]!=MAX_ITEM|| dsthead[i]!=MAX_ITEM){
      for(k=0;k<tot_path;k++){
        if(flow_vector[k]==-1){
          //int first=1;
          int flaff=0;
          int ptr = srchead[i];
          int ptr2 = dsthead[i];

          while (ptr != MAX_ITEM)
          {

            if (value[ptr] == k)flaff=1;
            ptr= next[ptr];

          }

          while (ptr2 != MAX_ITEM)
          {

            if (value[ptr2] == k) flaff=2;
            ptr2 = next[ptr2];

  	      }

          if(flaff==0) continue;
          else if (flaff ==1)
          {
            fprintf(ofd, "d%lld_%lld: ",k , cc++);  //using different constraint label 'd' instead of 'c'
            fprintf(ofd, " x ");
            for (j=0;j<MAX_DEGREE;j++)
            {
              if (graph[i][j] != -1){
                fprintf(ofd, " - x%lld_%lld_%d ",k,i,graph[i][j]);
                fprintf(ofd, " + x%lld_%d_%lld ",k,graph[i][j], i);
                if(j%3==0) fprintf(ofd, "\n");
              } 
            }
          }

          else if (flaff ==2 )
          {
            fprintf(ofd, "d%lld_%lld: ",k , cc++);  // different constraint label
            fprintf(ofd, " x ");
            for (j=0;j<MAX_DEGREE;j++)
            {
              if (graph[i][j] != -1) { 
                fprintf(ofd, " + x%lld_%lld_%d ",k,i,graph[i][j]);
                fprintf(ofd, " - x%lld_%d_%lld ",k,graph[i][j], i);
                if(j%3==0) fprintf(ofd, "\n");
              }
            }
          }


          fprintf(ofd, " <= 0\n");
        }
      }
    }
  }

  //printf("Printed 4th, 5th set of constraints\n");


  // Last constraint //

  for(i=0;i<totNode;i++ ){
    for(j=0;j<MAX_DEGREE;j++){
      if (graph[i][j] != -1){ 
        for (k = 0 ; k<tot_path;k++){

            fprintf(ofd, "c%lld: x%lld_%lld_%d ",cc++, k,i,graph[i][j]);
            fprintf(ofd, " >= 0\n");     
        }
      }
    }
  }

  //printf("Printed final set of constraints\n");



  fprintf(ofd, "End\n");
  if( ofd )fclose(ofd);   //finish write to LP file



  //STEP 3: Invoke CPLEX and obtain solution objective value
  runcmd("rm b.out");
  runcmd(" cplex -f cplex_command_generic.cmd  > b.out");
  //runcmd(" cplex -f cplex_command_generic.cmd  >> debug.out");
  //runcmd("echo   $\'\\n\' >> debug.out"); //add an extra newline character in dumped solution file

  if ((cplex_fd = fopen("b.out", "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }
  runcmd("echo   $\'\\n\' >> b.out"); //add an extra newline character in dumped solution file

  //runcmd(" cplex -f a.cmd >> d.out");
  long long int demand;
  double var;
  char dummy_buffer[100];
  do{

     if(fgets(dummy_buffer, sizeof(dummy_buffer), cplex_fd)==NULL) printf("error reading from solution file\n");
  }while(dummy_buffer[0]!='x');
  sscanf(dummy_buffer, "x %lf\n",&lp_solution );

  if(fgets(dummy_buffer, 100, cplex_fd)==NULL) printf("error reading from solution file");
  fputs(dummy_buffer, cplex_fd);
  while(fscanf(cplex_fd, "d%lld_%lld %lf\n",&demand,&i, &var  ))
  {
   if(var > 0 && flow_vector[demand]==-1){
     //printf("Flow #%lld has been saturated\n", demand);
     flow_vector[demand] = lp_solution;
     var=0;
   }
 }
  printf("Solution at the end of step %d: %lf\n", iterator, lp_solution);
  if( cplex_fd )fclose(cplex_fd);
  



  //Step 8: continue iteration if all variables are not set
  for (i=0;i<tot_path;i++){
    if(flow_vector[i]==-1) all_saturated=0;
  }




}  //end while(all_saturated)


 long long int tot_mmf_bw=0;
  strcpy(outputfilename, filename);
  strcat(outputfilename,".outGEN");
// strcpy(outputfilename, "OUT_GEN_");
// strcat(outputfilename, filename);
 if ((fd = fopen(outputfilename, "w")) != NULL) {
    for(i=0;i<tot_path;i++){
	     fprintf(fd, "%lld\n", flow_vector[i]);
       tot_mmf_bw+=flow_vector[i];
    }
    fprintf(fd, "Total BW: %lld\n", tot_mmf_bw);
    fclose(fd);
  }
  gettimeofday(&t2, NULL);
  *exec_time=timediff(t1,t2);
  *iteration_count=iterator;
  return tot_mmf_bw;
}



// Algorithm LP 
long long int xgft_mmf_cplex_from_trace_file(char *filename, int* iteration_count, double *exec_time){

  struct timeval t1,t2,t3,t4,t5,t6;
  FILE *fd, *ofd;
  char buf[1000], *ch;
  long long int i, j, k;
  long long int tot_path;
  long long int ii, jj, kk;
  long long int counter;
  int next_s, next_d;
  long long int cc=0, dd=0;

  /*Variables/storage specific to the iterative algorithm*/
  FILE *cplex_fd;			//file descriptor to read in LP solution
  char cmd[1000];			//command string to invoke CPLEX from shell
  double lp_solution=0;			// solution objective value
  int iterator=0;
  int all_saturated=0;			//boundary condition flag

  long long int tot_mmf_bw=0;

 // %% Step 1: initilizing bandwidth input to LP formulation
  gettimeofday(&t1, NULL);

  //initilizing bandwidth input to LP formulation
  for(i=0;i<MAX_NODE;i++){
      for(j=0;j<MAX_DEGREE;j++){
            bandwidth_mmf[i][j] = bandwidth[i][j];
            load[i][j] = MAX_ITEM;
    }
  }


  for(i=0;i<MAX_ITEM;i++){
    flow_vector[i]=-1;
    dirty_flow[i]=0;
  }

  static int flag = 0;
  if (flag == 0) {
    printf("This is xgft cplex_engine.\n");
    flag = 1;
  }
  
  // reading from the file and saving informtaion to data structures
  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }


  tot_path = 0;
  init_list();
  init_src_dst_list();

  // %% Step 2: Read from file, insert into local data structure
  gettimeofday(&t2, NULL);
  printf("Time to initialize: %lf\n", timediff(t1,t2));

   
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    // no need to route here, one variable for each SD pair
    insertsrc(tot_path, next_s);
    insertsrc(tot_path+1, next_s);
    insertdst(tot_path, next_d);
    insertdst(tot_path+1, next_d);

    tot_path += 1;


    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }


  for (i=0; i<totPE; i++) {
     int rule=0;
     long long int ptr1;
     ptr1 = srchead[i];
       while (ptr1 != MAX_ITEM)
          {
            rule++;
            ptr1 = next[next[ptr1]];
           }
	 }
  for (i=0; i<totPE; i++) {
     int rule=0;
     long long int ptr1;
     ptr1 = dsthead[i];
       while (ptr1 != MAX_ITEM)
          {
            rule++;
            ptr1 = next[next[ptr1]];
           }
	 }
 // tot_path is updated 

// Next up reading all the values from file again to read SW info 

   init_sw_list();
   tot_path = 0;
   fseek(fd, 0L, SEEK_SET);


  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    j = xgft_allpath_routing(next_s, next_d, sdpaths);
        //printf("%lld path\n", j);
        //print_allpath(next_s, next_d, sdpaths, j);
        //printf("--------------------------------\n");


   for (i=0; i<j; i++) {
      for (k=1; (sdpaths[i][k+1] != next_d) && (sdpaths[i][k+1] != -1); k++) {
        int port;
        for (port = 0; (port < MAX_DEGREE) && 
         (graph[sdpaths[i][k]][port] != sdpaths[i][k+1]); port ++);
        if (port >= MAX_DEGREE) {
          printf("Can't find node %d from %d\n", sdpaths[i][k+1],sdpaths[i][k]);
          exit(0);
        }
        insertsw((long long int)(tot_path), sdpaths[i][k] - totPE, port);
        insertsw((long long int)(j), sdpaths[i][k] - totPE, port);
      }
    }

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }

  if(fd) fclose(fd);



  //------------------------------------------
  //Step %% 3: begin iterative algorithm
  gettimeofday(&t3, NULL);
  printf("Time to read from input: %lf\n", timediff(t2,t3));
  printf("Paths processed during input %lld\n", tot_path);


while(!all_saturated){
  iterator++;
  lp_solution=0;
  // need to either open file traffic again or not close the fd in the bottom!!


  all_saturated=1;
  if ((ofd = fopen("xgft_cplex_gen0000.lp", "w")) == NULL) {
    printf("Can't open xgft_cplex_gen0000.lp for write.\n");
    exit(0);
  }

  /*STEP 1:  generate objectives and constraints at nodes */

  // STEP 2: generate optimization objective
  fprintf(ofd, "Maximize\nobj: ");
  fprintf(ofd, "x");
  fprintf(ofd, "\nSubject To\n");

  // source constraint
  for (i=0; i<totPE; i++) {
    if (srchead[i] != MAX_ITEM) { // constraints
      long long int ptr;
      int counter = 0;
        fprintf(ofd, "c%lld: ", cc++);
      ptr = srchead[i];

      while (ptr != MAX_ITEM) {
        

        if(flow_vector[value[ptr]]==-1)
          counter++;
        // only when the link is not set/saturated it is being considered to saturate another link/ 
        ptr = next[next[ptr]];
      }
      load[i][0]=counter;
      if(counter>0){
        //printf("I am the %lld th PE \n",i);							//to speed-up future lookup
        fprintf(ofd, "%d x", counter);
        fprintf(ofd, " <= %lld\n",bandwidth_mmf[i][0]);
      }
      else srchead[i]=MAX_ITEM;
    }
  }

   //destination constraint
   for (i=0; i<totPE; i++) {
    if (dsthead[i] != MAX_ITEM) { // constraints
      long long int ptr;
      int counter = 0;

      fprintf(ofd, "c%lld: ", cc++);
      ptr = dsthead[i];
      while (ptr != MAX_ITEM) {
        
        if(flow_vector[value[ptr]]==-1)counter++;
        // only when the link is not set/saturated it is being considered to saturate another link// 
        ptr = next[next[ptr]];
      }

    int min=0,max=0;int flaff=0;
    for(min=0; min <MAX_NODE; min++)
        {    for (max=0;max<MAX_DEGREE;max++)
            {
               if (graph[min][max]== i) {
                //printf("The desired destination node %lld is connected to port %d of switch %d\n", i, max, min);
                flaff=1;
		break;
               }

            }
         if (flaff==1)
              break;

        }
       load[min][max]=counter;
      if(counter>0){							//to speed-up future lookup
        fprintf(ofd, "%d x", counter);
        fprintf(ofd, " <= %lld\n",bandwidth_mmf[min][max]);
      }
      else dsthead[i]=MAX_ITEM;
    }
  }


  // compute and generate constraints for each link between switches
  //init_list();								// value list NOT CLEARED after writing S/D constraints


  // generate the constraints
  // destination constraint
  for (i=0; i<totNode - totPE; i++) {
    for (j=0; j<MAX_DEGREE; j++){

     

      if (swhead[i][j] != MAX_ITEM) { // constraints
        long long int ptr;
        int first=1;
        int counter = 0;

        fprintf(ofd, "c%lld: ", cc++);
        ptr = swhead[i][j];
        while (ptr != MAX_ITEM) {
          long long int v;
          double f;
          v= value[ptr];
          ptr = next[ptr];
          f = 1.0 / value[ptr];
          if(flow_vector[v]==-1){
            if (first) {fprintf(ofd, "%.12lf x", f); first = 0;}
            else fprintf(ofd, " + %.12lf x", f);
            if (counter % 3 == 0) fprintf(ofd, "\n");
            counter++;
          }
          ptr = next[ptr];
        }
        load[i+totPE][j] = counter;									//for future lookup

        if(counter>0){              //to speed-up future lookup
          fprintf(ofd, " <= %lld\n", bandwidth_mmf[i+totPE][j] );  //all links are the same for now    
         // iprintf("The BW here for switch %lld and port %d is %lld\n",i,j,bandwidth_mmf[i][j] );       }
           }
         else swhead[i][j]=MAX_ITEM;
      }
    }
  }

  fprintf(ofd, "End\n");
  if( ofd )fclose(ofd);
//  if ((ofd = fopen("xgft_cplex_gen0000.lp", "w")) == NULL) {

  //STEP 3: Invoke CPLEX and obtain solution objective value
  runcmd("rm -f b.out");
  runcmd(" cplex -f cplex.cmd|tail -n2 >b.out");
  if ((cplex_fd = fopen("b.out", "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }
  runcmd("echo   $\'\\n\' >> b.out"); 
  ch = fgets(buf, 1000, cplex_fd);
  i = sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s  %lf", &lp_solution);
  printf("Solution at the end of step %d: %lf\n", iterator, lp_solution);
  fflush(0);
  if( cplex_fd )fclose(cplex_fd);




  //STEP 4: Find PE link(s) saturated by the value of lp_solution

  //4.1: Source links (PEs)
  long long int ptr, ptr2;
  for (i=0; i<totPE; i++) {
    if (srchead[i] != MAX_ITEM) { // constraints
      ptr = srchead[i];
      double diff=0.0;
      diff = bandwidth_mmf[i][0] - load[i][0]*lp_solution;
           if(diff < BW_SATURATION_GAP ){		//if the link is saturated, update all flows going through it
        while (ptr != MAX_ITEM) {
          if(flow_vector[value[ptr]]==-1) {
              flow_vector[value[ptr]] =lp_solution;
              //DEBUG
              tot_mmf_bw+=lp_solution;

              dirty_flow[value[ptr]] = 1;
          }
          ptr = next[next[ptr]];
        }

        //printf("Link from source %lld to switch has been saturated\n", i);
        load[i][0] =0;                  //for future lookup
        srchead[i] = MAX_ITEM;								//saturated
      }
    }
  }

  //printf("AFTER 4.1: \n");
  //dumpvector(tot_path);

  //4.1.1: Destiation links (PEs)
  for (i=0; i<totPE; i++) {
    if (dsthead[i] != MAX_ITEM) { // constraints
      ptr = dsthead[i];
      double diff=0.0;
      int flaff=0;
      // find j, k for the corresponding i//
    for (j=totPE;j<totNode;j++)
     { for (k=0;k<MAX_DEGREE;k++)
       {
             if (graph[j][k] == i){
              flaff=1; 
              break;
             }

        }
        if (flaff==1 )
         break;
      }

      diff = bandwidth_mmf[j][k] - load[j][k]*lp_solution;
      if(diff < BW_SATURATION_GAP){		//if the link is saturated, update all flows going through it
        while (ptr != MAX_ITEM) {
          if(flow_vector[value[ptr]]==-1){
             flow_vector[value[ptr]] =lp_solution;
              //DEBUG
              tot_mmf_bw+=lp_solution;
             
             //printf ("path/flow %lld goes to destination link %lld\n",value[ptr],i);
             dirty_flow[value[ptr]] = 1;
          }
          ptr = next[next[ptr]];
        }

        //printf("Link to destination %lld from switch has been saturated\n", i);
        load[j][k] =0;                  //for future lookup
       // bandwidth[j][k]=MAX_ITEM;
        dsthead[i] = MAX_ITEM;								//saturated
      }
    }
  }

  //4.2: Switchports
  for (i=0; i<totNode - totPE; i++) {
    for (j=0; j<MAX_DEGREE; j++){
      if (swhead[i][j] != MAX_ITEM) { // constraints
        ptr = swhead[i][j];
         double ratio=0;
         double diff=0;
         //diff = bandwidth_mmf[i+totPE][j] - (1/load[i+totPE][j])*lp_solution;
        while (ptr != MAX_ITEM) {
                         long long int v;
                         double f;
                         v= value[ptr];
                         ptr = next[ptr];
                         f=1.0/value[ptr];
                         if (flow_vector[v]== -1) 
                          ratio += f;
                                 
                
                         ptr = next[ptr];
              }
        // Now if the new Bandwidth is less than GAP, we have a saturation and things need to change
        diff = bandwidth_mmf[i+totPE][j] - ratio*lp_solution;  
        ptr = swhead[i][j];
        
        if( diff < BW_SATURATION_GAP )
           {		//if the link is saturated, update all flows going through it
                bandwidth_mmf[i+totPE][j] -= ratio*lp_solution;
                while (ptr != MAX_ITEM) 
                {
                   if(flow_vector[value[ptr]]==-1){
                      
                      dirty_flow[value[ptr]] = 1;
                      flow_vector[value[ptr]] =lp_solution;
              //DEBUG
              tot_mmf_bw+=lp_solution;

                      //printf("Flow # %lld is tagged dirty\n", value[ptr]);
                    }
                   
	           
                   ptr = next[next[ptr]];
                }
 
           swhead[i][j] = MAX_ITEM;
           //printf("Link from switch %lld to port %lld has been saturated\n", i+totPE, j);
           }
                

       
      }
    }
  }

  //STEP 5: update link bandwidths containing flows saturated/partially saturated in step 4
  // First for PEs
  // Calculate for each link whether a recently saturated 
 //-----
    for (i=0; i<totPE; i++) {
        if (srchead[i] != MAX_ITEM) { // constraints
        long long int ptr;
        ptr = srchead[i];

      while (ptr != MAX_ITEM) {
        if(dirty_flow[value[ptr]]==1)
         {
          load[i][0]--;
          bandwidth_mmf[i][0]-=lp_solution;
        }
            ptr = next[next[ptr]];
       }

       if (load[i][0] == 0)

          srchead[i] = MAX_ITEM; // saturated 


    }

  }

//printf("AFTER 5: \n"); 
//dumpvector(tot_path);

//printf("AFTER 6: ");
//STEP 6: upda19te Switchports bandwidths containing flows saturated/partially saturated in step 4
  //  for outgoing links from all switchports   

for (i=0; i<totNode - totPE; i++) {
    for (j=0; j<MAX_DEGREE; j++){
      if (swhead[i][j] != MAX_ITEM) { // constraints
        ptr = swhead[i][j];
        //printf("[%lld][%lld]: ", i+totPE,j);
        while (ptr != MAX_ITEM) {
	  //printf(" ->%d", value[ptr]);
          if ((dirty_flow[value[ptr]]) == 1){
              //load[i+totPE][j] =0;
              ptr = next[ptr];
              bandwidth_mmf[i+totPE][j] -=  (1.0/value[ptr])*lp_solution;
              //printf("[%d]",bandwidth_mmf[i+totPE][j]);
          }
          else ptr = next[ptr];
          ptr = next[ptr];
        }
        //printf("\n");


  }
 }
}
//STep ??: Destination

    for (i=0; i<totPE; i++) {
        if (dsthead[i] != MAX_ITEM) { // constraints
        long long int ptr;
        ptr = dsthead[i];

      int flaff=0;
      // find j, k for the corresponding i//
      for (j=totPE;j<totNode;j++)
     { for (k=0;k<MAX_DEGREE;k++)
       {
             if (graph[j][k] == i){
              flaff=1;
              break;
             }

        }
        if (flaff==1 )
         break;
      }
      while (ptr != MAX_ITEM) {
        if(dirty_flow[value[ptr]]==1)
         {
           load[j][k]--;
           bandwidth_mmf[j][k]-=lp_solution;
           //printf("Updating BW at [%d][%d] connected to destination %d to %lld\n", j,k,i, bandwidth_mmf[j][k]);
         }
            ptr = next[next[ptr]];
       }

       if (load[j][k] == 0)

          dsthead[i] = MAX_ITEM; // saturated 

    }

  }

	//DEBUG
	/*  for(i=0;i<totNode;i++){
     		 if(i<totPE)
			  printf("The renewed of bandwidths of switch %lld and port %lld are %lld\n",i,j, bandwidth_mmf[i][0]);
    	         else
  	    		  for(j=0;j<6;j++)
         		 	printf("The renewed of bandwidths of switch %lld and port %lld are %lld\n",i,j, bandwidth_mmf[i][j]);

  	}
	*/

// Step 7
// clearing the dirty values

  for (i=0;i<MAX_ITEM;i++){
   dirty_flow[i] = 0;
  }

//Step 8: continue iteration if all variables are not set
  for (i=0;i<tot_path;i++){
    if(flow_vector[i]==-1) all_saturated=0;
  }


  //-------
  printf("tot_mmf_bw at the end of iteration %d %lld ",iterator, tot_mmf_bw );

}

//Step %% 4: end iterative algorithm  
  gettimeofday(&t4, NULL);
  printf("Time to run iterative algo: %lf\n", timediff(t3,t4));

  tot_mmf_bw=0;
  for(i=0;i<tot_path;i++){
     //fprintf(fd, "%lld\n", flow_vector[i]);
     tot_mmf_bw+=flow_vector[i];
    }

  //end while(all_saturated)
  //strcpy(outputfilename, "OUT_LP_");
  //strcpy(outputfilename, filename);
  //strcat(outputfilename,".outLP");
  //printf("Output filename: %s", outputfilename);
  //if ((fd = fopen(outputfilename, "w")) != NULL) {
    //fprintf(fd, "Total XGFT_VLB BW: %lld\n", tot_mmf_bw);
    //fclose(fd);
  //}
  //else printf("Error dumping flow vector!\n");
  gettimeofday(&t5, NULL);           
  printf("Time to calculate final mmf_bw: %lf\n", timediff(t4,t5));
  printf("total elapsed time at lp_fn: %lf\n", timediff(t1,t5));

  *exec_time = timediff(t2,t5);
  *iteration_count=iterator;
  return tot_mmf_bw;
}





/* NON LP VERSION*/
long long int xgft_mmf_nonlp_from_trace_file (char* filename, int* iteration_count, double *exec_time)
{

  struct timeval t1,t2,t3,t4,t5,t6;
  FILE *fd, *ofd;
  char buf[1000], *ch;
  long long int i, j, k;
  long long int tot_path;
  long long int ii, jj, kk;
  long long int counter;
  int next_s, next_d;

  /*Variables/storage specific to the NON-LP iterative algorithm*/
  int iterator=0;
  int all_saturated=0;                  //boundary condition flag
  long long int min_rate_limit_sw;
  double rate_limit_loadfactor, min_rate_limit;


  // %% Step 1: initilizing bandwidth input to LP formulation
  gettimeofday(&t1, NULL);

  //initilizing bandwidth input to LP formulation
  for(i=0;i<MAX_NODE;i++){
      for(j=0;j<MAX_DEGREE;j++){
            bandwidth_mmf[i][j] = bandwidth[i][j];
            load[i][j] = 0;
    }
  }

  for(i=0;i<MAX_ITEM;i++){
    temp_flow_vector[i]=0.0;	//initialized to 0 because  temp_flow_vector i always updated over its previous value
    final_flow_vector[i]=-1;	// -1 indicating unsaturated
  }



  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }

  tot_path = 0;
  init_list();
  init_src_dst_list();
  init_sw_list();

  // %% Step 2: Read from file, insert into local data structure
  gettimeofday(&t2, NULL);
  printf("Time to initialize: %lf\n", timediff(t1,t2));

  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    load[next_s][0]++;			// for links connecting PE to SEs
    insertsrc(tot_path, next_s);	// required for STEP 5

    j = xgft_allpath_routing(next_s, next_d, sdpaths); //j=pathcount in the demand [next_s, next_d]
    for (i=0; i<j; i++) {				//for each path
      for (k=1; sdpaths[i][k+1] != -1; k++) {           //for each hop in selected path

        //find the 'sw' and'port' component of the hop link
        int sw, port;

        sw = sdpaths[i][k]-totPE;
        for (port = 0; (port < MAX_DEGREE) &&
         (graph[sdpaths[i][k]][port] != sdpaths[i][k+1]); port ++);
        if (port >= MAX_DEGREE) {
          printf("Can't find node %d from %d\n", sdpaths[i][k+1],sdpaths[i][k]);
          exit(0);
        }
        else if(!all2allpattern_input && value[swtail[sw][port]]!=tot_path){		//if the current demand is not already inserted, insert it
          load[sdpaths[i][k]][port]++;
          insertsw(tot_path, sw, port);	// required in STEP 5
        }
      }
    }

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }

  if(fd) fclose(fd);

  if(all2allpattern_input){
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){    
        all2all_update_link(filename, i,j,1);	//call function to calculate link load/usage
      }
    }
  }


  //Step %% 3: begin iterative algorithm
  gettimeofday(&t3, NULL);
  printf("Time to read from input: %lf\n", timediff(t2,t3));
  printf("Paths processed during input %lld\n", tot_path);

  while(!all_saturated){
    iterator++;
    all_saturated=1;
    min_rate_limit = LLONG_MAX;		//temporary
    min_rate_limit_sw=-1;
    rate_limit_loadfactor = -1;

    //STEP 1: Find the most rate limiting link(L)
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1 && load[i][j]!=0){
          double  rate_limit;
          double  load_factor;
          if(node_level[i] < node_level[graph[i][j]])	//if link is an uplink
            load_factor = (load[i][j]*1.0)/fan_out(i);
          else						//f link is a downlink
            load_factor = (load[i][j]*1.0)/fan_out(graph[i][j]);
          rate_limit = bandwidth_mmf[i][j] / load_factor;

          //STEP 2: compute the smallest rate that every flow can increase to
          if(rate_limit < min_rate_limit){
            min_rate_limit = rate_limit;
            min_rate_limit_sw = i;		//for debugging
            rate_limit_loadfactor = load_factor;//for debugging
          }
        }
      }
    }

    // STEP 3: For all active SD pair (src, dst),CurrentRate[src][dst] += R;
    int first =1;
    for(i=0;i<tot_path;i++){
      if(final_flow_vector[i]==-1) temp_flow_vector[i]+=min_rate_limit;
      if(first){ 
        printf("At iteration %d, current rate limit: %lf , current saturation flow value: %lf\n",iterator, min_rate_limit, temp_flow_vector[i]); 
        fflush(0);
        first=0;
      }
    }
    

    //STEP 4: Compute the remaining BW on each link after the rate increase
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1 && load[i][j]!=0){
          if(node_level[i] < node_level[graph[i][j]])	//if link is an uplink
            bandwidth_mmf[i][j] -= (long long) (min_rate_limit *load[i][j]/fan_out(i));
          else
            bandwidth_mmf[i][j] -= (long long)(min_rate_limit *load[i][j]/fan_out(graph[i][j]));

        }
      }
    }

    //STEP 5: For all active SD pair (src, dst)
    //  if (a link from in the path src to dst has 0 remaining bandwidth) , update the final rate allocation and remove the SD pair from active list
    //5.1 PE to SE links
    for(i=0;i<totPE;i++){
      if(graph[i][0]!=-1 && bandwidth_mmf[i][0]<=load[i][0]  && srchead[i]!=MAX_ITEM){	//at the most loaded link, leftover bandwidth should be between 0 and load factor
        load[i][0]=0;
        long long int ptr = srchead[i];
        while (ptr != MAX_ITEM) {
          if(final_flow_vector[value[ptr]]==-1){
            final_flow_vector[value[ptr]]=temp_flow_vector[value[ptr]];
          }
          ptr = next[ptr]; 
        }
        srchead[i]=MAX_ITEM;
      }
    }

    //5.2 SWitchports
    for(i=totPE;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(!all2allpattern_input){
          if(graph[i][j]!=-1 && bandwidth_mmf[i][j]<=load[i][j]  && swhead[i-totPE][j]!=MAX_ITEM){
            load[i][j]=0;
            long long int ptr = swhead[i-totPE][j];
            int tmp=0;
            while (ptr != MAX_ITEM) {
              if(final_flow_vector[value[ptr]]==-1){
                final_flow_vector[value[ptr]]=temp_flow_vector[value[ptr]];
                tmp++;
              }
              ptr = next[ptr];
            }
            swhead[i-totPE][j]=MAX_ITEM;
          }
        }else if(graph[i][j]!=-1 && bandwidth_mmf[i][j]<=load[i][j]) {											//all to all patterns
           all2all_update_link(filename, i,j,0);
           if(i%1000==0)printf("All2All: finished scanning %lld switches for saturation", i);
        }

      }
    }

   //STEP 6: update the  link loads of saturated flows
   for(i=0;i<totPE;i++){
     load[i][0]=0;
     long long int ptr = srchead[i];
     while(ptr!=MAX_ITEM){
       if(final_flow_vector[value[ptr]]==-1) load[i][0]++;
       ptr = next[ptr];
    }
   }

  for(i=totPE;i<totNode;i++){
    for(j=0;j<MAX_DEGREE;j++){
      if(graph[i][j]!=-1){
        load[i][j]=0;
        if(!all2allpattern_input){
          long long int ptr = swhead[i-totPE][j];
          while(ptr!=MAX_ITEM){
            if(final_flow_vector[value[ptr]]==-1) load[i][j]++;
            ptr = next[ptr];
          }
        }else all2all_update_link(filename, i,j,1);
      }
    }
    //if(i%1000==0)printf("All2All: finished scanning %lld switches for load calculation", i);

  }

    //evaluate loop boundary condition
    for(i=0;i<tot_path;i++){
      if(final_flow_vector[i]==-1){
         all_saturated=0; 
         break;
      }
    }
  }//end while
  //Step %% 4: end iterative algorithm
  gettimeofday(&t4, NULL);
  printf("Time to run iterative algo: %lf, #iterations = %d\n", timediff(t3,t4), iterator);

  long long int tot_mmf_bw=0;
  for(i=0;i<tot_path;i++){
     //fprintf(fd, "%lld\n", flow_vector[i]);
     tot_mmf_bw+=final_flow_vector[i];
    }

  /*ofd = fopen("nonlpdump","w");
  for(i=0;i<tot_path;i++){
     fprintf(ofd, "Flow[%lld] allocation: %lld\n",i, final_flow_vector[i]);
  }*/

  //Dump result into file
 /* strcpy(outputfilename, filename);
  strcat(outputfilename,".outNONLP");
  if ((ofd = fopen(outputfilename, "w")) != NULL) {
    for(i=0;i<tot_path;i++){
        fprintf(ofd, "%lld\n", final_flow_vector[i]);
        tot_mmf_bw+=final_flow_vector[i];
    }
    fprintf(ofd, "Total NON_LP BW: %lld\n", tot_mmf_bw);
  }
  if(ofd) fclose(ofd);
  printf("NON_LP returned value: %lld\n", tot_mmf_bw);*/

  gettimeofday(&t5, NULL);
  printf("Time to calculate final mmf_bw: %lf\n", timediff(t4,t5));
  printf("total elapsed time at nonlp_fn: %lf\n", timediff(t1,t5));

  *iteration_count = iterator;
  *exec_time = timediff(t2,t5);  



  return tot_mmf_bw;
}




/* NON LP Approximate VERSION*/
long long int xgft_mmf_dmodk_from_trace_file (char* filename, int*  iteration_count, double *exec_time)
{

  struct timeval t1,t2,t3,t4,t5,t6;
  FILE *fd, *ofd;
  char buf[1000], *ch;
  long long int i, j, k;
  long long int tot_path;
  long long int ii, jj, kk;
  long long int counter;
  int next_s, next_d;

  /*Variables/storage specific to the NON-LP iterative algorithm*/
  int iterator=0;
  int all_saturated=0;                  //boundary condition flag
  long long int  min_rate_limit_sw;
  double min_rate_limit;
  int rate_limit_load;


  // %% Step 1: initilizing bandwidth input to LP formulation
  gettimeofday(&t1, NULL);
  for(i=0;i<MAX_NODE;i++){
      for(j=0;j<MAX_DEGREE;j++){
            bandwidth_mmf[i][j] = bandwidth[i][j];
            load[i][j] = 0;
    }
  }

  for(i=0;i<MAX_ITEM;i++){
    temp_flow_vector[i]=0;           //initialized to 0 because  temp_flow_vector i always updated over its previous value
    final_flow_vector[i]=-1;    // -1 indicating unsaturated
  }


  // reading from the file and saving informtaio to data structures
  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }


  tot_path = 0;
  init_list();
  init_src_dst_list();
  init_sw_list();


  // %% Step 2: Read from file, insert into local data structure
  gettimeofday(&t2, NULL);
  printf("Time to initialize: %lf\n", timediff(t1,t2));
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }


    load[next_s][0]++;                  // for links connecting PE to SEs
    insertsrc(tot_path, next_s);        // required for STEP 5
    xgft_dmodk_routing(next_s, next_d, path); //copy d-mod-k path between src and dst into path

    for (k=1; path[k+1] != -1; k++) {
      int sw, port;
      
      sw = path[k]-totPE;
      for (port = 0; (port < MAX_DEGREE) &&
       (graph[path[k]][port] != path[k+1]); port++);
      if (port >= MAX_DEGREE) {
        printf("Can't find node %d from %d\n", path[k+1],path[k]);
        exit(0);
      }
       else if(!all2allpattern_input){            //if the current demand is not already inserted, insert it
        load[path[k]][port]++;
        insertsw(tot_path,sw, port);      // required in STEP 5
      }
    }

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);

  } //end while

  if(fd) fclose(fd);



  if(all2allpattern_input){
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        printf("All2All flag is on\n");
        all2all_update_link(filename, i,j,1);	//call function to calculate link load/usage
      }
    }
  }


  //Step %% 3: begin iterative algorithm
  gettimeofday(&t3, NULL);
  printf("Time to read from input: %lf\n", timediff(t2,t3));
  printf("Paths processed during input %lld\n", tot_path);

  while(!all_saturated){
    gettimeofday(&t4, NULL);
    iterator++;
    all_saturated=1;
    min_rate_limit = LLONG_MAX;         //temporary
    min_rate_limit_sw=-1;
    rate_limit_load=-1;

    //STEP 1: Find the most rate limiting link(L)
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1 && load[i][j]!=0){
          double  rate_limit;
          rate_limit = bandwidth_mmf[i][j]*1.0 / load[i][j];
          //STEP 2: compute the smallest rate that every flow can increase to
          if(rate_limit < min_rate_limit){
            min_rate_limit = rate_limit;
            min_rate_limit_sw = i;
            rate_limit_load = load[i][j];
            //printf("min_rate_limit at [%d %d], current BW: %lld, load %d, rate_limit: %lf\n", i,j,bandwidth_mmf[i][j], load[i][j], min_rate_limit);
            //if(min_rate_limit<1) sleep(1);
          }
        }
      }
    }
    // STEP 3: For all active SD pair (src, dst),CurrentRate[src][dst] += R;
    int first=1;
    for(i=0;i<tot_path;i++){
      if(final_flow_vector[i]==-1){
        temp_flow_vector[i]+=min_rate_limit;
        if (first){
          printf("At iteration %d, current rate limit:%lf\n",iterator, min_rate_limit);
          fflush(0);
          first=0;

        }
      }
    }

    //STEP 4: Compute the remaining BW on each link after the rate increase
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1 && load[i][j]!=0)
          bandwidth_mmf[i][j] -=(long long) min_rate_limit *load[i][j];
      }
    }


    //int saturation_flag=1;
    //STEP 5: For all active SD pair (src, dst)
    //  if (a link from in the path src to dst has 0 remaining bandwidth) , update the final rate allocation and remove the SD pair from active list
    //5.1 PE to SE links
    for(i=0;i<totPE;i++){
      if(graph[i][0]!=-1 && bandwidth_mmf[i][0]<=rate_limit_load  && srchead[i]!=MAX_ITEM){        //at the most loaded link, leftover bandwidth should be between 0 and the corresponding load
        //if(saturation_flag) {printf("link [%lld 0] with BW %lld saturated\n", i, bandwidth_mmf[i][0]); saturation_flag=0;}
        load[i][0]=0;
        long long int ptr = srchead[i];
        while (ptr != MAX_ITEM) {
          if(final_flow_vector[value[ptr]]==-1){
            final_flow_vector[value[ptr]]=temp_flow_vector[value[ptr]];
          }
          ptr = next[ptr];
        }
        srchead[i]=MAX_ITEM;

      }
    }

    //5.2 SWitchports
    for(i=totPE;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(!all2allpattern_input){
          if(graph[i][j]!=-1 && bandwidth_mmf[i][j]<=rate_limit_load  && swhead[i-totPE][j]!=MAX_ITEM){
            //if(saturation_flag) {printf("link [%lld %lld] with BW %lld saturated\n", i,j,bandwidth_mmf[i][0]); saturation_flag=0;}
            load[i][j]=0;
            long long int ptr = swhead[i-totPE][j];
            while (ptr != MAX_ITEM) {
              if(final_flow_vector[value[ptr]]==-1){
                final_flow_vector[value[ptr]]=temp_flow_vector[value[ptr]];
              }
              ptr = next[ptr];
            }
            swhead[i-totPE][j]=MAX_ITEM;
          }
        }else if(graph[i][j]!=-1 && bandwidth_mmf[i][j]<=load[i][j]) {   //all2all pattern
           printf("all to all FLAG IS ##################################################################### %d\n",all2allpattern_input);
           all2all_update_link(filename, i,j,0);		//saturate link [i][j]
           if(i%1000==0)printf("All2All: finished scanning %lld switches for saturation", i);

        }

      }
    }
    /*if(saturation_flag==1) {
      printf("ERROR: At least one link must be saturated at each iteration!\n");
      printf("no saturated link found after %d iterations, limiting BW: %lf, limiting load: %d at sw: %lld \n",
              iterator, min_rate_limit, rate_limit_load, min_rate_limit_sw );
      exit(0);
    }*/

    //STEP 6: update the  link loads of saturated flows
    for(i=0;i<totPE;i++){
      load[i][0]=0;
      long long int ptr = srchead[i];
      while(ptr!=MAX_ITEM){
        if(final_flow_vector[value[ptr]]==-1) load[i][0]++;
        ptr = next[ptr];
      }
    }

    for(i=totPE;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1){
          load[i][j]=0;
          if(!all2allpattern_input){
            long long int ptr = swhead[i-totPE][j];
            while(ptr!=MAX_ITEM){
              if(final_flow_vector[value[ptr]]==-1) load[i][j]++;
              ptr = next[ptr];
            }
          } else all2all_update_link(filename, i,j,1);
        }
      }
      //if(i%1000==0)printf("All2All: finished scanning %lld switches for load calculation", i);
    }

    //evaluate loop boundary condition
    for(i=0;i<tot_path;i++){
      if(final_flow_vector[i]==-1){
         all_saturated=0;
         break;
      }
    }

   gettimeofday(&t5, NULL);
   //printf("time elapsed at iteration %d: %lf\n", iterator, timediff(t4,t5));
  }//end while

  long long int tot_mmf_bw=0;
    for(i=0;i<tot_path;i++){
        //fprintf(ofd, "%lld\n", final_flow_vector[i]);
        tot_mmf_bw+=final_flow_vector[i];
    }


  //Dump result into file
/*
  strcpy(outputfilename, filename);
  strcat(outputfilename,".outDMODK");
  if ((ofd = fopen(outputfilename, "w")) != NULL) {
    fprintf(ofd, "Total DMODK BW: %lld\n", tot_mmf_bw);
  }
  if(ofd) fclose(ofd);
*/

  gettimeofday(&t6, NULL);

  printf("Time to run iterative fn: %lf\n", timediff(t3,t5));
  printf("Time to calculate final mmf_bw: %lf\n", timediff(t5,t6));
  printf("total elapsed time at dmodk_fn: %lf\n", timediff(t1,t6));


  *iteration_count=iterator;
  *exec_time = timediff(t2,t6);

  return tot_mmf_bw;
}



//CPLEX version
double  xbar_from_trace_file(char* filename, int* flow_count ){
  FILE *fd, *ofd, *cplex_fd;
  char buf[1000], *ch;
  double lp_solution;
  long long int i,j,k, cc=0;
  long long int tot_path;
  long long int counter;

  int next_s, next_d;

  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }
  if ((ofd = fopen("crossbar_cplex_gen0000.lp", "w")) == NULL) {
    printf("Can't open crossbar_cplex_gen0000.lp for write.\n");
    exit(0);
  }

  tot_path=0;
  init_list();
  init_src_dst_list();  //init_sw_list() not required

  counter=0;
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    if (counter % 10000 == 0) {
      //printf("counter = %lld\n", counter++); 
      fflush(0);
    }

    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }
    insertsrc(tot_path, next_s);
    insertdst(tot_path, next_d);

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }
  *flow_count = tot_path;
  fprintf(ofd, "Maximize\nobj: ");
  for (i=0; i<tot_path-1; i++) {
    fprintf(ofd, "x%lld + ", i); 
    if (i % 3 == 0) fprintf(ofd, "\n");
  }
  fprintf(ofd, "x%lld\n", i);

  fprintf(ofd, "\nSubject To\n");

  // source constraint
  for (i=0; i<totPE; i++) {
    if (srchead[i] != MAX_ITEM) { // constraints
      long long int ptr;
      int first=1;

      fprintf(ofd, "c%lld: ", cc++);
      ptr = srchead[i];
      while (ptr != MAX_ITEM) {
        if (first) {
            fprintf(ofd, "x%lld", value[ptr]); 
            first = 0;
        } else fprintf(ofd, " + x%lld", value[ptr]);
        ptr = next[ptr];
      }
      fprintf(ofd, "<=%lld\n", bandwidth[1][0]);
    }
  }

  // destination constraint
  for (i=0; i<totPE; i++) {
    if (dsthead[i] != MAX_ITEM) { // constraints
      long long int ptr;
      int first=1;

      fprintf(ofd, "c%lld: ", cc++);
      ptr = srchead[i];
      while (ptr != MAX_ITEM) {
        if (first) {
            fprintf(ofd, "x%lld", value[ptr]); 
            first = 0;
          } else {
            fprintf(ofd, " + x%lld", value[ptr]);
          }
        ptr = next[ptr];
      }
      fprintf(ofd, "<=%lld\n", bandwidth[1][0]);
    }
  }

  fprintf(ofd, "End\n");
  fclose(ofd);

 //STEP 3: Invoke CPLEX and obtain solution objective value

  runcmd("rm -f cplex.out");
  runcmd(" cplex -f cplex_command_crossbar.cmd |tail -n2 > cplex.out");
  runcmd("echo   $\'\\n\' >> cplex.out");

  if ((cplex_fd = fopen("cplex.out", "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }
  ch = fgets(buf, 1000, cplex_fd);
  i = sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s  %le", &lp_solution);
  if( cplex_fd )fclose(cplex_fd);

  return lp_solution;
}



//Max-Min Fair Crossbar throughput calculation without using CPLEX

static int srclist[MAX_ITEM];	// stores source of each demand
static int dstlist[MAX_ITEM];	//stores destination of each demand

double  xbar_mmf_from_trace_file(char* filename, int* flow_count ){
  FILE *fd;
  char buf[1000], *ch;
  double solution=0;
  const double EPSILON=0.1;

  long long int i,j,k,s,d, cc=0;
  long long int tot_path;
  long long int counter;

  int src_load[MAX_NODE]={0};	// stores load(#of SD pairs) at each processing node
  int dst_load[MAX_NODE]={0};
  double bw_src[MAX_NODE];
  double bw_dst[MAX_NODE];
  double max_loaded_rate;



  for(i=0;i<MAX_NODE;i++){
    bw_src[i] = bandwidth[i][0];
    bw_dst[i] = bandwidth[i][0];
  }

  for(i=0;i<MAX_ITEM;i++)
    flow_vector[i]=-1;


  int next_s, next_d;

  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }

  tot_path=0;
  counter=0;

  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    if (counter % 10000 == 0) {
      fflush(0);
    }
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }
    src_load[next_s]++;
    dst_load[next_d]++;
    srclist[tot_path] = next_s;
    dstlist[tot_path] = next_d;

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }
  *flow_count = tot_path;

  long long int saturation_count=0;
  int iteration=0;
  do{
    iteration++;
    //find the smallest rate among all nodes
    //clarification: max_loaded_rate = MINIMUM fair rate;
    max_loaded_rate=LLONG_MAX;
    for(i=0;i<MAX_NODE;i++){
     if(src_load[i]!=0 && max_loaded_rate > (bw_src[i]/src_load[i]))
      max_loaded_rate = bw_src[i]/src_load[i];
     else if(dst_load[i]!=0 && max_loaded_rate > (bw_dst[i]/dst_load[i]))
      max_loaded_rate = bw_dst[i]/dst_load[i];
    }


    for(i=0;i<tot_path;i++){
      if(flow_vector[i]==-1){
        s = srclist[i];
        d = dstlist[i];
        if(fabs(bw_src[s]/src_load[s]-max_loaded_rate) <EPSILON || fabs(bw_dst[d]/dst_load[d] - max_loaded_rate) <EPSILON){
          flow_vector[i]=0;
          saturation_count++;
          solution+=max_loaded_rate;
          bw_src[s]-=max_loaded_rate;
          src_load[s]--;
          bw_dst[d]-=max_loaded_rate;
          dst_load[d]--;
        }
      }
    }
  }while(saturation_count<tot_path);
  return solution;
}



void all2all_update_link(char* filename, int sw,int port, int update_load){

  //read traffic file <filname>. for each traffic in <filename> if it uses the parameter switchport
  //then saturate that flow..as long as it is not already saturated
  FILE *fd;
  char buf[1000], *ch;
  long long int i, j, k;
  int next_s, next_d;
  long long int traffic_counter=0;

  load[sw][port]=0;

  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not existzzzz.\n", filename);
    exit(0);
  }

  //PRIMER
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);

  while (next_s > -1) {
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    j = xgft_allpath_routing(next_s, next_d, sdpaths);
    for (i=0; i<j; i++) {							//for each path
      for (k=1; (sdpaths[i][k+1] != next_d) && (sdpaths[i][k+1] != -1); k++){   //for each hop in each path
        if(sdpaths[i][k]==sw && graph[sdpaths[i][k]][port] == sdpaths[i][k+1]){	//if the hop jumps from [sw,port]
          if(update_load==0){							//0 ==SATURATE THE traffic going through the link
            if(final_flow_vector[traffic_counter]==-1)
              final_flow_vector[traffic_counter]=temp_flow_vector[traffic_counter];
          }
          else if (update_load==1) load[sw][port]++;				//1==UPDATE THE LINK LOAD
        }

      }
    }


    traffic_counter++;
    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }

}



//! This function computes the model index for Jellyfish topology using the traffic
//! from the trace_file
//need to keep track of path count for each SD pair for LLSKR, constant for KSP
static int pathcount[MAX_ITEM] = {_REAL_PATH_NUM}; //same value as jellyfish_K


double generic_cplex_from_trace_file(char* filename){
return -1.0;
}


double rate_limit[MAX_NODE][MAX_DEGREE];

/* NON LP VERSION partial progressive filling*/
long long int xgft_mmf_ppf_from_trace_file (char* filename, int* iteration_count, double *exec_time)
{

  struct timeval t1,t2,t3,t4,t5,t6;
  FILE *fd, *ofd;
  char buf[1000], *ch;
  long long int i, j, k;
  long long int tot_path;
  long long int ii, jj, kk;
  long long int counter;
  int next_s, next_d;
  char c;
  long long int max_concurrent_flow;
  /*Variables/storage specific to the NON-LP iterative algorithm*/
  int iterator=0;
  int all_saturated=0;                  //boundary condition flag
  long long int min_rate_limit_sw;
  double rate_limit_loadfactor, min_rate_limit;

  // %% Step 1: initilizing bandwidth input to LP formulation
  gettimeofday(&t1, NULL);

  //initilizing bandwidth input to LP formulation
  for(i=0;i<MAX_NODE;i++){
      for(j=0;j<MAX_DEGREE;j++){
            bandwidth_mmf[i][j] = bandwidth[i][j];
            load[i][j] = 0;
            rate_limit[i][j]=LLONG_MAX;
    }
  }


  for(i=0;i<MAX_ITEM;i++){
    //temp_flow_vector[i]=0.0;	//initialized to 0 because  temp_flow_vector i always updated over its previous value
    final_flow_vector[i]=-1;	// -1 indicating unsaturated
  }



  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }

  tot_path = 0;
  init_list();
  init_src_dst_list();
  init_sw_list();


  // %% Step 2: Read from file, insert into local data structure
  gettimeofday(&t2, NULL);
  printf("Time to initialize: %lf\n", timediff(t1,t2));

  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  while (next_s > -1) {
    if (next_s == next_d) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }

    load[next_s][0]++;			// for links connecting PE to SEs
    insertsrc(tot_path, next_s);	// required for STEP 5

    j = xgft_allpath_routing(next_s, next_d, sdpaths); //j=pathcount in the demand [next_s, next_d]
    for (i=0; i<j; i++) {				//for each path
      for (k=1; sdpaths[i][k+1] != -1; k++) {           //for each hop in selected path

        //find the 'sw' and'port' component of the hop link
        int sw, port;

        sw = sdpaths[i][k]-totPE;
        for (port = 0; (port < MAX_DEGREE) &&
         (graph[sdpaths[i][k]][port] != sdpaths[i][k+1]); port ++);
        if (port >= MAX_DEGREE) {
          printf("Can't find node %d from %d\n", sdpaths[i][k+1],sdpaths[i][k]);
          exit(0);
        }
        else if(!all2allpattern_input && value[swtail[sw][port]]!=tot_path){		//if the current demand is not already inserted, insert it
          load[sdpaths[i][k]][port]++;
          insertsw(tot_path, sw, port);	// required in STEP 5
        }
      }
    }

    tot_path += 1;

    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
  }

  printf("scan complete\n");


  if(fd) fclose(fd);

  if(all2allpattern_input){
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){    
        all2all_update_link(filename, i,j,1);	//call function to calculate link load/usage
      }
    }
  }


  //Step %% 3: begin iterative algorithm
  gettimeofday(&t3, NULL);
  printf("Time to read from input: %lf\n", timediff(t2,t3));
  printf("Paths processed during input %lld\n", tot_path);
  counter=0;

  while(!all_saturated){
    iterator++;
    all_saturated=1;
    min_rate_limit = LLONG_MAX;		//temporary
    min_rate_limit_sw=-1;
    rate_limit_loadfactor = -1;

    //STEP 1: Find the most rate limiting link(L)
    for(i=0;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(graph[i][j]!=-1 && load[i][j]!=0){
//          double  rate_limit;
          double  load_factor;
          if(node_level[i] < node_level[graph[i][j]])	//if link is an uplink
            load_factor = (load[i][j]*1.0)/fan_out(i);
          else						//f link is a downlink
            load_factor = (load[i][j]*1.0)/fan_out(graph[i][j]);
          
          rate_limit[i][j] = bandwidth_mmf[i][j] / load_factor;

          //STEP 2: compute the smallest rate that every flow can increase to
          if(rate_limit[i][j] < min_rate_limit){
            min_rate_limit = rate_limit[i][j];
            min_rate_limit_sw = i;		//for debugging
            rate_limit_loadfactor = load_factor;//for debugging
          }
        }
      }
    }

    printf("At iteration %d, most limiting rate %lf\n", iterator, min_rate_limit);


    //STEP 5: For all active SD pair (src, dst)
    //  if (a link from in the path src to dst has 0 remaining bandwidth) , update the final rate allocation and remove the SD pair from active list
    //5.1 PE to SE links
    for(i=0;i<totPE;i++){
      
      if(graph[i][0]!=-1 && fabs(rate_limit[i][0]-min_rate_limit)<0.0001  && srchead[i]!=MAX_ITEM){	//at the most loaded link, leftover bandwidth should be between 0 and load factor
        load[i][0]=0;
        long long int ptr = srchead[i];
        while (ptr != MAX_ITEM) {
          if(final_flow_vector[value[ptr]]==-1){
            final_flow_vector[value[ptr]]=min_rate_limit;
            counter++;
          }
          ptr = next[ptr]; 
        }
        srchead[i]=MAX_ITEM;
      }
    }
    printf("Number of flows leaf level flows saturated in itern %d is %lld\n",iterator,counter);
    //counter=0;    
    //5.2 SWitchports
    for(i=totPE;i<totNode;i++){
      for(j=0;j<MAX_DEGREE;j++){
        if(!all2allpattern_input){
          if(graph[i][j]!=-1 && fabs(rate_limit[i][j]-min_rate_limit)<0.0001  && swhead[i-totPE][j]!=MAX_ITEM){
            load[i][j]=0;
            long long int ptr = swhead[i-totPE][j];
            int tmp=0;
            while (ptr != MAX_ITEM) {
              if(final_flow_vector[value[ptr]]==-1){
                final_flow_vector[value[ptr]]=min_rate_limit; counter++;
                tmp++;
              }
              ptr = next[ptr];
            }
            swhead[i-totPE][j]=MAX_ITEM;
          }
        }else if(graph[i][j]!=-1 && bandwidth_mmf[i][j]<=load[i][j]) {											//all to all patterns
           all2all_update_link(filename, i,j,0);
           if(i%1000==0)printf("All2All: finished scanning %lld switches for saturation", i);
        }

      }
    }
      printf("Number of flows swport level flows saturated in itern %d is %lld\n",iterator,counter);


   //STEP 6: update the  link loads of saturated flows
   //STEP 6.1: Update the available bw at each link
    long long int used_bw;

   for(i=0;i<totPE;i++){
     load[i][0]=0;
     used_bw=0;
     long long int ptr = srchead[i];
     while(ptr!=MAX_ITEM){
       if(final_flow_vector[value[ptr]]==-1) load[i][0]++;
       else  used_bw+= final_flow_vector[value[ptr]];
       ptr = next[ptr];
    }

    bandwidth_mmf[i][0]= bandwidth[i][0]-used_bw;
   

   }

  for(i=totPE;i<totNode;i++){
    for(j=0;j<MAX_DEGREE;j++){
      if(graph[i][j]!=-1){
        load[i][j]=0;
        used_bw=0;
        if(!all2allpattern_input){
          long long int ptr = swhead[i-totPE][j];
          while(ptr!=MAX_ITEM){
            if(final_flow_vector[value[ptr]]==-1) load[i][j]++;
            else used_bw+=final_flow_vector[value[ptr]];

            ptr = next[ptr];
          }
        }else all2all_update_link(filename, i,j,1);
        if(node_level[i] < node_level[graph[i][j]])     
          bandwidth_mmf[i][j]= bandwidth[i][j]-(used_bw*1.0/fan_out(i));  
        else bandwidth_mmf[i][j]= bandwidth[i][j]-(used_bw*1.0/fan_out(graph[i][j]));  

      }
    }
    //if(i%1000==0)printf("All2All: finished scanning %lld switches for load calculation", i);

  }

    //evaluate loop boundary condition
    for(i=0;i<tot_path;i++){
      if(final_flow_vector[i]==-1){
         all_saturated=0; 
         break;
      }
    }
    
    if (iterator ==1 && MCF_FLAG==1) {
      max_concurrent_flow = min_rate_limit;
      break;					//early termination to return concurrent_flow instead of mmf flow
    }
  }//end while
  //Step %% 4: end iterative algorithm
  gettimeofday(&t4, NULL);
  printf("Time to run iterative algo: %lf, #iterations = %d\n", timediff(t3,t4), iterator);
 
  long long int tot_mmf_bw=0;
  if(MCF_FLAG){
    tot_mmf_bw = max_concurrent_flow*counter;
  }
  else{
    for(i=0;i<tot_path;i++){
     //fprintf(fd, "%lld\n", flow_vector[i]);
     tot_mmf_bw+=final_flow_vector[i];
    }
  }

  ofd = fopen("nonlpdump","w");
  for(i=0;i<tot_path;i++){
     fprintf(ofd, "Flow[%lld] allocation: %lld\n",i, final_flow_vector[i]);
  }

  //Dump result into file
 /* strcpy(outputfilename, filename);
  strcat(outputfilename,".outNONLP");
  if ((ofd = fopen(outputfilename, "w")) != NULL) {
    for(i=0;i<tot_path;i++){
        fprintf(ofd, "%lld\n", final_flow_vector[i]);
        tot_mmf_bw+=final_flow_vector[i];
    }
    fprintf(ofd, "Total NON_LP BW: %lld\n", tot_mmf_bw);
  }
  if(ofd) fclose(ofd);
  printf("NON_LP returned value: %lld\n", tot_mmf_bw);*/

  gettimeofday(&t5, NULL);
  printf("Time to calculate final mmf_bw: %lf\n", timediff(t4,t5));
  printf("total elapsed time at nonlp_fn: %lf\n", timediff(t1,t5));

  *iteration_count = iterator;
  *exec_time = timediff(t2,t5);  



  return tot_mmf_bw;
}

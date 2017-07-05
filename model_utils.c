#include "topology.h"

#ifdef _SMALL
#define MAX_ITEM 40000LL
#else
#define MAX_ITEM 400000000LL
#endif


long long int flow_vector[MAX_ITEM];

long long int value[MAX_ITEM];              
long long int next[MAX_ITEM];
long long int head = 0;


/** utility function
 *  returns elapsed time of the given interval in seconds
*/
double timediff(struct timeval start, struct timeval end)
{
 return (double)(((end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec))/1000000.0);
}

/**
  * Utility function         
  * @param cmd any shell command string
  * passes unix command to system() for execution
  */
void runcmd(char *cmd)   
{
  int ret;
  fflush(0);
  ret = system(cmd);      
  fflush(0);
  if ((ret == -1)) {      
    printf("Cmd: '%s' failed. %d, %d\n", cmd, ret, ret);
    exit(0);              
  }
}

/**
  * Linked list routine
  * init_list(): initialized the global arrays: value[] and next[]
  * Array implementation of link-list, same arrays(value and next) used to store multiple list
  */
void init_list(){
  long long int i;  
  head = 0;                  
  for (i=0; i<MAX_ITEM-1; i++) {
    next[i] = i+1;
    value[i] = MAX_ITEM;
  }
  next[MAX_ITEM-1] = MAX_ITEM;
  value[MAX_ITEM-1] = MAX_ITEM;
}  

/**
  * Linked list routine
  * LEGACY..
  */
void listfree(long long int index)
{
  next[index] = head;
  value[index] = MAX_ITEM;
  head = index;
}



/**
  * Linked list routine
  * reads up the global variable head (used like a 'free' pointer) and returns its value
  * @return index to the next available/usable data slot(lowest index from available slots)
  */
 
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


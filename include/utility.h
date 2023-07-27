#ifndef UTILITY_  

#define UTILITY_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h> 
#include <sys/time.h>
#include "tsp.h"

#define EPS 1e-5


double eucl_dist(instance *inst, int start_p , int final_p);
void print_error(const char *err);
double rand01();
double get_elapsed_time(struct timeval start, struct timeval end);
void print_sol(instance *inst , int start_n);
int xpos(int i, int j, instance *inst);
void copy_instance(instance *dst, instance *src);
int rand_choice(int from, int to);
double cut_violation( double *xstar , int *index , double *value , char sense , double rhs , int nnz );
void free_instance(instance *inst);
void print_instance(instance *inst);

#endif 
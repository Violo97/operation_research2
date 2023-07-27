#include "tsp.h"
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

double eucl_dist(instance *inst , int start_p , int final_p){
    double dx = inst->nodes[start_p].x - inst->nodes[final_p].x;
	double dy = inst->nodes[start_p].y - inst->nodes[final_p].y;
    if ( !inst->integer_cost ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; 					// nearest integer 
	return dis+0.0;
}

void free_instance(instance *inst)
{     
	free(inst->nodes);
    free(inst->solution.edges);
    //free(inst->solution.best_sol);
}  

void print_instance(instance *inst){
	printf("nnodes: %d\n" , inst->nnodes);
	for(int i = 0 ; i < inst->nnodes ; i++){
		printf("node %d : x->%f y->%f\n\n" , i+1 , inst->nodes[i].x , inst->nodes[i].y);
	}
}

double rand01(){
    return ((double) rand() / RAND_MAX);
}

void print_error(const char *err) { 
    printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); 
}   

double get_elapsed_time(struct timeval start, struct timeval end) {
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    return elapsed;
}

void print_sol(instance *inst , int start_n){

    for(int i = 0 ; i < inst->nnodes ; i++) printf("Start node : %d -> final node : %d \n" , inst->solution.edges[i].i , inst->solution.edges[i].j);

}

int xpos(int i, int j, instance *inst)      // to be verified                                           
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}

void copy_instance(instance *dst, instance *src) {
    memcpy(dst, src, sizeof(instance));
    if (src->nodes) {
        dst->nodes = calloc(src->nnodes, sizeof(point));
        memcpy(dst->nodes, src->nodes, sizeof(point) * src->nnodes);
    }
    if (src->solution.edges) {
        dst->solution.edges = calloc(src->nnodes, sizeof(edge));
        memcpy(dst->solution.edges, src->solution.edges, sizeof(edge) * src->nnodes);
    }
    dst->solution.obj_best = src->solution.obj_best;
}  


int rand_choice(int from, int to) {
    return from + ((int) ((((double) random()) / RAND_MAX) * (to - from)));
}

double cut_violation( double *xstar , int *index , double *value , char sense , double rhs , int nnz ){

    double lhs = 0.0;

    for(int i = 0 ; i < nnz ; i++){

        lhs += xstar[index[i]] * value[i];

    }

    if( sense == 'L' ){

        return lhs - rhs;

    }else if( sense == 'G' ){

        return rhs - lhs;

    }else{

        return fabs(lhs - rhs);
    
    }

}


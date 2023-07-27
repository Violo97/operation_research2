#ifndef TSP_  

#define TSP_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h> 

// Definition of Point
typedef struct {
    double x;
    double y;
} point;

// Edge that connects node i and node j.
typedef struct {
    int i; 
    int j;
} edge;

typedef struct {
	double obj_best;            // Stores the best value of the objective function
	double incumbent_val;
    edge *edges;            // List the solution's edges: list of pairs (i,j)
    double solve_time;   // Time used to solve the istance
    //int *best_sol;          // The best solution found in heuristics implementations
} solution;

//data structures  

typedef struct {   
	
	//input data
	int nnodes; 	   
	point *nodes;

	// parameters 
	int randomseed;
	int timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	int algo;                             //type of algorithm 
	solution solution;
	long ncols;
	int verbose;        // Verbose level of debugging printing : 0 - debug ; 10 - print all the incumbent sol ; 20 - print all the new better sol ; 30 - print the final sol ; >=100 - cplex log
    int active2opt;  // Used in incubement callbacks for 2opt refinement
	char algorithm[20]; 		//algorithm name
    int use_cplex;		
	double prob;  		//probability 
	int integer_cost;
	int preprocess_time; //define the time in which apply the preprocess algo for find an initial sol 
	char* preprocess_algo;//algorithm with which found a initial sol in any applied method 
	
} instance;



#endif 
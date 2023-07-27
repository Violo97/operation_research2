#include <cplex.h>
#include <concorde.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "tsp.h"
#include "benders.h"


static int CPXPUBLIC lazy_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, instance *inst ) {
    int ncols = inst->ncols;
	
    //initialyze the time of execution
	struct timeval start, end;
    gettimeofday(&start, 0);

	pthread_mutex_init(&mutex, NULL);

	double* xstar = (double*) malloc(ncols * sizeof(double));  
	double objval = CPX_INFBOUND; 
    //printf("numbers of cols : %d \n" , ncols);
	if ( CPXcallbackgetcandidatepoint(context, xstar, 0, ncols-1, &objval) ) print_error("CPXcallbackgetcandidatepoint error");
	
	// get some random information at the node (as an example for the students)
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread); 
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode); 
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent); 
	//if ( VERBOSE >= 100 ) 
    printf(" callback at node %5d thread %2d incumbent %10.2lf\n", mynode , mythread , incumbent );
	
	//int nnz = 0; 
	// if xstart is infeasible, find a violated cut and store it in the usual Cplex data structute (rhs, sense, nnz, index and value)
    int *succ = calloc(inst->nnodes , sizeof(int));
	int *comp = calloc(inst->nnodes , sizeof(int));
	int ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
    build_sol(xstar, inst,succ, comp, &ncomp);

	double *xheu = (double *) calloc(inst->ncols, sizeof(double));

	if ( ncomp > 1 ) // means that the solution is infeasible and a violated cut has been found
	{	
		//pthread_mutex_lock(&mutex);
		int* test_succ = calloc(inst->nnodes , sizeof(int));
		memcpy(test_succ, succ, inst->nnodes * sizeof(int));
		heur_patch(inst, test_succ , comp, &ncomp , xheu);
		constr_gen( NULL , NULL , context, xstar, inst, succ, comp, &ncomp);
		
		int *ind = (int *) malloc(inst->ncols * sizeof(int));
		for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;
		
		//set the solution into cplex
		if(CPXcallbackpostheursoln(context, inst->ncols, ind , xheu , inst->solution.obj_best, CPXCALLBACKSOLUTION_NOCHECK)){
			print_error("An error occured on CPXcallbackpostheursoln");
		}

		//if ( CPXcallbackrejectcandidate(context, 0, NULL, NULL, NULL, NULL, NULL, NULL) ) print_error("CPXcallbackrejectcandidate() error"); // just reject the solution without adding cuts (less effective)
	}else if(ncomp == 1){
		updateSol(inst , succ , xheu , 1);

        gettimeofday(&end, 0);
        double elapsed = get_elapsed_time(start, end);
		inst->solution.solve_time = elapsed;

		int *ind = (int *) malloc(inst->ncols * sizeof(int));
		for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;

		//set the solution into cplex
		if(CPXcallbackpostheursoln(context, inst->ncols, ind , xheu , inst->solution.obj_best, CPXCALLBACKSOLUTION_NOCHECK)){
			print_error("An error occured on CPXcallbackpostheursoln");
		}
		
		
		
	}

	free(xheu);
    free(succ);
    free(comp);
	free(xstar);
	
	return 0; 
}    

int CPXPUBLIC cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) {
    instance *inst = (instance *) userhandle;
	// Initialize the mutex
    //pthread_mutex_init(&mutex, NULL);
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        return lazy_callback(context, contextid, inst);
    }
    if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        //return relaxation_callback(context, contextid, inst);
    }
	// Destroy the mutex
    //pthread_mutex_destroy(&mutex);
    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////
//user_cut////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


#include "tsp.h"
#include "heuristic.h"
#include "utility.h"
#include "debug.h"
#include "benders.h"
#include "vns.h"
#include "callbacks.h"
#include <math.h>
#include <sys/time.h>
#include <cplex.h>


int hardfixing(instance *inst){

    // open CPLEX model
	int error;
	//cplex environment pointer 
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	//cplex parameters pointer 
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model"); 
	if ( error ) print_error("CPXcreateprob() error");
	build_model(inst, env, lp);
	
	//initialyze the time of execution
	struct timeval start, end;
    gettimeofday(&start, 0);

	// Cplex's parameter setting
	//setParam(inst , env);
			
	inst->ncols =  CPXgetnumcols(env, lp);  // ncols == n variables into the model , ncol don't return error only the number of the vriables
    
	if( strncmp(inst->algorithm , "HARDFIX_CALL" , 12) == 0 ){
    	//call the callbacks when we have a incumbent solution and when there is a fractional solution (like gomory cut)
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION; // both lazy and usercuts
		//whit he upper call we obtain a flag for context id that 0 correspond to the first operator and 1 to the second 
		if ( CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst) ) print_error("CPXcallbacksetfunc() error"); 
    }

	//initialize a cplex initila solution for the 20% of time 
	int total_time = inst->timelimit;
	inst->timelimit = total_time*0.2;

	//grasp_algorithm(inst , 0);
	//twoopt_algorithm(inst);
	vns(inst , 0);

    //set remaining time 
	inst->timelimit = total_time - inst->timelimit;

	////initialize solution of cplex/////////////////////////////////////////
	initSol( inst , env, lp );
    setParam(inst , env);

    //set the best obj function found so far
    double best_obj = inst->solution.obj_best;

	double solve_time = inst->solution.solve_time;

    //init value of probability
    double prob = 0.6;

    inst->ncols = CPXgetnumcols(env, lp);
    
	//continue to resolve model until we reach the time limit 

    //printf("best obj : %f\n" , best_obj); 
	while(1){
        //double curr_sol = inst->solution.obj_best;
        //control remaining time
    	gettimeofday(&end, 0);
    	double time_diff = get_elapsed_time(start, end);
		//control if remain other time for execute another execution 
		if(time_diff > inst->timelimit ) print_error("Not sufficient time for find a solution!\n");
		//set the remaining time for find a new solution of the model 
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - time_diff); 
        
        

        //random_fix(env, lp, prob, &ncols_fixed, indexes, xstar);
        node_fix(env , lp  , inst , prob);
        

       //apply benders
	    if( strncmp( inst->algorithm , "HARDFIX_CALL" , 12) == 0 ) {
		    error = CPXmipopt(env,lp);
		    if ( error ) {

		        printf("CPX error code %d\n", error);
		        print_error("CPXmipopt() error");

		    }
	    }else{	
		    benders_solver(env , lp , inst);
	    }

        //printf("dopo , obj_best : %f , %f\n" , inst->solution.obj_best , best_obj);

        if( best_obj <= inst->solution.obj_best){

             break;

        }else{
            
			gettimeofday(&end, 0);
    		double elapsed = get_elapsed_time(start, end);
			best_obj = inst->solution.obj_best;
			solve_time = elapsed;
        } 

        reset_fix( env  , lp , inst , prob);
    
    }
    //close cplex
    closeCplex(env , lp);

	inst->solution.solve_time = solve_time;
	
	return 0;

}

int branching(instance *inst){

	// open CPLEX model
	int error;
	//cplex environment pointer 
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	//cplex parameters pointer 
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model"); 
	if ( error ) print_error("CPXcreateprob() error");
	build_model(inst, env, lp);
	
	//initialyze the time of execution
	struct timeval start, end;
    gettimeofday(&start, 0);

	// Cplex's parameter setting
	//setParam(inst , env);
			
	inst->ncols =  CPXgetnumcols(env, lp);  // ncols == n variables into the model , ncol don't return error only the number of the vriables
    
	if( strncmp(inst->algorithm , "BRANCHING_CALL" , 12) == 0 ){
    	//call the callbacks when we have a incumbent solution and when there is a fractional solution (like gomory cut)
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION; // both lazy and usercuts
		//whit he upper call we obtain a flag for context id that 0 correspond to the first operator and 1 to the second 
		if ( CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst) ) print_error("CPXcallbacksetfunc() error"); 
    }

	//initialize a cplex initila solution for the 20% of time 
	int total_time = inst->timelimit;
	inst->timelimit = total_time*0.2;

	//grasp_algorithm(inst , 0);
	//twoopt_algorithm(inst);
	vns(inst , 0);

    //set remaining time 
	inst->timelimit = total_time - inst->timelimit;

	////initialize solution of cplex/////////////////////////////////////////
	initSol( inst , env, lp );
    setParam(inst , env);

	double actual_cost;
	int K = 10;
	int iter = 0 ;

	//variable for constraint
	int *index = (int *)calloc(inst->nnodes , sizeof(int));
	double *value = (double *)calloc(inst->nnodes , sizeof(double));
	char **cname = ( char **)calloc(1 , sizeof(char *));
	cname[0] = (char *)calloc(100 , sizeof(char));
	sprintf(cname[0] , "branching constraint");
	double *xstar = (double *)calloc(inst->ncols , sizeof(double));
	double best_obj = inst->solution.obj_best;
	double solve_time;

	while(1){

		//save the actually best value 
		//actual_cost = inst->solution.obj_best;
		initSol( inst , env, lp );

		double time_diff = get_elapsed_time(start, end);
		//control if remain other time for execute another execution 
		//if(time_diff > inst->timelimit ) break;
		//set the remaining time for find a new solution of the model 
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - time_diff); 

		//add violation constraint 
		double rhs = inst->nnodes - K;
		char sense = 'G';
		int nnz = 0 ;
		bzero(xstar , inst->ncols*sizeof(double));

		//change the obtained solution in cplex form
		for ( int i = 0; i < inst->nnodes; i++ ){
			xstar[xpos(i,inst->solution.edges[i].j,inst)] = 1.0;
		} 

		//find the nonzero variables 
		for(int i = 0 ; i < inst->ncols ; i++){

			if(xstar[i]<1.0) continue;
			index[nnz] = i;
			value[nnz] = 1.0;
			nnz++;

		}

		int izero = 0;
		//add the constraint into cplex
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");

		//solve the problem
		if( strncmp( inst->algorithm , "BRANCHING_CALL" , 14) == 0 ) {
		    error = CPXmipopt(env,lp);
		    if ( error ) {

		        printf("CPX error code %d\n", error);
		        print_error("CPXmipopt() error");

		    }
	    }else{	
		    benders_solver(env , lp , inst);
	    }

		//remove the constraint added 
		int nrows = CPXgetnumrows(env , lp);
		if( CPXdelrows(env , lp , nrows -1 , nrows - 1 )) print_error(" Error in dleting rows");

		if(inst->solution.obj_best >= best_obj){
			K+=10;
			if(K > 40) break;
		}else{
			gettimeofday(&end, 0);
    		double elapsed = get_elapsed_time(start, end);
			best_obj = inst->solution.obj_best;
			solve_time = elapsed;
		} 


	}

	inst->solution.solve_time = solve_time;

	free(xstar);
	free(cname[0]);
	free(cname);
	free(value);
	free(index);
	
	return(0);

}



 void node_fix(CPXENVptr env, CPXLPptr lp, instance *inst, double prob){

    double rand_num;

    char lb = 'L'; //change only the lower bound
    double lower_val = 1.0; //value f lower bound
    
    //control if the probability is valid 
    //if(prob < 0 || prob > 1) { print_error("probability must be in [0,1]"); }
    //set the new lower bound for a random variables 
    for(int i = 0; i < inst->nnodes; i++) {
		rand_num = rand01();
		if(rand_num < prob) {
            //printf("set to 1 node : %d , %d\n" , i , inst->solution.edges[i].j);
            //fix the edge
            int pos = xpos(i , inst->solution.edges[i].j , inst);
            //set the new lower bound 
            if(CPXchgbds(env, lp, 1, &pos, &lb, &lower_val)) print_error("error in setting new lower bound");
		}
	}

 }


 void reset_fix(CPXENVptr env, CPXLPptr lp, instance *inst, double prob){

    double rand_num;

    char lb = 'L'; //change only the lower bound
    double lower_val = 0.0; //value f lower bound

    //set the new lower bound for a random variables 
    for(int i = 0; i < inst->ncols; i++) {
		
        int pos = i;
		//if(xstar[i] > 0.5) {
            //set the new lower bound 
            if(CPXchgbds(env, lp, 1, &pos, &lb, &lower_val)) print_error("error in setting new lower bound");
		//}
	}

 }


 
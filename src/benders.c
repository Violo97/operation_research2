#include <cplex.h> 
#include <pthread.h>
#include "utility.h"
#include "heuristic.h"
#include "tsp.h"   
#include "callbacks.h"  

#define COMPARE_VAL 0.999
#define FLAG_CALLBACK 1

pthread_mutex_t mutex;

//update solution 
void updateSol(instance *inst , int *succ , double *xheu , int update){

	

	instance tmp_inst;
	// init current best sol

	tmp_inst.nodes = calloc(inst->nnodes, sizeof(point));
    memcpy(tmp_inst.nodes, inst->nodes, sizeof(point) * inst->nnodes);

	tmp_inst.nnodes = inst->nnodes;
	tmp_inst.solution.obj_best = 0.0;


	//update path 
	double incumbent = 0.0;
	int next_node;
	tmp_inst.solution.edges = calloc(inst->nnodes , sizeof(edge));
	
	int node = succ[0];
	while(1){

		tmp_inst.solution.edges[node].i = node;
		tmp_inst.solution.edges[node].j = succ[node];
		tmp_inst.solution.obj_best += eucl_dist(&tmp_inst , node , succ[node]);
		node = succ[node];
		if(node == succ[0]) break;
	}
	
	//apply a 2-opt refinement
	twoopt_algorithm(&tmp_inst);
	//printf("=======> obj _ best : %f , incumbent : %f\n" , inst->solution.obj_best , tmp_inst.solution.obj_best );
	//control if the new best solution is better than the previous
	if( incumbent < inst->solution.obj_best && update){

		
		inst->solution.obj_best = tmp_inst.solution.obj_best;
		memcpy(inst->solution.edges, tmp_inst.solution.edges, inst->nnodes * sizeof(edge));
		//twoopt_algorithm(inst);

	}

	if(strncmp(inst->algorithm , "CALLBACK" , 8) == 0 || strncmp(inst->algorithm , "HARDFIX_CALL" , 12) == 0){
		for ( int i = 0; i < inst->nnodes; i++ ){
			xheu[xpos(i,tmp_inst.solution.edges[i].j,&tmp_inst)] = 1.0;
		} 

	} 
	
	free(tmp_inst.solution.edges);
	free(tmp_inst.nodes);
	
	
}


void benders_solver(CPXENVptr env , CPXLPptr lp ,  instance *inst){
	int error;
    double objval;

	//initialyze the time of execution
	struct timeval start, end;
    gettimeofday(&start, 0);

	//continue to resolve model until we reach a tsp solution or the time is end 
	while(1){

		error = CPXmipopt(env,lp);
		if ( error ) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		//obtain the model results into xstar
		double *xstar = (double *) calloc(inst->ncols, sizeof(double));
		if ( CPXgetx(env, lp, xstar, 0, inst->ncols-1) ) print_error("CPXgetx() error");

		
		/*for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				if ( xstar[xpos(i,j,inst)] > 0.5 ) printf("  ... x(%3d,%3d) = 1\n", i+1,j+1);
			}
		}*/

		int *succ = calloc(inst->nnodes , sizeof(int));
		int *comp = calloc(inst->nnodes , sizeof(int));
		//number of component into cplex solution 
		int ncomp = 0;

		for ( int i = 0; i < inst->nnodes; i++ )
		{
			succ[i] = -1;
			comp[i] = -1;
		}
		//insert verbose
        //printf("ncomp : %d\n" , ncomp);
	
		//build a solution for apply new constraints 
		build_sol(xstar, inst,succ, comp, &ncomp);

		//control if obtain a tsp solution
		if( ncomp == 1) {
			gettimeofday(&end, 0);
			double elapsed = get_elapsed_time(start, end);
			updateSol(inst, succ , NULL , 1);
			inst->solution.solve_time = elapsed;
			break;
		}	
		//printf("numbers of the components : %d \n" , ncomp);
		//for(int i = 0 ; i < inst->nnodes ; i++) printf("succ value : %d , comp value : %d\n" , succ[i] , comp[i]);

		//apply patching heuristic 
		heur_patch(inst, succ, comp, &ncomp);

		//define new constraints
		constr_gen( env , lp , NULL , xstar, inst, succ, comp, &ncomp);


		//control remaining time
    	gettimeofday(&end, 0);
    	double time_diff = get_elapsed_time(start, end);
		CPXwriteprob(env, lp, "log/partial_model.lp", NULL); 
		//getchar(); 
		//control if remain other time for execute another execution if we dont apply patching heur return error otherwhise 
		//return the best solution found with patching heur 
		if(time_diff > inst->timelimit ) print_error("Not sufficient time for find a solution!\n");
		//set the remaining time for find a new solution of the model 
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - time_diff); 


		free(succ);
		free(comp);
		free(xstar);
	}
    CPXgetobjval(env, lp, &objval);

}


//apply benders' loop method
int benders(instance *inst){

	// open CPLEX model
	int error;
	//cplex environment pointer 
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	//cplex parameters pointer 
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model"); 
	if ( error ) print_error("CPXcreateprob() error");

	// Initialize the mutex
    pthread_mutex_init(&mutex, NULL);

	build_model(inst, env, lp);

	//set cplex parameters
	setParam(inst  , env );

	//control if active callbacks
	if( strncmp( inst->algorithm , "CALLBACK" , 8) == 0 ) {

		//call the callbacks when we have a incumbent solution and when there is a fractional solution (like gomory cut)
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // both lazy and usercuts
		//whit he upper call we obtain a flag for context id that 0 correspond to the first operator and 1 to the second 
		if ( CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst) ) print_error("CPXcallbacksetfunc() error"); 

	}
	if( strncmp( inst->algorithm , "USER_CALLBACK" , 13) == 0 ) {

		//call the callbacks when we have a incumbent solution and when there is a fractional solution (like gomory cut)
		CPXLONG contextid = CPX_CALLBACKCONTEXT_RELAXATION; // both lazy and usercuts
		//whit he upper call we obtain a flag for context id that 0 correspond to the first operator and 1 to the second 
		if ( CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst) ) print_error("CPXcallbacksetfunc() error"); 

	}

	//set the number of columns
	inst->ncols =  CPXgetnumcols(env, lp);  // ncols == n variables into the model , ncol don't return error only the number of the vriables

	//initialize a cplex initila solution
	grasp_algorithm(inst , 0);
	twoopt_algorithm(inst);
	//set solution into cplex
	initSol(inst , env , lp);
	
	//we can move it in other function
	//apply benders
	if( strncmp( inst->algorithm , "CALLBACK" , 8) == 0 ) {
		error = CPXmipopt(env,lp);
		if ( error ) {

		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error");

		}
	}else{	
		benders_solver(env , lp , inst);
	}

	// Destroy the mutex
    pthread_mutex_destroy(&mutex);

	//close cplex model
	closeCplex(env , lp);

}	


void setParam(instance *inst , CPXENVptr env){

	//CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	//if ( VERBOSE >= 60 ) 
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randomseed);	
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit); 

}

void closeCplex(CPXENVptr env, CPXLPptr lp){

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

}

void initSol(instance *inst , CPXENVptr env, CPXLPptr lp ){

	////initialize solution of cplex/////////////////////////////////////////
	//----------------------------------------------------------------------//
	//convert a tsp solution into a cplex solution
	double *xheu = (double *) calloc(inst->ncols, sizeof(double));  // all zeros, initially
	for ( int i = 0; i < inst->nnodes; i++ ) xheu[xpos(i,inst->solution.edges[i].j,inst)] = 1.0;

	//adding a mip start solution
	int *ind = (int *) malloc(inst->ncols * sizeof(int));
	for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;
	int effortlevel = CPX_MIPSTART_NOCHECK;  
	int beg = 0; 						
	if ( CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");

	free(ind);
	free(xheu);
	//----------------------------------------------------------------------//


}

void heur_patch(instance *inst, int *succ, int *comp, int *ncomp , double *xheu){
	//generate a temporal comp array 
	int *tmp_comp = calloc(inst->nnodes , sizeof(int));
	//memcpy(tmp_comp, comp, inst->nnodes*sizeof(int));
	for(int i = 0 ; i< inst->nnodes ; i++) tmp_comp[i] = comp[i]; 
	//printf("value tmpo : %d\n" , tmp_comp[i]);


	//generate a tsp solution fuse all components to comp[1]
	for ( int n = 2 ; n <= *ncomp ; n++ ){
		double delta = INFINITY;
		int node_a = -1;
		int node_b = -1;
		for(int i = 0 ; i < inst->nnodes ; i++){
			//control il the point i is into the first component 
			if(tmp_comp[i] == 1){
				//control al the component of the nth componente for find the best delta 
				for(int j = i+1 ; j < inst->nnodes ; j++){

					if(tmp_comp[j] == n){

						//compute delta value 
						double delta_candidate = eucl_dist(inst , succ[i] , succ[succ[j]]) + eucl_dist(inst , succ[j] , succ[succ[i]])
												- eucl_dist(inst , succ[i] , succ[succ[i]]) - eucl_dist(inst ,  succ[j] , succ[succ[j]]);
						//printf("i value : %d j value : %d\n" , i , j);						
						//printf("delta candidate : %f\n" , delta_candidate);						
						//update delta 
						if( delta_candidate < delta ){
							delta = delta_candidate;
							node_a = i;
							node_b = j;
						}

					}

				}

			}

		}
		//update path
		if(node_a == -1 || node_b == -1) print_error("swapping node not found");
		int tmp = succ[succ[node_a]];
		succ[succ[node_a]] = succ[succ[node_b]];
		succ[succ[node_b]] = tmp;

		//update all the element into the merged component
		for(int j = 0 ; j < inst->nnodes ; j++){

			if(tmp_comp[j] == n) tmp_comp[j] = 1;

		}

			
	}
	
	updateSol(inst , succ , xheu , 0);

	free(tmp_comp);


}


void constr_gen( CPXENVptr env, CPXLPptr lp ,CPXCALLBACKCONTEXTptr context, double *xstar, instance *inst, int *succ, int *comp, int *ncomp){

	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

	//position of the coefficients
	int *index = (int *) calloc(inst->ncols, sizeof(int));
	//value of the coefficience
	double *value = (double *) calloc(inst->ncols, sizeof(double));

	char sense = 'L'; // L less or equals 

	//control if the the number of components is 2
	if(*ncomp == 2 ) *ncomp = 1;

	//for each components generate a strong constraint
	for( int n = 1 ; n <= *ncomp ; n++){
		sprintf(cname[0], "new constraint(%d)" , n);
		//define value of right hand side 
		double rhs = 0.0;
		//define value of non-zero elements 
		int nnz = 0;
		
		//for each component find a strong constraint to apply 
		for( int i = 0 ; i < inst->nnodes ; i++){

			if( comp[i] != n ) continue;
			rhs+=1.0;

			for( int j = 0 ; j < inst->nnodes ; j++){

				if( j>i && comp[j] == n){

					index[nnz] = xpos(i , j , inst);
					value[nnz] = 1.0;
					nnz++;
			
				}

			}

		}
		rhs-=1.0;
		int izero = 0;
		//printf("%f\n" , rhs);

		//control if is a valid cut 
		if(cut_violation( xstar , index , value , sense , rhs , nnz ) < 0){ print_error("Cut not valid \n");}

		//add new constraint
		if(context == NULL){
			if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
		}else{
			if ( CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value) ) print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
		}
	}

}
	
//pass all the relevant information 
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = eucl_dist(inst,i,j); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

// add the degree constraints 
	//position of the coefficients
	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	//value of the coefficience
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint  
		sprintf(cname[0], "degree(%d)", h+1); 
		//number of non-zeros  
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	//if ( VERBOSE >= 100 ) 
	CPXwriteprob(env, lp, "log/model.lp", NULL);   

}

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
{   

	//debug phase
	control_cplex_sol(xstar,inst);
	*ncomp = 0;
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
		
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}



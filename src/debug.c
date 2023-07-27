#include "tsp.h"
#include "utility.h"

void visit_all_node(instance *inst , int start_node);
void control_obj(instance *inst , int start_node);
void control_obj_2opt(instance *inst , int start_node);

//change while with a for because cant found a loop that not return into start point 
void visit_all_node(instance *inst , int start_node){

    //initialize the vector for control the node that the solution had visit
    int *visited_node = calloc(inst->nnodes , sizeof(int));

    //for(int i = 0 ; i < inst->nnodes ; i++) printf("init node -> %d , final node -> %d \n" , inst->solution.edges[i].i , inst->solution.edges[i].j);


    //set visited the starting node
    visited_node[start_node] ++;

    //set the next visited node 
    int node = inst->solution.edges[start_node].j;
    //control if there is multiple loop
    int count = 0;
    
    //analyze the node until we reach the initial node
    while(node!=start_node){
        if(count == inst->nnodes) print_error(" Path error ");
        count ++;
        //printf("dentro");
        visited_node[node] ++;
        //printf("init node -> %d , final node -> %d\n" , node , inst->solution.edges[node].j);
        //set the next solution node 
        node = inst->solution.edges[node].j;
        


    }

    //control if all the node are visited
    for(int i = 0 ; i < inst->nnodes ; i++){

        //if the node is visited 0 or more tan 1 times return an error 
        if(visited_node[i] != 1){
            //printf("missing node : %d\n" , i );
            print_error("Infisible solution");
        }

    }



}

void control_obj(instance *inst , int start_node){
    int count = 0;
    //value of obj_function
    double control_obj = 0; 
    int s_node = start_node;
    while(1){
        count++;
        control_obj += eucl_dist(inst , s_node , inst->solution.edges[s_node].j);
        //printf("valori dei nodi : %d -> %d , dist : %f , valore di obj = %f\n" , s_node , inst->solution.edges[s_node].j , eucl_dist(inst , s_node , inst->solution.edges[s_node].j) , control_obj);
        //printf("%d , " , inst->solution.edges[s_node].j);
        s_node = inst->solution.edges[s_node].j;
        if(s_node == start_node) break;

    }
    //printf("valore di onj = %f  , %f\n" , control_obj , eucl_dist_root(control_obj));
    //control_obj = eucl_dist_root(control_obj);
    //printf("\nf%f , %f , %f , %d , %d" , fabs(control_obj - inst->solution.obj_best) , control_obj , inst->solution.obj_best , count , start_node);
    if( fabs(control_obj - inst->solution.obj_best) > EPS  ) print_error(" objective function value error");
    //if( inst->solution.obj_best == nan || inst->solution.obj_best == inf)
}

void control_cplex_sol(const double *xstar, instance *inst){

    int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0) > EPS ) print_error(" wrong xstar in build_sol()");
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);

}



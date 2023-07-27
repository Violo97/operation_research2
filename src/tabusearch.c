#include "heuristic.h"
#include "utility.h"
#include "tsp.h"
#include "debug.h"
#include "float.h" 
#include <sys/time.h>

#define HIGH_TENURE 50
#define MIN_TENURE 10
#define NUM_ITER 100 // It' the number of iterations where the tenure changes

void tabu(instance *inst , int start_n);

void tabu( instance *inst , int start_n){

    //initialize the tabu list
    int *tabu_list = calloc(inst->nnodes , sizeof(int));
    //init current best sol
    edge *best_sol = calloc(inst->nnodes , sizeof(edge)); 

    for( int i = 0 ; i < inst->nnodes ; i++){
        tabu_list[i] = -(HIGH_TENURE + 1);
    }

    //init time
    struct timeval start, end;
    gettimeofday(&start, 0);

    //generate a solution 
    grasp_algorithm(inst , start_n);

    if( inst->verbose <= 20 ) printf("initial obj func : %f\n" , inst->solution.obj_best);

    double best_obj = inst->solution.obj_best;
    memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));

    //number of iteration 
    int iter = 0;
    int tenure = MIN_TENURE;

    while(1){

        //best node
        int besti = -1;
        int bestj = -1;

        //control if apply step policy
        //calc tenure value - costant method
        if( strncmp(inst->algorithm , "STEP_TABU" , 9) == 0 ){
            tenure = (iter/NUM_ITER)%2 ? HIGH_TENURE/5 : HIGH_TENURE;
        }

        if( strncmp(inst->algorithm , "LIN_TABU" , 8) == 0 ){

            if(tenure < HIGH_TENURE) tenure++;
            else tenure = MIN_TENURE;

        }

        //calc tenure value - random method
        if( strncmp(inst->algorithm , "RAND_TABU" , 9) == 0 ){

             if ( iter == 1 || iter % NUM_ITER == 0 ) {
                tenure = rand_choice(MIN_TENURE, HIGH_TENURE + 1); 
            } 

        } 

        //calc tenure value - dynamic change 

        double best_delta = DBL_MAX;

        //find the best node to swap 
        for( int i = 0 ; i < inst->nnodes ; i++){
            for( int j = i+1 ; j < inst->nnodes ; j++){
                //control if the two node is consecutive
                if( inst->solution.edges[j].j == i || inst->solution.edges[i].j == j) continue;

                //calc the delta value
                double delta = eucl_dist(inst , i, j) + eucl_dist(inst , inst->solution.edges[i].j, inst->solution.edges[j].j) 
                            - eucl_dist(inst ,i, inst->solution.edges[i].j) - eucl_dist(inst , j, inst->solution.edges[j].j); 

                //control if the nodes are a tabu
                if( iter - tabu_list[i] <= tenure || iter - tabu_list[j] <= tenure){

                    continue;

                }
                //control if the obj func is better than the other
                if(delta < 0 && best_delta > delta){
                    best_delta = delta;
                    besti = i;
                    bestj = j;
                }
            }
        }

        if(besti == -1 || bestj == -1) break; 

        //swap nodes
        int *partial_sol = calloc(inst->nnodes , sizeof(int)); 
        for (int i = 0; i < inst->nnodes; i++) {
            partial_sol[inst->solution.edges[i].j] = i;
        }
        //memcpy(partial_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
        int besti_f = inst->solution.edges[besti].j;
        int bestj_f = inst->solution.edges[bestj].j;
        inst->solution.edges[besti].j = bestj;
        inst->solution.edges[besti_f].j = bestj_f;

        inverse_path(inst , bestj , besti_f , partial_sol);

        free(partial_sol);
        //update obj function
        inst->solution.obj_best += best_delta;
        //print_sol(inst , start_n);
        visit_all_node(inst , start_n);
        control_obj(inst , start_n);

        //control elapsed time
        gettimeofday(&end, 0);
        double elapsed = get_elapsed_time(start, end);
        
        //control which is the best solution
        if( inst->solution.obj_best < best_obj ){
            best_obj = inst->solution.obj_best;
            if( inst->verbose <= 20 ) printf("update obj func : %f at iteration : %d at time: %f\n" , inst->solution.obj_best , iter , elapsed);
            memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
            inst->solution.solve_time = elapsed;
        }
        

        //set tabu node 
        tabu_list[besti] = iter;
        tabu_list[bestj] = iter;

        iter++;

        if (elapsed > inst->timelimit) {
            break;
        }

    }
    //debug 
    visit_all_node(inst , start_n);
    //save best solution
    inst->solution.obj_best = best_obj;
    memcpy(inst->solution.edges, best_sol, inst->nnodes * sizeof(edge));
    free(best_sol);

    
}
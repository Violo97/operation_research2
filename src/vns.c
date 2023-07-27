#include "utility.h"
#include "heuristic.h"
#include "debug.h"
#include "plot.h"
#include <sys/time.h>

void vns(instance *inst, int start_n);
void kick(instance *inst , int start_n);

void vns(instance *inst, int start_n)
{
    //save the better solution found
    edge *best_sol = calloc(inst->nnodes , sizeof(edge));
    double best_obj = inst->solution.obj_best;

    // Start counting time
    struct timeval start, end;
    gettimeofday(&start, 0);

    // test with greedy nearest neighbour
    grasp_algorithm(inst , start_n);
    twoopt_algorithm(inst);

    if( inst->verbose <= 10){
        printf("Initialize a solution with obj : %f \n" , inst->solution.obj_best);
    }

    // copy the current best sol into curr_sol
    memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));

    while (1)
    {

        // make a solution perturbation
        kick(inst, start_n);

        // apply 2-opt algorithm
        twoopt_algorithm(inst);
        
        if( inst->verbose <= 10){
            printf("--Obj After 2-opt Refinement : %f \n" , inst->solution.obj_best);
        }

        // control elapsed time
        gettimeofday(&end, 0);
        double elapsed = get_elapsed_time(start, end);

        // control if the new solution is a minimum better than the old
        if (inst->solution.obj_best < best_obj)
        {

            // update objective func 
            best_obj = inst->solution.obj_best;

            inst->solution.solve_time = elapsed;
            // update the path
            memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));

            if( inst->verbose <= 20){
                printf("---Update Better Obj : %f , at time : %f s\n" , inst->solution.obj_best , elapsed);
            }
            
        }else{

            // restore best solution
            inst->solution.obj_best = best_obj;
            memcpy(inst->solution.edges, best_sol, inst->nnodes * sizeof(edge));

        }

        if (elapsed > inst->timelimit)
        {
            break;
        }
        
    }

    free(best_sol);
}

void kick(instance *inst , int start_n)
{

    int *candidate_n = calloc(5, sizeof(int));
    int *cand_succ = calloc(5 , sizeof(int));
    int n_cand = 0;
    
    //generate the 5 node to swap
    while (n_cand < 5){
    
        candidate_n[n_cand] = rand() % inst->nnodes;
        for (int i = 0; i < n_cand; i++){
            // check if the random value is valid
            if (candidate_n[n_cand] == candidate_n[i] || candidate_n[n_cand] == inst->solution.edges[candidate_n[i]].j || inst->solution.edges[candidate_n[n_cand]].j == candidate_n[i]){
                n_cand--;
                break;
            }
        }
        n_cand++;
    }

    //reorder the nodes 
    n_cand = 0 ;
    int node = start_n;
    while(n_cand < 5){

        for(int i = 0 ; i < 5 ; i++){

            if(candidate_n[i] == node){
                int tmp = candidate_n[n_cand];
                candidate_n[n_cand] = node;
                candidate_n[i] = tmp;
                cand_succ[n_cand] = inst->solution.edges[node].j;
                n_cand++;
                break;
            }

        }
        node = inst->solution.edges[node].j;
        

    }

    inst->solution.edges[candidate_n[3]].j = cand_succ[4];
    inst->solution.edges[candidate_n[1]].j = cand_succ[2];
    inst->solution.edges[candidate_n[4]].j = cand_succ[0];
    inst->solution.edges[candidate_n[2]].j = cand_succ[3];
    inst->solution.edges[candidate_n[0]].j = cand_succ[1];
    
    if( inst->verbose == 0){
        visit_all_node(inst , 0);
    }    

    // recompute the cost of the path
    double obj_cost = 0;
    // set the initial node
    node = start_n;
    while (1)
    {
        obj_cost += eucl_dist(inst, node, inst->solution.edges[node].j);
        node = inst->solution.edges[node].j;
        if (node == start_n) break;
        
    }
    // set the new value
    inst->solution.obj_best = obj_cost;
    if( inst->verbose <= 10){
        printf("--Obj After 5-Kick : %f \n" , inst->solution.obj_best);
    }
    free(candidate_n);
    free(cand_succ);
}
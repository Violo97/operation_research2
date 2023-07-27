#include "heuristic.h"
#include "utility.h"
#include "tsp.h"
#include "debug.h"
#include <sys/time.h>

#define ALFA 0.95
#define T_MAX 100
#define T_MIN 5
#define MAX_ITER 1000
#define N_ANNEALING 10
#define EULER 2.718

void annealing(instance *inst, int start_n);
void annealing2opt(instance *inst, double temperature, int iter);
double manhattan(instance *inst, double temperature, int iter, double delta);

void annealing(instance *inst, int start_n)
{

    // init current best sol
    edge *best_sol = calloc(inst->nnodes, sizeof(edge));

    // init time
    struct timeval start, end;
    gettimeofday(&start, 0);

    // generate a solution
    grasp_algorithm(inst, start_n);
    twoopt_algorithm(inst);

    if( inst->verbose <= 20 ) printf("initial obj val : %f\n" , inst->solution.obj_best);

    //save initial solution
    double best_obj = inst->solution.obj_best;
    memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));

    // execute n annealing
    for (int i = 0; i < N_ANNEALING; i++)
    {
        printf("INIT ANNEALING\n");
        // define the initial value of temperature
        double temperature = T_MAX;

        // number of iteration
        int iter = 0;
        while (1)
        {   

            // choose the two random node to switch
            annealing2opt(inst, temperature, iter);
            // update temperature value
            if (temperature > T_MIN)
                temperature *= ALFA;

            // control elapsed time
            gettimeofday(&end, 0);
            double elapsed = get_elapsed_time(start, end);

            // control the new solution if it is the best sol
            if (inst->solution.obj_best < best_obj)
            {
                if(inst->verbose <= 20) printf("new obj : %f\n", inst->solution.obj_best);
                best_obj = inst->solution.obj_best;
                memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
                inst->solution.solve_time = elapsed;
            }

            iter++;

            if (elapsed > inst->timelimit || iter == MAX_ITER)
            {
                break;
            }
        }
    }

    //restore the best solution
    inst->solution.obj_best = best_obj;
    memcpy(inst->solution.edges , best_sol , inst->nnodes * sizeof(edge));
}

void annealing2opt(instance *inst, double temperature, int iter)
{
    
    int first_node;
    int second_node;
    double delta;

    while (1)
    {

        // choose the first node
        first_node = rand() % inst->nnodes;

        // choose the second node
        second_node = rand() % inst->nnodes;

        // control if the second node is valid
        while (second_node == first_node || second_node == inst->solution.edges[first_node].j || inst->solution.edges[second_node].j == first_node)
        {
            // choose another node
            second_node = rand() % inst->nnodes;
        }

        // calc the delta value
        delta = eucl_dist(inst, first_node, second_node) + eucl_dist(inst, inst->solution.edges[first_node].j, inst->solution.edges[second_node].j) - eucl_dist(inst, first_node, inst->solution.edges[first_node].j) - eucl_dist(inst, second_node, inst->solution.edges[second_node].j);

        // control if the delta value is negative
        if (delta < 0)
            break;
        // if delta is positive use metropolis formula
        double manhattan_value = manhattan(inst, temperature, iter, delta);
        //printf("manhattan value : %f \n", manhattan_value);
        if (rand01() < manhattan_value)
            break;
    }

    // update obj value
    inst->solution.obj_best += delta;
    if(inst->verbose <= 10){printf("new incumbent : %f\n" , inst->solution.obj_best);}

    // update path
    // swap nodes
    int *partial_sol = calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++)
    {
        partial_sol[inst->solution.edges[i].j] = i;
    }
    int besti_f = inst->solution.edges[first_node].j;
    int bestj_f = inst->solution.edges[second_node].j;
    inst->solution.edges[first_node].j = second_node;
    inst->solution.edges[besti_f].j = bestj_f;
    inverse_path(inst, second_node, besti_f, partial_sol);

    if( inst->verbose == 0){
        control_obj(inst , 0);
        visit_all_node(inst , 0);
    }


    free(partial_sol);
}

double manhattan(instance *inst, double temperature, int iter, double delta)
{

    // compute power value
    double power = -1 * (delta / temperature);

    // compute the probability value for choose a bad delta
    double prob = pow(EULER, power);

    return prob;
}
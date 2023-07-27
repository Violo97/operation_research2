#include "heuristic.h"
#include "utility.h"
#include "tsp.h"
#include "debug.h"
#include <float.h>
#include <sys/time.h>
#include <stdlib.h>

#define POPULATION_SIZE 1000
#define PARENT_RATE 0.5
#define MUTATION_RATE 0.1
#define GRASP_PROB 0.5
#define TOURNAMENT_SIZE 20
#define TEST 1

void genetic(instance *inst);

//struct for save all the information of the population 
typedef struct {   

	// parameters 
	double fitness;
    int *genes;
	
	
} population;

int compare (const void* first, const void* second)
{
  const population* left_p = first;
  const population* right_p = second;

  return right_p->fitness - left_p->fitness;

}

//find the people with the best fitness
void calc_best_fitness(const population* people, double* best, int *best_idx) {
    *best = DBL_MAX;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        double fitness = people[i].fitness;
        if (fitness < *best) {
            *best = fitness;
            *best_idx = i;
        }
    }

}


void genetic_algo(instance *inst ){

    //Start counting time from now
    struct timeval start, end;
    gettimeofday(&start, 0);

    //initialize the population
    population *people = calloc(POPULATION_SIZE, sizeof(population));
    //initialize all the solution population 
    generate_population(inst , people);

    //Allocate memory for parents and child
    const int parent_size = (int) (POPULATION_SIZE * PARENT_RATE);
    int* parents = calloc(parent_size, sizeof(int)); 
    population* children = calloc(parent_size, sizeof(population));
    for (int i = 0; i < parent_size; i++) {
        children[i].genes = calloc(inst->nnodes, sizeof(int));
    }

    double best_fitness = DBL_MAX;
    int best_idx = 0;
    double incumbent = best_fitness;

    int generation = 1;
    //save best individual
    population best_individual;
    //Repeat until time limit is reached
    while (1) {
        //Check elapsed time
        gettimeofday(&end, 0);
        double elapsed = get_elapsed_time(start, end);
        if (elapsed > inst->timelimit) {
            break;
        }


        //Compute the best fitness value and the best individual's index in the population
        calc_best_fitness(people , &best_fitness, &best_idx);
        //If there is a tour in the current population that is better than the one seen so far, save it
        if (best_fitness < incumbent) {
            incumbent = best_fitness;
            best_individual = people[best_idx];
            inst->solution.obj_best = best_fitness;
            inst->solution.solve_time = elapsed;
            //update_sol(inst, best_individual); //Update best solution
            if(inst->verbose == 2) printf("The new best solution is : %f at generation : %d\n" , best_fitness , generation);
        }

        //if(inst->verbose <= 1) 
        printf("At generation : %d , the best solution is : %f  and the incumbent solution is : %f\n" ,generation , best_fitness , incumbent );
        
        //select the parents for the new generation
        printf("init select parent\n");
        select_parents(people, parents, parent_size);
        printf("finish select parent \n");
        //generate children by combining two parents
        procreate(inst, people, parents, parent_size, children);
        printf("finish to generate children \n");

        //apply the mutation into the population
        mutation(inst, children, parent_size);
        printf("finish to make mutation \n");

        //Replace the individuals of the current populations with the children that has better fitness
        choose_survivors(inst, people, children, parent_size , best_idx);
        printf("finish to generate new generation \n");

        generation++;

    }    

    //save best sol 
    update_sol(inst, best_individual); //Update best solution
    inst->solution.obj_best = incumbent;
}



void generate_population(instance *inst , population *people){
    
    //set the probability for grasp algorithm 
    //inst->prob = GRASP_PROB;

    for(int i = 0 ; i < POPULATION_SIZE ; i++){
        
        //choose at random the initial node
        int start_node = rand() % inst->nnodes;

        //chang the prob of grasp for each people
        //inst->prob = random01();

        //initialize the parents istance 
        people[i].genes = calloc(inst->nnodes, sizeof(int));

        //generate a solution 
        grasp_algorithm(inst , start_node);

        //save the solution of the instance
        int node_idx = start_node;
        int node_iter = 0;
        while (node_iter < inst->nnodes) {
            people[i].genes[node_iter++] = inst->solution.edges[node_idx].i;
            node_idx = inst->solution.edges[node_idx].j;
        }

        //set the fitness of this individual
        people[i].fitness = inst->solution.obj_best;


    }


}

void update_sol(instance* inst, population people) {
    for (int i = 0; i < inst->nnodes - 1; i++) {

        int index = people.genes[i];
        inst->solution.edges[index].i = people.genes[i];
        inst->solution.edges[index].j = people.genes[i + 1];

    }

    int index = people.genes[inst->nnodes - 1];
    inst->solution.edges[index].i = people.genes[inst->nnodes - 1];
    inst->solution.edges[index].j = people.genes[0];
}

//apply the tournament method for choose the parents 
void select_parents(population* people, int* parents, const int parent_size) {


    memset(parents, -1, parent_size * sizeof(int));
    int count = 0;
    int t_size = TOURNAMENT_SIZE;
    //int* visited = (int*)malloc(POPULATION_SIZE);
    int* visited = calloc(POPULATION_SIZE , sizeof(int));
    for(int i = 0 ; i<POPULATION_SIZE ; i++) visited[i] = 0;
    //printf("here\n");
    for( int i = 0 ; i < parent_size ; i++){
        //printf("%d\n" , i);
        double best_obj = INFINITY;
        int pos_best = -1;
        count = 0;
        while(count < t_size){
            int candidate = rand_choice(0 , POPULATION_SIZE);
            if( visited[candidate]!=1 ){
                //printf("%f , %d , $%d\n" , people[candidate].fitness , count , parent_size);
                if(best_obj > people[candidate].fitness) {
                    best_obj = people[candidate].fitness;
                    pos_best = candidate;
                }
                count++;
            }


        }

        visited[pos_best] = 1;
        parents[i] = pos_best;

    }

    free(visited);
    

    
}


void crossover(instance* inst, const population *people, const int parent1, const int parent2) {
    population p1 = people[parent1];
    population p2 = people[parent2];

    inst->solution.obj_best = 0;

    int *visited = calloc(inst->nnodes, sizeof(int));
    
    //select the random node from which define the truncation of the genes    
    int truncation_index = rand() % inst->nnodes;
    int index = 1;
    //set the first node of the child
    visited[p1.genes[0]] = 1;
    int prev_node = p1.genes[0];
    //chromosome[0] = p1.genes[0];
    inst->solution.edges[p1.genes[0]].i = p1.genes[0];
    for (int i = 1; i < inst->nnodes; i++) {
        int node;
        if (i <= truncation_index) {
            node = p1.genes[i];
            visited[node] = 1;
            inst->solution.edges[node].i = node;
            inst->solution.edges[prev_node].j = node;
            
            //chromosome[index] = node;
        } else {
            node = p2.genes[i];
            if (visited[node]) { continue; }
            inst->solution.edges[node].i = node;
            inst->solution.edges[prev_node].j = node;
            //chromosome[index] = node;
        }
        //*fitness_value += eucl_dist(inst , chromosome[index] , chromosome[index-1]);
        inst->solution.obj_best += eucl_dist(inst , inst->solution.edges[prev_node].i , inst->solution.edges[prev_node].j);
        prev_node = node;
        index++;
    }
    if (index < inst->nnodes) {
        for (int i = 0; i <= truncation_index; i++) {
            int node = p2.genes[i];
            if (visited[node]) { continue; }

            inst->solution.edges[node].i = node;
            inst->solution.edges[prev_node].j = node;
            inst->solution.obj_best += eucl_dist(inst , inst->solution.edges[prev_node].i , inst->solution.edges[prev_node].j);
            prev_node = node;
            index++;

            //chromosome[index++] = node;
            //*fitness_value += eucl_dist(inst , chromosome[index] , chromosome[index-1]);
        }
    }
    inst->solution.edges[prev_node].j = inst->solution.edges[p1.genes[0]].i;
    inst->solution.obj_best += eucl_dist(inst , inst->solution.edges[prev_node].i , inst->solution.edges[prev_node].j);
    //if(inst.verbose == 0)
    visit_all_node(inst , 0);
    control_obj(inst , 0);
    //make a 2opt refinement 
    twoopt_algorithm(inst);

    //*fitness_value += eucl_dist(inst , chromosome[index] , chromosome[0]);
    

    free(visited);
}


void procreate(instance* inst, const population* people, const int* parents, const int parent_size, population* children) {
    //save the best solution found so far
    //edge *best_sol = calloc(inst->nnodes , sizeof(edge)); 
    //double best_obj = inst->solution.obj_best;
    //memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));

    int counter = 0;
    double fitness_value;

    for (int i = 0; i < parent_size; i++) {
        fitness_value = 0;
        int j = (i + 1) % parent_size;
        int parent1 = parents[i];
        int parent2 = parents[j];
        //printf("parents 1 : %d , parents 2 : %d\n" , parent1 , parent2);
        //generate the child
        //printf("start crossover \n");
        crossover(inst, people, parent1, parent2);
        //printf("crossover complete \n");
        //copy the genes of child 
        //memcpy(children[counter].genes, chromosome, sizeof(int) * inst->nnodes);
        int node_idx = 0;
        int node_iter = 0;
        while (node_iter < inst->nnodes) {
            children[i].genes[node_iter++] = inst->solution.edges[node_idx].i;
            node_idx = inst->solution.edges[node_idx].j;
        }

        //update fitness
        children[i].fitness = inst->solution.obj_best;

        counter++;
    }

    //restore best sol
    //inst->solution.obj_best = best_obj;
    //memcpy(inst->solution.edges, best_sol, inst->nnodes * sizeof(edge));

    //free(best_sol);
}


void mutation(instance* inst, population* children , int parent_size) {

    //for each child control if the prob make a mutation
    for( int index = 0 ; index < parent_size ; index++){

        //control the prob
        if( rand01() < MUTATION_RATE){

            //choose two random node to swap 
            int rand1 = rand_choice(0 , inst->nnodes);
            int rand2 = rand_choice(0 , inst->nnodes);

            //control if the two node are the same 
            while(rand1 == rand2) rand2 = rand_choice(0 , inst->nnodes);

            int temp = children[index].genes[rand1];
            children[index].genes[rand1] = children[index].genes[rand2];
            children[index].genes[rand2] = temp;
            calc_fitness(inst, &(children[index]));

        }

        //if(inst.verbose == 0)
    visit_all_node(inst , 0);
    control_obj(inst , 0);



    }
}



void calc_fitness(instance* inst, population* people) {
    int prev_node = people->genes[0];
    people->fitness = 0;
    for (int i = 1; i < inst->nnodes; i++) {
        int node = people->genes[i];
        people->fitness += eucl_dist(inst , prev_node, node);
        prev_node = node;
    }
    people->fitness += eucl_dist(inst , prev_node, people->genes[0]);
}


void choose_survivors(instance* inst, population* people, const population* children, int parent_size , int best_idx) {

    //int* survivor = (int*)malloc(inst->nnodes);
    int* survivor = calloc(POPULATION_SIZE , sizeof(int));
    int count = 0;
    for(int i = 0 ; i<POPULATION_SIZE ; i++) survivor[i] = 0;

    for( int i = 0 ; i < POPULATION_SIZE - parent_size ; i++){

        double best_obj = INFINITY;
        int pos_best = -1;
        count = 0;
        while(count < TOURNAMENT_SIZE){

            int candidate = rand_choice(0 , POPULATION_SIZE);
            if( !survivor[candidate] ){
                if(best_obj > people[candidate].fitness) {
                    best_obj = people[candidate].fitness;
                    pos_best = candidate;
                }
                count++;
            }


        }

        survivor[pos_best] = 1;

    }
    printf("choose survivor\n");
    count = 0;
    for(int index = 0 ; index < POPULATION_SIZE ; index++){

        if(!survivor[index]){

            memcpy(people[index].genes, children[count].genes, sizeof(int) * inst->nnodes);
            people[index].fitness = children[count].fitness;
            count++;

        }

    }

    free(survivor);
    
}

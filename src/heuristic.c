#include "tsp.h"
#include "heuristic.h"
#include "utility.h"
#include "debug.h"
#include "float.h"
#include <math.h>
#include <sys/time.h>

void grasp_algorithm(instance *inst , int start_n);
void extramileage_algorithm(instance *inst);
void twoopt_algorithm(instance *inst);
void inverse_path(instance *inst , int final_p , int init_p , int *prev_sol);
void Grasp( instance *inst , int start_n);
void Extra_Mileage( instance *inst);




void grasp_algorithm(instance *inst , int start_n){
    //control if the starting node or path length is valid
    if(start_n > inst->nnodes - 1) print_error("NODE NOT VALID");

    //initialize the vector of visited node 
    int *visited = calloc(inst->nnodes , sizeof(int));
    //total distance
    double obj = 0;

    //set visited the starting node
    visited[start_n] = 1;
    //set the current node for analyze the nearest nodes
    int curr_node = start_n;
    
    double obj_update;
    double best_dist;
    double second_dist;
    //found a path with every node  
    while(1){
        best_dist = DBL_MAX;
        second_dist = DBL_MAX;
        //save the best nearest node found
        int best_near = -1;
        int second_near = -1;
        //for each non visited node control the nearest neighbor
        for(int i = 0 ; i < inst->nnodes ; i++){
            if(visited[i] != 1){
                double dist = eucl_dist(inst , curr_node , i);
                if(dist < second_dist){
                    if(dist < best_dist){
                        second_near = best_near;
                        second_dist = best_dist;
                        best_near = i;
                        best_dist = dist;
                    }else{
                        second_near = i;
                        second_dist = dist;
                    }
                    
                }
            }
        }

        // if we visited all nodes
        if (best_near == -1) {
            // Closing the cycle 
            inst->solution.edges[curr_node].i = curr_node;
            inst->solution.edges[curr_node].j = start_n;
            //update the objective func
            obj += eucl_dist(inst , curr_node , start_n);
            break;
        }
        
        if( rand01()<inst->prob && second_near!=-1){
            //set the edge between the current node and the next node
            inst->solution.edges[curr_node].i = curr_node;
            inst->solution.edges[curr_node].j = second_near;

            //set the new node
            curr_node = second_near;
            visited[curr_node] = 1;

            //update the objective function
            obj += second_dist;
        }else{
            //set the edge between the current node and the next node
            inst->solution.edges[curr_node].i = curr_node;
            inst->solution.edges[curr_node].j = best_near;

            //set the new node
            curr_node = best_near;
            visited[curr_node] = 1;

            //update the objective function
            obj += best_dist;
        }

    }
     
    //update obj funct
    inst->solution.obj_best = obj;

    //debug control 
    if( inst->verbose == 0){
        control_obj(inst , start_n);
        visit_all_node(inst , start_n);
    }    

    free(visited);

}

void extramileage_algorithm(instance *inst){

    //found the two starting node with the gratest distance between them 
    int start_n;
    int end_n;
    int initial_point;
    double obj = 0;
    for(int i = 0 ; i < inst->nnodes ; i++){
        for(int j = i+1; j < inst->nnodes ; j++){
            double dist = eucl_dist(inst , i , j);
            if(dist > obj){
                start_n = i;
                end_n = j;
                obj = dist;
            }
        }
    }

    //update the obj funct 
    obj = 2*obj;

    //set the two initial edge 
    inst->solution.edges[start_n].i = start_n;
    inst->solution.edges[start_n].j = end_n;
    inst->solution.edges[end_n].i = end_n;
    inst->solution.edges[end_n].j = start_n;

    initial_point = start_n;
    //create a vector for maintain the visited nodes
    int *visited = calloc(inst->nnodes , sizeof(int));

    //set vidited the two two nodes
    visited[start_n] = 1;
    visited[end_n] = 1;

    while(1){

        double min_path = DBL_MAX;
        for(int i = 0 ; i < inst->nnodes ; i++){

            if(visited[i]){ continue;}

            for(int j = 0 ; j < inst->nnodes ; j++){

                if(visited[j] == 1){

                    double delta = eucl_dist(inst , j , i) 
                                        + eucl_dist(inst , i , inst->solution.edges[j].j) 
                                        - eucl_dist(inst , inst->solution.edges[j].i , inst->solution.edges[j].j);
                    if(min_path > delta){
                        start_n = j;
                        end_n = i;
                        min_path = delta;

                    }                    


                }

            }

        }
        
        if(min_path == DBL_MAX){
            break;
        }
            
        //printf("%d \n" , end_n);
        visited[end_n] = 1;
        int final_point = inst->solution.edges[start_n].j;
        inst->solution.edges[start_n].j = end_n;
        inst->solution.edges[end_n].i = end_n;
        inst->solution.edges[end_n].j = final_point;
        obj += min_path;

        

    }

    //set new objective func
    inst->solution.obj_best = obj;
    //debug control 
    control_obj(inst , start_n);
    visit_all_node(inst , start_n);


    free(visited);

}



////////////////////////////////////////////////////////////////////////
////////////////  2-opt algorithm /////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void twoopt_algorithm(instance *inst){
    double best_obj = inst->solution.obj_best;
    int *prev_sol = calloc(inst->nnodes , sizeof(int));
    //set the edges into prev_sol
    for (int i = 0; i < inst->nnodes; i++) {
        prev_sol[inst->solution.edges[i].j] = i;
    }

    while(1){
        //for each pair of nodes
        for(int i = 0 ; i < inst->nnodes - 1 ; i++){
            for(int j = i+1 ; j < inst->nnodes ; j++){

                //save the two node that we want to switch
                int s_node_1 = i;
                int s_node_2 = j;
                int e_node_1 = inst->solution.edges[s_node_1].j;
                int e_node_2 = inst->solution.edges[s_node_2].j;

                // Skip non valid configurations 
                if (e_node_1 == e_node_2 || s_node_1 == e_node_2 || s_node_2 == e_node_1) {continue;}

                // Compute the delta 
                double delta = eucl_dist(inst , s_node_1, s_node_2) + eucl_dist(inst , e_node_1, e_node_2) - eucl_dist(inst ,s_node_1, e_node_1) - eucl_dist(inst , s_node_2, e_node_2); 
                //delta <0 crossing
                if(delta < 0){
                    //swap the two edges
                    inst->solution.edges[s_node_1].j = s_node_2;
                    inst->solution.edges[e_node_1].j = e_node_2;
                    //update path
                    inverse_path(inst , s_node_2 , e_node_1 , prev_sol);
                    inst->solution.obj_best += delta;
                    //debug 
                    if( inst->verbose == 0){
                        visit_all_node(inst , 0);
                        control_obj(inst , 0);
                    }    
                    //printf("two-opt update : %f \n" , inst->solution.obj_best);
                }                


            }
        }

        //if not found a new better solution stop the execution
        if(inst->solution.obj_best >= best_obj){break;}

        //update the cost 
        best_obj = inst->solution.obj_best;
        

    }

    free(prev_sol);
}

void inverse_path(instance *inst , int start_node , int end_node , int *prev_sol){
    int currnode = start_node;
    while (1) {
        int node = prev_sol[currnode];
        inst->solution.edges[currnode].j = node;
        currnode = node;
        if (node == end_node) {
            break;
        }
    }

    for (int i = 0; i < inst->nnodes; i++) {
        prev_sol[inst->solution.edges[i].j] = i;
    }

}

////////////////////////////////////////////////////////
//// Function for call heuristic ///////////////////////
////////////////////////////////////////////////////////



void Grasp( instance *inst , int start_n){

    
    int iter = 1;
    //save the best solution found
    edge *best_sol = calloc(inst->nnodes, sizeof(edge));
    memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
    double best_obj = inst->solution.obj_best;
    //start timer 
    struct timeval start, end;
    gettimeofday(&start, 0);
    while(1){
        //call freedy algorithm 
        grasp_algorithm(inst , start_n);

        if( inst->verbose <= 10){
            printf("---------------------------------\n");
            printf("-------- Iteration : %d --------\n" , iter);
            printf("---------------------------------\n\n");
            printf("--- Returned Value : %f \n\n" , inst->solution.obj_best);
            printf("---------------------------------\n");
            printf("---------------------------------\n\n");
        }

        if( strncmp(inst->algorithm , "GRASP_2_OPT" , 11) == 0  || strncmp(inst->algorithm , "GREEDY_2_OPT" , 12) == 0){
            //call two opt algorithm
            twoopt_algorithm(inst);
            if( inst->verbose <= 10){
                printf("--- Returned Value After Refinement : %f \n\n" , inst->solution.obj_best);
            }

        }

        //control the timer
        gettimeofday(&end, 0);
        double time_diff = get_elapsed_time(start, end);

        //control which is the better solution 
        if(best_obj > inst->solution.obj_best){

            best_obj = inst->solution.obj_best;
            inst->solution.solve_time = time_diff;
            memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
                
            if( inst->verbose <= 20){
                printf("--- Update Better Obj : %f \n\n" , inst->solution.obj_best);
            }

        }
            
        if(time_diff > inst->timelimit){

            memcpy(inst->solution.edges, best_sol, inst->nnodes * sizeof(edge));
            inst->solution.obj_best = best_obj;
            break;

        }

        iter++;
            
        }

}


void Extra_Mileage( instance *inst){

    
        int iter = 0;
        //save the best solution found
        edge *best_sol = calloc(inst->nnodes, sizeof(edge));
        memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
        double best_obj = inst->solution.obj_best;
        //start timer 
        struct timeval start, end;
        gettimeofday(&start, 0);
        while(1){
            
            

            //call algorithm 
            extramileage_algorithm(inst);

            if( inst->verbose <= 10){
                printf("---------------------------------\n");
                printf("-------- Iteration : %d --------\n" , iter);
                printf("---------------------------------\n\n");
                printf("--- Returned Value : %f \n\n" , inst->solution.obj_best);
                printf("---------------------------------\n");
                printf("---------------------------------\n\n");
            }

            if( strncmp(inst->algorithm , "EXTRAMILL_2_OPT" , 15) == 0  ){
                //call two opt algorithm
                twoopt_algorithm(inst);

                if( inst->verbose <= 10){
                    printf("--- Returned Value After Refinement : %f \n\n" , inst->solution.obj_best);
                }

            }

            //control the timer
            gettimeofday(&end, 0);
            double time_diff = get_elapsed_time(start, end);
    
            //control which is the better solution 
            if(best_obj > inst->solution.obj_best){

                best_obj = inst->solution.obj_best;
                inst->solution.solve_time = time_diff;
                memcpy(best_sol, inst->solution.edges, inst->nnodes * sizeof(edge));
                
                if( inst->verbose <= 20){
                    printf("--- Update Better Obj : %f \n\n" , inst->solution.obj_best);
                }

            }
            
            if(time_diff > inst->timelimit){

                memcpy(inst->solution.edges, best_sol, inst->nnodes * sizeof(edge));
                inst->solution.obj_best = best_obj;
                break;

            } 
            iter++;
            
        }

   
}


#include "tsp.h"
#include "heuristic.h"
#include "plot.h"
#include "benders.h"
#include "matheu.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>




//void gen_plot();
void read_input(instance *inst);
void generate_inst(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst , int *start_node , char* file_name);
//call the input algorithm 
void call_algo(instance *inst , int start_node ); 



int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); } 

	//file in ehich save the result
	char file_name[100];

	//set the default starting node to 0 
	int start_node = 0 ;
	instance inst;

	parse_command_line(argc,argv, &inst , &start_node , &file_name);

	//printf("name of file : %s \n" , file_name);
	
	//set a default seed value 
	if(inst.randomseed == 0){
		clock_t c = clock();
		inst.randomseed = c;	
	}    

	//generate or read the coordinates of the points  
	if(strncmp(inst.input_file, "NULL", 4) == 0){
        generate_inst(&inst);
    }else{
        read_input(&inst);
    } 

	//execute the algorithm 
	call_algo(&inst , start_node);

	FILE *file;
    
    // Open the file in append mode
    //file = fopen(file_name, "a");

    
    //if (file == NULL) {
    //    print_error("Unable to open the file.\n");
	//}
    
    // Write the number to the file
    //fprintf(file, "%lf\n", inst.solution.solve_time);
    
    // Close the file
    //fclose(file);

	if( inst.verbose <= 30){
        printf("Best Obj Found : %f , with : %f s\n" , inst.solution.obj_best , inst.solution.solve_time);
    }

	point_file(&inst);
	solution_file(&inst , start_node);
	
	gen_plot();
	//print_instance(&inst);
	free_instance(&inst);
	return 0; 
} 

void generate_inst(instance *inst){
 
    inst->nodes = calloc(inst->nnodes , sizeof(point));
	inst->solution.edges = calloc((inst->nnodes) , sizeof(edge));
	//inst->solution.best_sol = calloc((inst->nnodes+1) , sizeof(int));

	//generate a random graph 
	srand(inst->randomseed);
    for(int i = 0 ; i < inst->nnodes  ; i++){
        inst->nodes[i].x = (int)(((double) rand() / RAND_MAX) * 1000);
        inst->nodes[i].y = (int)(((double) rand() / RAND_MAX) * 1000);
    }

}

void call_algo(instance *inst , int start_node ){
	
	if( strncmp(inst->algorithm , "GREEDY" , 6) == 0 || strncmp(inst->algorithm , "GREEDY_2_OPT" , 12) == 0){

		inst->prob = 0.0;
		Grasp(inst , start_node);

	}else if(strncmp(inst->algorithm , "GRASP" , 5) == 0 || strncmp(inst->algorithm , "GRASP_2_OPT" , 11) == 0){

		Grasp(inst , start_node);

	}else if(strncmp(inst->algorithm , "EXTRAMILL" , 9) == 0 || strncmp(inst->algorithm , "EXTRAMILL_2_OPT" , 15) == 0){

		Extra_Mileage(inst);

	}else if(strncmp(inst->algorithm , "VNS" , 3) == 0 ){

		vns(inst , start_node);

	}else if(strncmp(inst->algorithm , "LIN_TABU" , 8) == 0 || strncmp(inst->algorithm , "STEP_TABU" , 9) == 0 || strncmp(inst->algorithm , "RAND_TABU" , 9) == 0){

		tabu( inst , start_node );

	}else if(strncmp(inst->algorithm , "ANNEALING" , 9) == 0 ){

		annealing( inst, start_node );

	}else if(strncmp(inst->algorithm , "GENETIC" , 7) == 0 ){

		genetic_algo( inst );

	}else if(strncmp(inst->algorithm , "BENDERS" , 7) == 0 || strncmp(inst->algorithm , "CALLBACK" , 8) == 0 || strncmp(inst->algorithm , "USER_CALLBACK" , 13) == 0 ){

		benders(inst);

	}else if(strncmp(inst->algorithm , "HARDFIX" , 7) == 0 || strncmp(inst->algorithm , "HARDFIX_CALL" , 12) == 0){

		hardfixing(inst);

	}else if(strncmp(inst->algorithm , "BRANCHING" , 9) == 0 || strncmp(inst->algorithm , "BRANCHING_CALL" , 14) == 0){

		branching(inst);

	}else{
		print_error(" Algorithm not valid \n");
	}

}

void read_input(instance *inst) {
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");
	
	inst->nnodes = -1;

    // Open file
    FILE *fp = fopen(inst->input_file, "r");
   
    char line[128];          // 1 line of the file
    char *par_name;          // name of the parameter in the readed line
    char *token1;            // value of the parameter in the readed line
    char *token2;            // second value of the parameter in the readed line (used for reading coordinates)
    int active_section = 0;  // 0=reading parameters, 1=NODE_COORD_SECTION, 2=EDGE_WEIGHT_SECTION
    char sep[] = " :\n\t\r"; // separator for parsing

    // Read the file line by line
    while(fgets(line, sizeof(line), fp) != NULL) {
        if (strlen(line) <= 1 ) continue; // skip empty lines
        par_name = strtok(line, sep);

        if(strncmp(par_name, "NAME", 4) == 0){
			active_section = 0;
            //future use
			continue;
		}

		if(strncmp(par_name, "COMMENT", 7) == 0){
			active_section = 0;    
			continue;
		}   

        if(strncmp(par_name, "TYPE", 4) == 0){
            active_section = 0;
            continue;
		}

        if(strncmp(par_name, "DIMENSION", 9) == 0 ){
            token1 = strtok(NULL, sep);
            inst->nnodes = atoi(token1);
            inst->nodes = calloc(inst->nnodes, sizeof(point));
			inst->solution.edges = calloc((inst->nnodes) , sizeof(edge));
            active_section = 0;  
            continue;
		}

        if(strncmp(par_name, "EOF", 3) == 0 ){
			active_section = 0;
			break;
		}

        if(strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ){
			active_section = 0;  
            //future implement
            continue;
		}

        if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0){
            active_section = 1;
            continue;
        }

        if (strncmp(par_name, "EDGE_WEIGHT_SECTION", 19) == 0 || strncmp(par_name , "TOUR_SECTION" , 12) == 0){
            active_section = 2;
            continue;
        }

        // NODE_COORD_SECTION
        if(active_section == 1){ 
            int i = atoi(par_name) - 1; // Nodes in problem's file start from index 1
			if ( i < 0 || i >= inst->nnodes) print_error(" unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, sep);
			token2 = strtok(NULL, sep);
			point p = {atoi(token1), atoi(token2)};
			inst->nodes[i] = p;
            continue;
        }
        
        
    }

	fclose(fin);    
	
}

void parse_command_line(int argc, char** argv, instance *inst , int *start_node , char* file_name) 
{ 
		
	// default   
	strcpy(inst->input_file, "NULL");
	inst->integer_cost = 0;
	inst->randomseed = 0; 
	inst->timelimit = 0;
	inst->solution.obj_best = INFINITY;
	//inst->algorithm = GRASP;
	//set default value for the probability into grasp and greedy for choose the best new node or the second 
	inst->prob = 0.0;
	//in default 2opt refinement id deactivated
	inst->active2opt = 0;
	//dafault verbose set to debug
	inst->verbose = 0;
	strcpy( file_name, "results/" ); 

	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); strcat( file_name , argv[i] ); continue; } 			// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atoi(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-n") == 0 ) { inst->nnodes = atoi(argv[++i]); strcat( file_name , argv[i] ); continue; }                  //number of nodes
		//if ( strcmp(argv[i],"-l") == 0 ) { inst->path_length = atoi(argv[++i]); continue; } 			// path length of sol
		//if ( strcmp(argv[i],"-d") == 0 ) { inst->debug = 1; continue; }
		if ( strcmp(argv[i],"-start_n") == 0 ) { *start_node = atoi(argv[++i]); continue; } 	   // start node
		//if ( strcmp(argv[i],"-end_n") == 0 ) { inst->end_node = atoi(argv[++i]); continue; }         // end node
		//type of algorithm input
		if (strcmp(argv[i] , "-algo" ) == 0){ strcpy(inst->algorithm,argv[++i]); strcat( file_name , argv[i] ); strcat( file_name , "_" );continue; }
		if (strcmp(argv[i] , "-prob" ) == 0){inst->prob = atof(argv[++i]); continue; }
		if (strcmp(argv[i] , "-verbose" ) == 0){inst->verbose = atoi(argv[++i]); continue; }

    }    

	strcat( file_name, ".txt" );  

}    
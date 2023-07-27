#include "tsp.h"

void point_file(instance *inst);
void solution_file(instance *inst , int start_n);
void gen_plot();

void point_file(instance *inst){
	FILE *fin = fopen("../plot/points.dat", "w");
	
	//generate file for gnuplot
	for(int i=0 ; i<inst->nnodes ; i++){
		fprintf(fin , "%f %f\n" , inst->nodes[i].x , inst->nodes[i].y);
	}
	
	fclose(fin);

}

void solution_file(instance *inst , int start_n){

	FILE *sol = fopen("../plot/solution.dat", "w");
	//generate file that contain the solution sequence of points
	int curr_node = start_n;
	while(1){
		fprintf(sol , "%f %f\n" , inst->nodes[inst->solution.edges[curr_node].i].x , inst->nodes[inst->solution.edges[curr_node].i].y);	
		curr_node = inst->solution.edges[curr_node].j;
		if(curr_node == start_n){
			fprintf(sol , "%f %f\n" , inst->nodes[inst->solution.edges[curr_node].i].x , inst->nodes[inst->solution.edges[curr_node].i].y);	
			break;
		}
	}
	fclose(sol);

}

void gen_plot(){
	system("gnuplot ../plot/commands.txt");
}
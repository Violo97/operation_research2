#ifndef DEBUG_  

#define DEBUG_

 
#include "tsp.h"

void control_obj(instance *inst , int start_node);
void visit_all_node(instance *inst , int start_node);
void control_cplex_sol(const double *xstar, instance *inst);

#endif 
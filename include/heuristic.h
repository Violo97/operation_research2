#ifndef HEURISTIC_  

#define HEURISTIC_

 
#include "tsp.h"

void grasp_algorithm(instance *inst , int start_n);
void extramileage_algorithm(instance *inst);
void twoopt_algorithm(instance *inst);
void Grasp( instance *inst , int start_n);
void Extra_Mileage( instance *inst );


#endif 
#ifndef CALLBACKS_  

#define CALLBACKS_
 
#include "utility.h"
#include "tsp.h" 
#include <cplex.h>

int CPXPUBLIC cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle );



#endif 
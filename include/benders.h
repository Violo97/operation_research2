#ifndef BENDERS_  

#define BENDERS_

#include <cplex.h> 
#include <pthread.h>
#include "utility.h"
#include "tsp.h" 

extern pthread_mutex_t mutex;

int TSPopt(instance *inst);
int benders(instance *inst);
int benders_lu(instance *inst);
int benders_patching(instance *inst);
int benders_callbaks(instance *inst);
void constr_gen( CPXENVptr env, CPXLPptr lp ,CPXCALLBACKCONTEXTptr context, double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void heur_patch(instance *inst, int *succ, int *comp, int * , double *xheu);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void updateSol(instance *inst , int *succ , double *xheu , int update);
void benders_solver(CPXENVptr env , CPXLPptr lp ,  instance *inst);
//utility for cplex
void setParam(instance *inst , CPXENVptr env);
void closeCplex(CPXENVptr env, CPXLPptr lp);
void initSol(instance *inst , CPXENVptr env, CPXLPptr lp );


#endif 
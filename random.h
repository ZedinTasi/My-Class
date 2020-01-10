#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <stdlib.h>
#include <time.h>
#include <math.h>

double random_uniform();
double random_uniform(double a, double b);

double random_normal();
double random_normal(double mean, double var);





#endif //_RANDOM_H_
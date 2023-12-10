#include <R.h>
#include <Rmath.h>
#include "wquantile.h"
#include "constants.h"

#ifndef _HUBER2_H
#define _HUBER2_H
void whuber2(double* restrict, double* restrict, double* restrict, double*,
    double*, double*, int*, int*, const double*, int*, int*);
#endif

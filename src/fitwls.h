#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#ifndef _FITWLS_H
#define _FITWLS_H

// prototypes for the functions
void fitwls(double* restrict, double* restrict, double* restrict,
    double* restrict, double* restrict, double* restrict, double* restrict,
    int*, int*, double*, int*, int*);
#endif

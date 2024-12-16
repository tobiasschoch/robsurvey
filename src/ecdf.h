#include <R.h>
#include <Rmath.h>

#ifndef _ECDF_H
#define _ECDF_H

void ecdf_cd(double* restrict res_sample, double* restrict linpred_nonsample,
    double* restrict sd_nonsample, double* restrict at, int *at_length, int *n,
    int *N_minus_n, double *result);
#endif

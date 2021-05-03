#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h"
#include "psifunctions.h"

// #include "mallows.h"

#ifndef _ROBSURVEY_H
#define _ROBSURVEY_H

// prototypes for the functions
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*,
    double*, double*, double*, double*, int*, int*, int*);
void cov_rwlslm(double*, double*, double*, double*, double*, double*, double*,
    double*, int*, int*, int*, int*);
#endif

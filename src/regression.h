// define before any R headers
#ifndef USE_FC_LEN_T
    #define USE_FC_LEN_T
#endif
// include Rconfig.h (possibly via R.h)
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
    #define FCONE
#endif

#include <Rmath.h>
#include "robsurvey_error.h"
#include "wquantile.h"
#include "psifunctions.h"
#include "mallows.h"
#include "regression_scale.h"
#include "regression_data.h"

#ifndef _REGRESSION_H
#define _REGRESSION_H

// prototypes for the functions
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*,
    double*, double*, double*, double*, int*, int*, int*, int*, int*, int*,
    int*);
void wlslm(double*, double*, double*, double*, int*, int*, double*, double*);
#endif

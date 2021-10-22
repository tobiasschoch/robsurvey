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
#include "constants.h"

#ifndef _ROBSURVEY_H
#define _ROBSURVEY_H

// prototypes for the functions
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*,
    double*, double*, double*, double*, int*, int*, int*, int*, int*);
void cov_reg_model(double*, double*, double*, double*, double*, double*,
    double*, double*, int*, int*, int*, int*, int*);
void cov_reg_design(double*, double*, double*, double*, double*, double*, int*,
    int*, int*, int*, int*, double*);
#endif

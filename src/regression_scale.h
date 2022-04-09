#include <Rmath.h>
#include "robsurvey_error.h"
#include "wquantile.h"
#include "constants.h"
#include "regression_data.h"

#ifndef _REGRESSION_SCALE_H
#define _REGRESSION_SCALE_H

// prototypes for the functions
robsurvey_error_type wmad(regdata*, workarray*, double* restrict, int*,
    double, double*);
robsurvey_error_type wiqr(regdata*, workarray*, double* restrict, double*);
#endif

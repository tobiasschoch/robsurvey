/******************************************************************************\
|* robsurvey                                                                  *|
|* -------------------------------------------------------------------------- *|
|* PROJECT  robsurvey library                                                 *|
|* SUBEJCT  header file for basic statistics functions                        *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January, 2020              *|
|* LICENSE  GPL >= 2                                                          *|
|* COMMENT  [none]                                                            *|
\******************************************************************************/
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h"
#include "mallows.h"
#include "huberm.h"
#include "fitwls.h"

#ifndef _ROBSURVEY_H
#define _ROBSURVEY_H

// prototypes for the functions
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*,
    double*, double*, double*, double*, int*, int*, int*);
void cov_rwlslm(double*, double*, double*, double*, double*, double*, double*,
    double*, int*, int*, int*, int*);
#endif

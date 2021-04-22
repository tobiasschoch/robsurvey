/******************************************************************************\
|* fitwls                                                                     *|
|* -------------------------------------------------------------------------- *|
|* PROJECT  robsurvey library                                                 *|
|* SUBEJCT  header file for weighted least square                             *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), Dec 6, 2020                *|
|* LICENSE  GPL >= 2                                                          *|
|* COMMENT  [none]                                                            *|
\******************************************************************************/
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#ifndef _FITWLS_H
#define _FITWLS_H

// prototypes for the functions
void fitwls(double*, double*, double*, double*, double*, double*, double*, int*,
    int*, double*, int*, int*);
#endif

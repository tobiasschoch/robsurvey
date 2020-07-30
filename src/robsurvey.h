/*****************************************************************************\
|* robsurvey								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  robsurvey library						     *|
|* SUBEJCT  header file for basic statistics functions	 		     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h" 

#ifndef _ROBSURVEY_H 
#define _ROBSURVEY_H
// prototypes for the functions
void wtrimmedmean(double*, double*, double*, double*, double*, int*);
void wwinsorizedmean(double*, double*, double*, double*, double*, int*);
void huberm(double*, double*, double*, double*, double*, double*, int*, 
   int*, const double*);
void rwlslm(double*, double*, double*, double*, double*, int*, int*, double*, 
   double*, double*, int*, double*, int*, double*, double*);
void wkwinsorizedmean(double*, double*, int*, double*, int*, double*);
#endif 

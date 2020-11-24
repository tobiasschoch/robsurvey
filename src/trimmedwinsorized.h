/*****************************************************************************\
|* robsurvey								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  robsurvey library						     *|
|* SUBEJCT  header file for trimmed and winsorized mean	 		     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/

#include "wquantile.h"

#ifndef _TRIMMEDWINSORIZED_H 
#define _TRIMMEDWINSORIZED_H 

// prototypes for the functions
void wtrimmedmean(double*, double*, double*, double*, double*, int*);
void wwinsorizedmean(double*, double*, double*, double*, double*, int*);
void wkwinsorizedmean(double*, double*, int*, double*, int*, double*);
#endif 

/*****************************************************************************\
|* wquantile: weighted quantiles 					     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sctbase 							     *|
|* SUBEJCT  header file for weighted quantile				     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), February 10, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include <R.h>

#ifndef _WQUANTILE_H
#define _WQUANTILE_H

// prototypes for the functions
void wquantile(double*, double*, int*, double*, double*);

#endif 

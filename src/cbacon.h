/*****************************************************************************\
|* bacon								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sct library							     *|
|* SUBEJCT  header file for weighted BACON				     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h" 

#ifndef _CBACON_H 
#define _CBACON_H

void wbacon(double*, double*, double*, double*, double*, int*, int*, double*, 
   int*, double*, int*, int*);

#endif 

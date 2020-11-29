#include <R.h>

#ifndef _WQUANTILE_H
#define _WQUANTILE_H
void wselect0(double*, double*, int, int, int);
void wquantile(double*, double*, int*, double*, double*);
void wquantile_noalloc(double*, double*, double*, int*, double*, double*);
#endif

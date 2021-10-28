#include <R.h>
#include <Rmath.h>

#ifndef _PSIFUNCTIONS_H
#define _PSIFUNCTIONS_H

// return type for wgt, psi, and psi-prime functions
typedef double (*f_ptr)(double, const double);

// prototypes for functions whose return value is a function pointer
f_ptr get_wgt_function(int);
f_ptr get_psi_function(int);
f_ptr get_psi_prime_function(int);

// prototypes for function callable from R
void psi_function(double*, double*, int*, int*, double*);
#endif

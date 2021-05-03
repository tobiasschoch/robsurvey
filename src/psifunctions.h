#include <R.h>
#include <Rmath.h>

#ifndef _PSIFUNCTIONS_H
#define _PSIFUNCTIONS_H

// prototypes for the functions
double huber_wgt(double, const double);
double huber_psi(double, double);
double huber_psi_prime(double, const double);

double huber_wgt_asym(double, const double);
double huber_psi_asym(double, const double);
double huber_psi_prime_asym(double, const double);

double tukey_wgt(double, const double);
double tukey_psi(double, const double);
double tukey_psi_prime(double, const double);
#endif

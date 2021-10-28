/* Psi-, weight- and psi-prime functions for M-estimators

   Copyright (C) 2020-21 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, a copy is available at
   https://www.gnu.org/licenses/

   Notes:
    Currently, the following psi-functions (and their respective weight-
    and psi-prime functions) are implemented:
       - psi = 0: Huber psi-function (see Huber, 1981)
       - psi = 1: asymmetric Huber psi-function (see Hulliger, 1995)
       - psi = 2: Tukey biweight psi-function (see Huber, 1981)
    If you want to add other psi-functions, you must include them in the
    switch case statements of the functions 'get_wgt_function',
    'get_psi_function', and 'get_psi_prime_function'; see below.

  References:
    - Huber, P.J. (1981). Robust Statistics, New York: John Wiley & Sons.
    - Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimator, Survey
      Methodology 21, pp. 79-87.
*/

#include "psifunctions.h"
#define _POWER2(_x) ((_x) * (_x))

double huber_wgt(double, const double);
double huber_psi(double, double);
double huber_psi_prime(double, const double);

double huber_wgt_asym(double, const double);
double huber_psi_asym(double, const double);
double huber_psi_prime_asym(double, const double);

double tukey_wgt(double, const double);
double tukey_psi(double, const double);
double tukey_psi_prime(double, const double);

/******************************************************************************\
|* Obtain the weight (wgt), psi- or psi-prime function associated with a      *|
|* psi-function; this wrapper function has its own return type 'f_ptr'        *|
\******************************************************************************/
f_ptr get_wgt_function(int psi)
{
    switch (psi) {
    case 0: // weight of Huber psi-function
        return huber_wgt;
    case 1: // weight of Huber asymmetric psi-function
        return huber_wgt_asym;
    case 2: // weight of Tukey biweight psi-function
        return tukey_wgt;
    default:
        return huber_wgt;
    }
}
f_ptr get_psi_function(int psi)
{
    switch (psi) {
    case 0: // Huber psi-function
        return huber_psi;
    case 1: // Huber asymmetric psi-function
        return huber_psi_asym;
    case 2: // Tukey biweight psi-function
        return tukey_psi;
    default:
        return huber_psi;
    }
}
f_ptr get_psi_prime_function(int psi)
{
    switch (psi) {
    case 0: // Huber psi-function
        return huber_psi_prime;
    case 1: // Huber asymmetric psi-function
        return huber_psi_prime_asym;
    case 2: // Tukey biweight psi-function
        return tukey_psi_prime;
    default:
        return huber_psi_prime;
    }
}

/******************************************************************************\
|* psi-function callable from R                                               *|
|*                                                                            *|
|* x     vector, array[n]                                                     *|
|* k     tuning constant                                                      *|
|* n     dimension                                                            *|
|* psi   0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight                  *|
|* res   on return: psi(x), array[n]                                          *|
\******************************************************************************/
void psi_function(double *x, double *k, int *n, int *psi, double *res)
{
    double (*f_psi)(double, const double);
    f_psi = get_psi_function(*psi);
    for (int i = 0; i < *n; i++)
        res[i] = (*f_psi)(x[i], *k);
}

/******************************************************************************\
|* Huber psi-, psi-prime- and wgt-functions                                   *|
\******************************************************************************/
double huber_psi(double x, const double k)
{
    return (x <= -k) ? -k : ((x < k) ? x : k);
}

double huber_psi_prime(double x, const double k)
{
    return (fabs(x) <= k) ? 1.0 : 0.0;
}

double huber_wgt(double x, const double k)
{
    double z = fabs(x);
    if (z > DBL_EPSILON)
        return (z >= k) ? k / z : 1.0;
    else
        return 0.0;
}

/******************************************************************************\
|* Huber asymmetric psi-, psi-prime- and wgt-functions                        *|
\******************************************************************************/
double huber_psi_asym(double x, const double k)
{
    return (x <= k) ? x : k;
}

double huber_psi_prime_asym(double x, const double k)
{
    return (x <= k) ? 1.0 : 0.0;
}

double huber_wgt_asym(double x, const double k)
{
    double z = fabs(x);
    if (z > DBL_EPSILON)
        return (x <= k) ? 1.0 : k / x;
    else
        return 0.0;
}

/******************************************************************************\
|* Tukey biweigt psi-, psi-prime- and wgt-functions                           *|
\******************************************************************************/
double tukey_psi(double x, const double k)
{
    if (fabs(x) > k) {
        return 0.0;
    } else {
        double z = x / k;
        double u = 1.0 - _POWER2(z);
        return x * _POWER2(u);
    }
}

double tukey_psi_prime(double x, const double k)
{
    if (fabs(x) > k) {
        return 0.0;
    } else {
        x /= k;
        double z = _POWER2(x);
        return (1.0 - z) * (1.0 - 5.0 * z);
    }
}

double tukey_wgt(double x, const double k)
{
    if (fabs(x) > k) {
        return 0.0;
    } else {
        double z = x / k;
        z = (1.0 - z) * (1.0 + z);
        return _POWER2(z);
    }
}
#undef _POWER2

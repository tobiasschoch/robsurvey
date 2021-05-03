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
*/

#include "psifunctions.h"
#define _POWER2(_x) ((_x) * (_x))

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

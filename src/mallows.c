/* Functions used to compute the Mallows generalized regression M-estimator

   Copyright (C) 2020-2026 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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

#include "mallows.h"

// data structure used in the utility function zeroin_mallows_mad
typedef struct zeroin_data_struct {
    int n;          // sample size
    double *xwgt;   // ptr to the weights in the design space, xwgt
} zeroin_data;

// declaration
double zeroin_mallows_mad(double, void*);

/******************************************************************************\
|* Compute normalization constant (wrapper function for Richard Brent's root  *|
|* finding algorithm 'zeroin')                                                *|
|*                                                                            *|
|*  xwgt     weights in the design space, array[n]                            *|
|*  n        dimension                                                        *|
|*  constant on return: normalization constant                                *|
\******************************************************************************/
robsurvey_error_type mallows_mad_normalization(double* restrict xwgt, int *n,
    double *constant)
{
    // initialize and populate the typedef struct zeroin_data
    zeroin_data data;
    zeroin_data *dat = &data;
    dat->n = *n;
    dat->xwgt = xwgt;

    // compute the root with Richard Brent's zeroin algorithm ('zeroin.c')
    int Maxit = zeroin_MAXIT;
    double ax = zeroin_LEFTBOUND, bx = zeroin_RIGHTBOUND, Tol = zeroin_TOL;
    double fa = zeroin_mallows_mad(ax, dat);
    double fb = zeroin_mallows_mad(bx, dat);
    double result = R_zeroin2(ax, bx, fa, fb, zeroin_mallows_mad, dat, &Tol,
        &Maxit);

    if (Maxit > 0) {
        *constant = (result > DBL_EPSILON) ? 1.0 / result : mad_NORM_CONSTANT;
        return ROBSURVEY_ERROR_OK;
    } else {
        return ROBSURVEY_ERROR_MALLOWS_NOT_CONVERGED;
    }
}

/******************************************************************************\
|* Estimating equation whose root is the normalization constant of the        *|
|* weighted MAD for Mallows GM-estimator of regression.                       *|
|*                                                                            *|
|*  constant argument (root)                                                  *|
|*  ptr      void pointer (called with an instance of typedef zeroin_data)    *|
\******************************************************************************/
double zeroin_mallows_mad(double constant, void *ptr)
{
    // cast ptr to typedef zeroin_data
    zeroin_data *dat;
    dat = (zeroin_data*)ptr;
    int n = dat->n;
    double* restrict xwgt = dat->xwgt;

    double tmp = 0.0;
    for (int i = 0; i < n; i++)
        tmp += pnorm(constant, 0.0, sqrt(xwgt[i]), 1, 0);

    tmp /= (double)n;
    return tmp - 0.75;
}

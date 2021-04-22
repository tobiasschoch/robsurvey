/* Functions to compute weighted (generalized) regression M-estimators,
   winsorized, and trimmed estimators of location

   Copyright (C) 2020 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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

#include "wquantile.h"

/******************************************************************************\
|* weighted trimmed mean (scalar)                                             *|
|*                                                                            *|
|*  x     data, array[n]                                                      *|
|*  w     weights, array[n]                                                   *|
|*  lo    lower bound [0, 1]                                                  *|
|*  hi    upper bound [0, 1] s.t. lo < hi                                     *|
|*  mean  on return: weighted trimmed mean                                    *|
|*  n     dimension                                                           *|
\******************************************************************************/
void wtrimmedmean(double *x, double *w, double *lo, double *hi, double *mean,
    int *n)
{
    double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;
    double *work_2n;
    work_2n = (double*) Calloc(2 * *n, double);

    // quantiles
    wquantile_noalloc(x, w, work_2n, n, lo, &quantile_lo);
    wquantile_noalloc(x, w, work_2n, n, hi, &quantile_hi);

    // trimmed mean
    for (int i = 0; i < *n; i++) {
        if (quantile_lo <= x[i] && x[i] <= quantile_hi) {
            sum_x += x[i] * w[i];
            sum_w += w[i];
        }
    }

    if (sum_w > DBL_EPSILON) {
        *mean = sum_x / sum_w;
    } else {
        *mean = 0.0;
        error("Error: trimmed mean: division by zero\n");
    }
    Free(work_2n);
}

/******************************************************************************\
|* weighted winsorized mean (scalar)                                          *|
|*                                                                            *|
|*  x     data, array[n]                                                      *|
|*  w     weights, array[n]                                                   *|
|*  lo    lower bound [0, 1]                                                  *|
|*  hi    upper bound [0, 1] s.t. lo < hi                                     *|
|*  mean  on return: weighted trimmed mean                                    *|
|*  n     dimension                                                           *|
\******************************************************************************/
void wwinsorizedmean(double *x, double *w, double *lo, double *hi, double *mean,
    int *n)
{
    double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;
    double *work_2n;
    work_2n = (double*) Calloc(2 * *n, double);

    // quantiles
    wquantile_noalloc(x, w, work_2n, n, lo, &quantile_lo);
    wquantile_noalloc(x, w, work_2n, n, hi, &quantile_hi);

    // winsorized mean
    for (int i = 0; i < *n; i++) {
        if (x[i] < quantile_lo) {
            sum_x += quantile_lo * w[i];
        } else {
            if (x[i] < quantile_hi)
                sum_x += x[i] * w[i];
            else
                sum_x += quantile_hi * w[i];
        }
        sum_w += w[i];
    }

    *mean = sum_x / sum_w;
    Free(work_2n);
}

/******************************************************************************\
|* weighted k one-sided winsorized mean (scalar)                              *|
|*                                                                            *|
|*  x      data, array[n]                                                     *|
|*  w      weights, array[n]                                                  *|
|*  k      k-th largest element (zero index)                                  *|
|*  mean   on return: weighted trimmed mean                                   *|
|*  n      dimension                                                          *|
|*  prob   on return: estimated probability                                   *|
\******************************************************************************/
void wkwinsorizedmean(double *x, double *w, int *k, double *mean,
    int *n, double *prob)
{
    double sum_xw = 0.0, below_sum_w = 0.0, above_sum_w = 0.0;

    // determine k-th largest element (it puts element k into its final
    // position in the sorted array)
    wselect0(x, w, 0, *n - 1, *k);

    // partial sums of x[0..k] * w[0..k] and w[0..k]
    for (int i = 0; i < *k + 1; i++) {
        sum_xw += w[i] * x[i];
        below_sum_w += w[i];
    }

    // partial sum of w[(k+1)..(n-1)]
    for (int i = *k + 1; i < *n; i++)
        above_sum_w += w[i];

    sum_xw += above_sum_w * x[*k];
    *mean = sum_xw / (below_sum_w + above_sum_w);       // winsorized mean
    *prob = below_sum_w / (below_sum_w + above_sum_w);  // estimate of prob
}

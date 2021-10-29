/* Functions for winsorized and trimmed estimators of location

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

   References:
   Higham (2002). Accuracy and Stability of Numerical Algorithms, 2nd ed.,
       Philadelphia: SIAM.
*/

#include "trimmedwinsorized.h"

/******************************************************************************\
|* weighted trimmed mean (scalar)                                             *|
|*                                                                            *|
|*  x       data, array[n]                                                    *|
|*  w       weights, array[n]                                                 *|
|*  lo      lower bound [0, 1]                                                *|
|*  hi      upper bound [0, 1] s.t. lo < hi                                   *|
|*  mean    on return: weighted trimmed mean                                  *|
|*  n       dimension                                                         *|
|*  success on return: 1: successful; 0: failure (division by zero)           *|
\******************************************************************************/
void wtrimmedmean(double* restrict x, double* restrict w, double *lo,
    double *hi, double *mean, int *n, int *success)
{
    *success = 1;

    if (*n == 1) {
        *mean = x[0];
        return;
    }

    // quantiles
    double quantile_lo, quantile_hi;
    double *work_2n = (double*) Calloc(2 * *n, double);
    wquantile_noalloc(x, w, work_2n, n, lo, &quantile_lo);
    wquantile_noalloc(x, w, work_2n, n, hi, &quantile_hi);
    Free(work_2n);

    // Kahan compensated (weighted) summation; see e.g., Higham (2002, ch. 4.3)
    double sum_w = 0.0, sum_x = 0.0, comp = 0.0, tmp, a;
    for (int i = 0; i < *n; i++) {
        if (quantile_lo <= x[i] && x[i] <= quantile_hi) {
            tmp = sum_x;
            a = w[i] * x[i] + comp;
            sum_x = tmp + a;
            comp = (tmp - sum_x) + a;

            sum_w += w[i];
        }
    }
    sum_x += comp;

    if (sum_w > DBL_EPSILON) {
        *mean = sum_x / sum_w;
    } else {
        *mean = 0.0;
        *success = 0;
    }
}

/******************************************************************************\
|* weighted winsorized mean (scalar)                                          *|
|*                                                                            *|
|*  x       data, array[n]                                                    *|
|*  w       weights, array[n]                                                 *|
|*  lo      lower bound [0, 1]                                                *|
|*  hi      upper bound [0, 1] s.t. lo < hi                                   *|
|*  mean    on return: weighted trimmed mean                                  *|
|*  n       dimension                                                         *|
\******************************************************************************/
void wwinsorizedmean(double* restrict x, double* restrict w, double *lo,
    double *hi, double *mean, int *n)
{
    if (*n == 1) {
        *mean = x[0];
        return;
    }

    // quantiles
    double quantile_lo, quantile_hi;
    double *work_2n = (double*) Calloc(2 * *n, double);
    wquantile_noalloc(x, w, work_2n, n, lo, &quantile_lo);
    wquantile_noalloc(x, w, work_2n, n, hi, &quantile_hi);
    Free(work_2n);

    // Kahan compensated (weighted) winsorized summation; see e.g., Higham
    // (2002, ch. 4.3)
    double sum_w = 0.0, sum_x = 0.0, comp = 0.0, tmp, a;
    for (int i = 0; i < *n; i++) {
        tmp = sum_x;

        if (x[i] < quantile_lo) {
            a = quantile_lo * w[i];
        } else if (x[i] < quantile_hi) {
            a = x[i] * w[i];
        } else {
            a = quantile_hi * w[i];
        }
        a += comp;
        sum_x = tmp + a;
        comp = (tmp - sum_x) + a;
        sum_w += w[i];
    }

    *mean = sum_x / sum_w;
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
void wkwinsorizedmean(double* restrict x, double* restrict w, int *k,
    double *mean, int *n, double *prob)
{
    // determine k-th largest element
    *k = *n - *k - 1;
    wselect0(x, w, 0, *n - 1, *k);
    double cutoff = x[*k];

    // Kahan compensated (weighted) one-sided winsorized summation; see e.g.,
    // Higham (2002, ch. 4.3)
    double sum_x = 0.0, sum_w = 0.0, comp = 0.0, tmp, a, below_sum_w = 0.0;
    for (int i = 0; i < *n; i++) {
        tmp = sum_x;

        if (x[i] <= cutoff) {
            a = x[i] * w[i] + comp;
            below_sum_w += w[i];
        } else {
            a = cutoff * w[i] + comp;
        }

        sum_x = tmp + a;
        comp = (tmp - sum_x) + a;
        sum_w += w[i];
    }

    *mean = sum_x / sum_w;

    // weighted ecdf(x[k])
    *prob = below_sum_w / sum_w;
}

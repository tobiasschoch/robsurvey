/* Model-based and model-assisted estimators of the empirical distribution
   function.

   Copyright (C) 2024 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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

#include "ecdf.h"

// variadic arguments in macros are supported by gcc (>=3.0), clang (all
// versions), visual studio (>=2005); the version with ##__VA_ARGS__ (which
// will silently eliminate the trailing comma) is not portable (clang complains)
#define PRINT_OUT(...) Rprintf(__VA_ARGS__)

double cd_sum(double*, double*, int, int);

/******************************************************************************\
|* Chambers-Dunstan estimate of the empirical distribution function           *|
|*                                                                            *|
|*  res_sample          sample residuals, array[n]                            *|
|*  linpred_nonsample   linear predictor of non-sampled part,                 *|
|*                      array[N_minus_n]                                      *|
|*  sd_nonsample        std. dev. in non-sampled part (heterscedastic model), *|
|*                      array[N_minus_n]                                      *|
|*  at                  value at which the ECDF is evaluated                  *|
|*  at_length           number of points where ECDF should be evaluated       *|
|*  n                   sample size                                           *|
|*  N_minus_n           size of the non-sampled part                          *|
|*  on return:          estimate                                              *|
\******************************************************************************/
void ecdf_cd(double* restrict res_sample, double* restrict linpred_nonsample,
    double* restrict sd_nonsample, double* restrict at, int *at_length, int *n,
    int *N_minus_n, double *result)
{
    // allocate memory
    double* restrict res_nonsample = (double*) R_Calloc(*N_minus_n, double);
    if (res_nonsample == NULL) {
        PRINT_OUT("Error: Cannot allocate memory\n");
        return;
    }
    int* restrict index_nonsample = (int*) R_Calloc(*N_minus_n, int);
    if (index_nonsample == NULL) {
        PRINT_OUT("Error: Cannot allocate memory\n");
        return;
    }

    // standardized residuals in the non-sampled part
    for (int i = 0; i < *N_minus_n; i++)
        res_nonsample[i] = (at[0] - linpred_nonsample[i]) / sd_nonsample[i];

    // sort residuals
    R_rsort(res_sample, *n);
    rsort_with_index(res_nonsample, index_nonsample, *N_minus_n);
    // sort sd_nonsample along the index of res_nonsample
    Memcpy(linpred_nonsample, sd_nonsample, *N_minus_n);
    for (int i = 0; i < *N_minus_n; i++)
        sd_nonsample[i] = linpred_nonsample[index_nonsample[i]];

    int evals = 0;
    while (evals < *at_length) {
        // re-compute std. residuals in the non-sampled part (not need to sort)
        if (evals > 0)
            for (int i = 0; i < *N_minus_n; i++)
                res_nonsample[i] = (res_nonsample[i] * sd_nonsample[i] +
                    at[evals] - at[evals -1]) / sd_nonsample[i];

        // estimate for the sample
        int count = 0;
        for (int i = 0; i < *n; i++)
            if (res_sample[i] <= at[evals])
                count++;

        result[evals] = (double)count;

        // estimate over the non-sampled part
        result[evals] += cd_sum(res_sample, res_nonsample, *n, *N_minus_n);

        result[evals] /= (*n + *N_minus_n);
        evals++;
    }

    // housekeeping
    R_Free(res_nonsample); Free(index_nonsample);
}

/******************************************************************************\
|* Compute the double summation in the Chambers-Dunstan estimator             *|
|*                                                                            *|
|*  res_sample      sample residuals, array[n]                                *|
|*  res_nonsample   residuals of non-sampled part, array[N_minus_n]           *|
|*  n               sample size                                               *|
|*  N_minus_n       size of the non-sampled part                              *|
\******************************************************************************/
double cd_sum(double* restrict res_sample, double* restrict res_nonsample,
    const int n, const int N_minus_n)
{
    int at_sample = 0, at_nonsample = 0;
    unsigned long count = 0, cumsum = 0;

    // trivial cases
    if (res_sample[n - 1] < res_nonsample[0])
        return (double)N_minus_n;
    if (res_nonsample[N_minus_n - 1] < res_sample[0])
        return 0.0;

    // fast-forward
    while (res_nonsample[at_nonsample] < res_sample[0]) {
        at_nonsample++;
    }

    // regular counts
    for (int i = at_nonsample; i < N_minus_n; i++) {

        if (at_sample == n)
            count = n;

        while (at_sample < n && (res_sample[at_sample] <= res_nonsample[i])) {
            count++;
            at_sample++;
        }

        // cumulative sum
        cumsum += count;
    }

    // return
    return (double)cumsum / (double)n;
}
#undef PRINT_OUT

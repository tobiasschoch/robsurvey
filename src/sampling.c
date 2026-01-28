/* Sampling with probability proportional to size (pps)

   Copyright (C) 2024-2026 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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

#include "sampling.h"

static inline void set_value(double*, double, int)
    __attribute__((always_inline));

/******************************************************************************\
|* pps_prob: first-order sample inclusion probabilities for pps design        *|
|*                                                                            *|
|* size     measure of size, array[N]                                         *|
|* N        population size                                                   *|
|* n        sample size                                                       *|
|* pik      first-order inclusion probabilities, array[N]                     *|
\******************************************************************************/
void pps_prob(double *size, int *N, int *n, double *pik)
{
    double n_double = (double)*n;

    double* is_regular = (double*) R_Calloc(*N, double);
    set_value(is_regular, 1.0, *N);

    // total of size measure
    double sum_size = 0.0;
    for (int i = 0; i < *N; i++)
        sum_size += size[i];

    // degenerate case
    if (sum_size < DBL_EPSILON) {
        set_value(pik, 0.0, *N);
        return;
    }

    // ordinary case
    int exc_count, exc_count_new = 0;       // number of exceptions: pik > 1
    double exc_size, exc_size_new = 0.0;    // size associated with exception

    do {
        exc_count = exc_count_new;
        exc_count_new = 0;
        exc_size = exc_size_new;
        exc_size_new = 0.0;

        for (int i = 0; i < *N; i++) {
            pik[i] = is_regular[i] * (n_double - (double)exc_count) *
                size[i] / (sum_size - exc_size) + (1.0 - is_regular[i]);
            if (pik[i] >= 1.0) {
                exc_count_new += 1;
                exc_size_new += size[i];
                pik[i] = 1.0;
                is_regular[i] = 0.0;
            }
        }

    } while (exc_count_new != exc_count);

    // housekeeping
    R_Free(is_regular);
}

/******************************************************************************\
|* pps_brewer: draw a pps-sample (without replacement) by Brewer's method     *|
|*                                                                            *|
|* pik      first-order inclusion probabilities, array[N]                     *|
|* N        population size                                                   *|
|* n        sample size                                                       *|
|* sample   index of sampled elements, array[n]; must be zero on entry        *|
\******************************************************************************/
void pps_brewer(double *pik, int *N, int *n, int *sample)
{
    double u, p_sum, tmp;
    double n_drawn = (double)*n;

    // allocate memory for the N-vector of probabilities
    double* restrict p = (double*) R_Calloc(*N, double);
    if (p == NULL) {
        Rprintf("Error: Cannot allocate memory\n");
        return;
    }
    // initialize the R random number generator
    GetRNGstate();

    // draw elements w.p. 1
    int wp1 = 0;
    for (int i = 0; i < *N; i++) {

        if (pik[i] < 1.0)
            continue;

        // sampling w.p. 1
        sample[wp1] = i + 1;            // one-based indexing
        n_drawn -= pik[i];              // update remaining sample size
        pik[i] = 0.0;                   // not eligible for drawing
        wp1++;

        if (wp1 == *n) {
            Rprintf("Warning: All elements are sampled w.p. 1\n");
            return;
        }
    }

    // draw the rest of the sample
    for (int i = wp1; i < *n; i++) {

        // compute selection probabilities
        p[0] = pik[0] * (n_drawn - pik[0]) / (n_drawn - pik[0] *
            (double)(*n - i));
        p_sum = p[0];
        for (int j = 1; j < *N; j++) {
            tmp = pik[j] * (n_drawn - pik[j]) / (n_drawn - pik[j] *
                (double)(*n - i));
            p_sum += tmp;               // sum(p)
            p[j] = p[j - 1] + tmp;      // cumsum(p)
        }

        // select next element
        u = unif_rand() * p_sum;
        for (int j = 0; j < *N; j++) {
            if (p[j] > u) {
                sample[i] = j + 1;      // one-based indexing
                n_drawn -= pik[j];      // update remaining sample size
                pik[j] = 0.0;           // not eligible for drawing
                break;
            }
        }
    }

    // housekeeping
    PutRNGstate();
    R_Free(p);
}

/******************************************************************************\
|* set_value: set elements of array to some value                             *|
|*                                                                            *|
|* x        array[n]                                                          *|
|* value    value                                                             *|
|* n        array dimension                                                   *|
\******************************************************************************/
static inline void set_value(double *x, double value, int n)
{
    for (int i = 0; i < n; i++)
        x[i] = value;
}

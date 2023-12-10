/* Weighted Huber proposal 2 estimator of location and scale, initialized
   by, respectively, the weighted median and the weighted interquartile range

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

   References
   Higham (2002). Accuracy and Stability of Numerical Algorithms, 2nd ed.,
       Philadelphia: SIAM.
   Huber (1981). Robust Statistics, New York: John Wiley & Sons.
*/

#include <R.h>
#include "huber2.h"

// some macros
#define _WGT_HUBER(_x, _k) ((fabs(_x) >= _k) ? _k / fabs(_x) : 1.0)
#define _POWER2(_x) ((_x) * (_x))
#define _MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#define _MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))

// declaration of 'local' functions (inline imperative is GCC specific)
static inline double kappa_huber(const double) __attribute__((always_inline));

/******************************************************************************\
|* huber mean and scale estimator: proposal 2 estimator (using weights;       *|
|* the same tuning constant 'k' used for location and scale; the scale        *|
|* is normalized to be a Fisher consistent estimator at the Gaussian          *|
|* core model)                                                                *|
|*                                                                            *|
|*  x         array                                                           *|
|*  w         array of weights                                                *|
|*  robwgt    on return: array of robustness weights                          *|
|*  k         robustness tuning constant (location and scale)                 *|
|*  loc       on return: robust location                                      *|
|*  scale     on return: robust scale                                         *|
|*  n         array dimension                                                 *|
|*  maxit     maximum number of iteration (termination rule); on return:      *|
|*            effective  number of iterations                                 *|
|*  tol       numerical convergence tolerance for iterative refinement        *|
|*            iterations                                                      *|
|*  df_cor    toggle: 1: degrees of freedom (scale estimate) are n-1; 0: the  *|
|*            degrees of freedom are n (m.l.e.)                               *|
|*  success   on return: 1: successful; 0: failure (when initial scale        *|
|*            estimate = 0)                                                   *|
\******************************************************************************/
void whuber2(double* restrict x, double* restrict w, double* restrict robwgt,
    double *k, double *loc, double *scale, int *n, int *maxit,
    const double *tol, int* df_cor, int *success)
{
    int iter;
    double loc0, scale0;

    *success = 1;

    // sample size n = 1
    if (*n == 1) {
        *loc = *x;
        *scale = 0.0;
        *robwgt = 1.0;
        *maxit = 1;
        return;
    }

    // initialize location (weighted median)
    double p50 = 0.5;
    double* restrict work_2n = (double*) Calloc(2 * *n, double);
    wquantile_noalloc(x, w, work_2n, n, &p50, &loc0);

    // initialize variable for winsorized x-variable
    double* restrict x_wins = (double*) Calloc(*n, double);

    // initialize scale estimate by (weighted) interquartile range (IQR is
    // normalized to be a Fisher consistent estimate of the scale at the
    // Gaussian distr.)
    double p25 = 0.25, x25 = 0.0;
    double p75 = 0.75, x75 = 0.0;
    wquantile_noalloc(x, w, work_2n, n, &p25, &x25);
    wquantile_noalloc(x, w, work_2n, n, &p75, &x75);

    // initial scale estimate (it can be 0.0; then, the iterative updating
    // algorithm does not converge, hence fails gently)
    scale0 = (x75 - x25) * iqr_NORM_CONSTANT;

    // compute the weight total
    double wtotal = 0.0;
    for (int i = 0; i < *n; i++)
        wtotal += w[i];

    // compute consistency correction term
    double kappa = kappa_huber(*k);

    // loop
    double a, s, total_wins, ssq, tmp, comp;
    for (iter = 0; iter < *maxit; iter ++) {
        s = *k * scale0;

        // update location
        total_wins = 0.0; comp = 0.0;
        for (int i = 0; i < *n; i++) {
            // winsorized x-variable
            x_wins[i] = _MIN(_MAX(loc0 - s, x[i]), loc0 + s);
            // Kahan compensated weighted sum, see e.g., Higham (2002, ch. 4.3)
            tmp = total_wins;
            a = w[i] * x_wins[i] + comp;
            total_wins = tmp + a;
            comp = (tmp - total_wins) + a;
        }
        *loc = total_wins / wtotal;

        // update scale
        ssq = 0.0; comp = 0.0;
        for (int i = 0; i < *n; i++) {
            // Kahan compensated weighted sum
            tmp = ssq;
            a = w[i] * _POWER2(x_wins[i] - *loc) + comp;
            ssq = tmp + a;
            comp = (tmp - ssq) + a;
        }

        if (*df_cor)
            *scale = ssq / (wtotal - 1.0);
        else
            *scale = ssq / wtotal;

        // normalize the variance/ scale
        *scale = sqrt(*scale / kappa);

        // termination rule
        if (fabs(*loc - loc0) < *tol * scale0 &&
                fabs(*scale / scale0 - 1.0) < *tol) {
            break;
        } else {
            loc0 = *loc;
            scale0 = *scale;
        }
    }
    // on return: maxit = effective number of iterations
    *maxit = iter;

    // compute robustness weights
    for (int i = 0; i < *n; i++)
        robwgt[i] = _WGT_HUBER((x[i] - *loc) / *scale, *k);

    Free(x_wins); Free(work_2n);
}

/******************************************************************************\
|* This function computes the multiplicative consistency correction term      *|
|* for the Huber 'Proposal 2' type of scale estimator; Huber (1981, ch. 5.2)  *|
\******************************************************************************/
static inline double kappa_huber(const double k)
{
    if (k < 10.0) {
        double pdf_k = dnorm(k, 0.0, 1.0, 0);       // std. normal p.d.f.
        double cdf_k = pnorm(k, 0.0, 1.0, 1, 0);    // std. normal c.d.f.
        return 2.0 * (_POWER2(k) * (1.0 - cdf_k) + cdf_k - 0.5 - k * pdf_k);
    } else {
        return 1.0;
    }
}
#undef _WGT_HUBER
#undef _POWER2
#undef _MIN
#undef _MAX

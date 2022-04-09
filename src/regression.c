/* Functions to compute weighted (generalized) regression M-estimators

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

#include "regression.h"

// some macros
#define _POWER2(_x) ((_x) * (_x))

// variadic arguments in macros are supported by gcc (>=3.0), clang (all
// versions), visual studio (>=2005); the version with ##__VA_ARGS__ (which
// will silently eliminate the trailing comma) is not portable (clang complains)
#define PRINT_OUT(...) Rprintf(__VA_ARGS__)

// declaration
robsurvey_error_type compute_scale(regdata*, workarray*, double* restrict,
    int*, double, double*, int*, int*);
robsurvey_error_type rfitwls(regdata*, workarray*, double* restrict,
    double* restrict, double* restrict);
static inline double norm(const double*, const double*, const int)
    __attribute__((always_inline));
static inline void weighting_scheme(regdata*, double (*f_wgt_psi)(double,
    const double), double* restrict, double*, double*, int*, double* restrict);
robsurvey_error_type initialize(regdata*, workarray*, double* restrict,
    double* restrict, double*, int*, int*, int*, int*);

/******************************************************************************\
|* wlslm: Weighted least squares                                              *|
|*                                                                            *|
|* x          design matrix, array[n * p];                                    *|
|* y          response vector, array[n]                                       *|
|* w          sampling weight, array[n]                                       *|
|* resid      residuals, array[n]                                             *|
|* n, p       dimensions                                                      *|
|* beta0      on return: regression coefficients                              *|
|* scale      estimate of scale                                               *|
\******************************************************************************/
void wlslm(double *x, double *y, double *w, double *resid, int *n, int *p,
    double *beta0, double *scale)
{
    // initialize and populate structure with regression-specific data
    regdata data;
    regdata *dat = &data;
    dat->n = *n;
    dat->p = *p;
    dat->x = x;
    dat->y = y;
    dat->w = w;

    // initialize and populate structure with work arrays
    double* restrict work_x = (double*) Calloc(*n * *p, double);
    double* restrict work_y = (double*) Calloc(*n, double);
    workarray wwork;
    workarray *work = &wwork;
    work->work_x = work_x;
    work->work_y = work_y;

    // determine work array for 'dgels' (and allocate 'work_lapack')
    work->lwork = -1;
    robsurvey_error_type status = rfitwls(dat, work, w, beta0, resid);
    double* restrict work_lapack = (double*) Calloc(work->lwork, double);
    work->work_lapack = work_lapack;

    // compute least squares fit
    status = rfitwls(dat, work, w, beta0, resid);
    if (status != ROBSURVEY_ERROR_OK)
        PRINT_OUT("Error: %s\n", robsurvey_error(status));

    // compute scale
    double sum_w = 0.0;
    *scale = 0.0;
    for (int i = 0; i < *n; i++) {
        sum_w += w[i];
        *scale += w[i] * _POWER2(resid[i]);
    }
    *scale /= sum_w - (double)*p;
    *scale = sqrt(*scale);

    Free(work_x); Free(work_y); Free(work_lapack);
}
/******************************************************************************\
|* rwlslm: regression M- and GM-estimator and estimate of scale (w MAD)       *|
|*                                                                            *|
|* x          design matrix, array[n * p];                                    *|
|* y          response vector, array[n]                                       *|
|* w          sampling weight, array[n]                                       *|
|* resid      residuals, array[n]                                             *|
|* robwgt     robustness weight, array[n]                                     *|
|* xwgt       weights in design space, array[n]                               *|
|* n, p       dimensions                                                      *|
|* k          robustness tuning constant                                      *|
|* beta0      on return: regression coefficients                              *|
|* scale      estimate of scale (weighted MAD)                                *|
|* tol        numerical tolerance criterion to stop the iterations            *|
|* maxit      maximum number of iterations to use                             *|
|* psi        0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight             *|
|* type       0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.           *|
|* init       1 = method is initialized by weighted least squares; 0 = beta0  *|
|*            is taken as initial estimate of regression                      *|
|* mad_center 1 = mad is centered about the median, 0 = centered about zero   *|
|* verbose    1 = verbose, 0 = quiet                                          *|
|* used_iqr   on return: 1 = scale estimated by iqr, 0 = mad (default)        *|
\******************************************************************************/
void rwlslm(double *x, double *y, double *w, double *resid, double *robwgt,
    double *xwgt, int *n, int *p, double *k, double *beta0, double *scale,
    double *tol, int *maxit, int *psi, int *type, int *init, int *mad_center,
    int *verbose, int *used_iqr)
{
    // STEP 0: general preparations
    *used_iqr = 0;
    robsurvey_error_type status;
    double* restrict beta1 = (double*) Calloc(*p, double);

    // initialize and populate structure with regression-specific data
    regdata data;
    regdata *dat = &data;
    dat->n = *n;
    dat->p = *p;
    dat->x = x;
    dat->y = y;
    dat->w = w;

    // initialize and populate structure with work arrays
    double* restrict work_x = (double*) Calloc(*n * *p, double);
    double* restrict work_y = (double*) Calloc(*n, double);
    double* restrict work_2n = (double*) Calloc(2 * *n, double);
    workarray wwork;
    workarray *work = &wwork;
    work->work_x = work_x;
    work->work_y = work_y;
    work->work_2n = work_2n;

    // determine work array for 'dgels' (and allocate 'work_lapack')
    work->lwork = -1;
    status = rfitwls(dat, work, w, beta0, resid);
    double* restrict work_lapack = (double*) Calloc(work->lwork, double);
    work->work_lapack = work_lapack;

    // STEP 1: estimator type-specific preparations
    double (*f_wgt_psi)(double, const double);  // wgt-function
    f_wgt_psi = get_wgt_function(*psi);

    double mad_const = mad_NORM_CONSTANT;       // mad consistency correction

    switch (*type) {
    case 1:                                     // Mallows GM-estimator
        // consistency correction for mad
        status = mallows_mad_normalization(xwgt, n, &mad_const);
        if (status != ROBSURVEY_ERROR_OK) {
            PRINT_OUT("Error: %s\n", robsurvey_error(status));
            *maxit = 0;
            goto clean_up;
        }
        dat->xwgt = xwgt;
        break;
    case 2:                                     // Schweppe GM-estimator
        // turn xgwt into a multiplicative weight
        for (int i = 0; i < *n; i++) {
            if (fabs(xwgt[i]) < DBL_EPSILON)
                xwgt[i] = 0.0;
            else
                xwgt[i] = 1.0 / xwgt[i];
        }
        dat->xwgt = xwgt;
        break;
    }

    // STEP 2: initialize
    status = initialize(dat, work, resid, beta0, scale, init, mad_center,
        verbose, used_iqr);
    if (status != ROBSURVEY_ERROR_OK) {
        PRINT_OUT("Error: %s\n", robsurvey_error(status));
        *maxit = 0;
        goto clean_up;
    }

     // STEP 3: irwls updating
    int iterations = 0, converged = 0;
    while (iterations++ < *maxit) {
        // robustness weights: robwgt
        weighting_scheme(dat, f_wgt_psi, resid, scale, k, type, robwgt);

        // update beta and residuals
        status = rfitwls(dat, work, robwgt, beta1, resid);
        if (status != ROBSURVEY_ERROR_OK) {
            PRINT_OUT("Error: %s\n", robsurvey_error(status));
            *maxit = 0;
            goto clean_up;
        }

        // update estimate of scale
        if (*type == 1) {                       // Mallows GM
            double* restrict dummy_resid = work_x;
            for (int i = 0; i < *n; i++)
                dummy_resid[i] = resid[i] * sqrt(xwgt[i]);
            status = compute_scale(dat, work, dummy_resid, mad_center,
                mad_const, scale, verbose, used_iqr);

        } else {                                // otherwise
            status = compute_scale(dat, work, resid, mad_center, mad_const,
                scale, verbose, used_iqr);
        }

        if (status != ROBSURVEY_ERROR_OK) {
            PRINT_OUT("Error: %s\n", robsurvey_error(status));
            *maxit = 0;
            goto clean_up;
        }

        // check for convergence
        converged = (norm(beta0, beta1, *p) < *tol * *scale) ? 1: 0;
        if (converged)
            break;

        // prepare the next while run
        Memcpy(beta0, beta1, *p);
    }
    *maxit = (converged) ? iterations : 0;

    // final robustness weights
    for (int i = 0; i < *n; i++)
        robwgt[i] /= w[i];

clean_up:
    Free(beta1); Free(work_x); Free(work_y); Free(work_2n); Free(work_lapack);
}

/******************************************************************************\
|* initialization of the regression estimator (beta0 and scale)               *|
|*                                                                            *|
|* data       typedef struct regdata                                          *|
|* work       typedef struct workarray                                        *|
|* resid      on return: residuals, array[n]                                  *|
|* beta0      on return: weighted least squares estimates, array[p]           *|
|* scale      on return: weighted mad                                         *|
|* init       type of initialization                                          *|
|* mad_center mad centered about the median (1) or about zero (0)             *|
|* used_iqr   on return: 1 = scale estimated by iqr, 0 = mad (default)        *|
\******************************************************************************/
robsurvey_error_type initialize(regdata *dat, workarray *work,
    double* restrict resid, double* restrict beta0, double *scale, int *init,
    int *mad_center, int *verbose, int *used_iqr)
{
    robsurvey_error_type status;
    if (*init) {
        // compute least squares estimate of 'beta'
        status = rfitwls(dat, work, dat->w, beta0, resid);
        if (status != ROBSURVEY_ERROR_OK)
            return status;
    }
    // compute the residuals (BLAS::dgemv)
    const int n = dat->n;
    const int p = dat->p;
    const int int_1 = 1;
    const double double_minus1 = -1.0, double_1 = 1.0;
    Memcpy(resid, dat->y, n);
    F77_CALL(dgemv)("N", &n, &p, &double_minus1, dat->x, &n, beta0, &int_1,
        &double_1, resid, &int_1 FCONE);

    // compute 'scale' by the weighted mad (or iqr)
    status = compute_scale(dat, work, resid, mad_center, mad_NORM_CONSTANT,
        scale, verbose, used_iqr);

    return status;
}

/******************************************************************************\
|* Weighting scheme for iteratively reweighted least squares; weighting with  *|
|* respect to sampling, residual and design space outlyingness                *|
|*                                                                            *|
|* dat       typedef struct regdata                                           *|
|* f_wgt_psi function pointer to the wgt/psi-function                         *|
|* resid     residuals, array[n]                                              *|
|* scale     scale                                                            *|
|* k         tuning constant of the psi-function                              *|
|* type      type of estimator                                                *|
|* robwgt    on return: combined robustness (and sampling) weight             *|
\******************************************************************************/
static inline void weighting_scheme(regdata *dat,
    double (*f_wgt_psi)(double, const double), double* restrict resid,
    double *scale, double *k, int *type, double* restrict robwgt)
{
    const int n = dat->n;
    double dummy_resid;
    double* restrict w = dat->w;
    double* restrict xwgt = dat->xwgt;

    switch (*type) {
    case 0: // M-estimator
        for (int i = 0; i < n; i++)
            robwgt[i] = w[i] * f_wgt_psi(resid[i] / *scale, *k);
        break;
    case 1: // Mallows GM-estimator
        for (int i = 0; i < n; i++)
            robwgt[i] = w[i] * xwgt[i] * f_wgt_psi(resid[i] / *scale, *k);
        break;
    case 2: // Schweppe GM-estimator
        for (int i = 0; i < n; i++) {
            dummy_resid = resid[i] * xwgt[i];
            robwgt[i] = w[i] * f_wgt_psi(dummy_resid / *scale, *k);
        }
        break;
    }
}

/******************************************************************************\
|* Weighted least squares estimate (coefficients and residuals)               *|
|*                                                                            *|
|* dat    typedef struct regdata                                              *|
|* work   typedef struct work array                                           *|
|* w      weights vector, array[n]                                            *|
|* beta0  on return: coefficient vector, array[p]                             *|
|* resid  on return: residuals vector, array[n]                               *|
\******************************************************************************/
robsurvey_error_type rfitwls(regdata *dat, workarray *work, double* restrict w,
    double* restrict beta, double* restrict resid)
{
    const int n = dat->n;
    const int p = dat->p;
    int lwork = work->lwork;
    double* restrict x = dat->x;
    double* restrict y = dat->y;

    // define constants for the call of 'dgels'
    const int int_1 = 1;
    int info_dgels = 1;

    // STEP 0: determine the optimal size of array 'work_dgels'
    if (lwork < 0) {
        double dummy_work_dgels;
        F77_CALL(dgels)("N", &n, &p, &int_1, x, &n, y, &n, &dummy_work_dgels,
            &lwork, &info_dgels FCONE);
        work->lwork = (int) dummy_work_dgels;
        return ROBSURVEY_ERROR_OK;

	// STEP 1: compute least squares fit
	} else {
        double* restrict work_x = work->work_x;
        double* restrict work_y = work->work_y;
        double* restrict work_dgels = work->work_lapack;

        // pre-multiply the design matrix and the response vector by sqrt(w)
        double tmp;
        for (int i = 0; i < n; i++) {
            tmp = sqrt(w[i]);
            work_y[i] = y[i] * tmp;

            for (int j = 0; j < p; j++)
                work_x[n * j + i] = x[n * j + i] * tmp;
        }

        // compute the (weighted) least squares estimate (LAPACK::dgels),
        // solves minimize |B - A*X| for X (using QR factorization)
        F77_CALL(dgels)("N", &n, &p, &int_1, work_x, &n, work_y, &n,
            work_dgels, &lwork, &info_dgels FCONE);

        // dgels is not well suited as a rank-revealing procedure; i.e., INFO<0
        // iff a diagonal element of the R matrix is exactly 0. This is not
        // helpful; hence, we check the diagonal elements of R separately and
        // issue and error flag if any(abs(diag(R))) is close to zero
        for (int i = 0; i < p; i++)
            if (fabs(work_x[(n + 1) * i]) < sqrt(DBL_EPSILON))
                return ROBSURVEY_ERROR_RANK_DEFICIENT;

        // retrieve 'betacoefficients'
        Memcpy(beta, work_y, p);

        // compute the residuals (BLAS::dgemv)
        const double double_minus1 = -1.0, double_1 = 1.0;
        Memcpy(resid, y, n);
        F77_CALL(dgemv)("N", &n, &p, &double_minus1, x, &n, beta, &int_1,
            &double_1, resid, &int_1 FCONE);
        return ROBSURVEY_ERROR_OK;
    }
}
/******************************************************************************\
|* Scale estimator (mad or IQR)                                               *|
|*                                                                            *|
|* dat          typedef struct regdata                                        *|
|* work         typedef struct workarray                                      *|
|* resid        residuals, array[n]                                           *|
|* mad_center   1 = center mad about median, 0 = center about zero            *|
|* mad_constant normalization constant of the mad                             *|
|* scale        on return: scale estimat                                      *|
|* verbose      1 = verbose, 0 = quiet                                        *|
|* used_iqr     on return: 1 = scale estimated by iqr, 0 = mad (default)      *|
\******************************************************************************/
robsurvey_error_type compute_scale(regdata *dat, workarray *work,
    double* restrict resid, int *mad_center, double mad_constant,
    double *scale, int *verbose, int *used_iqr)
{
    // compute mad
    robsurvey_error_type status = wmad(dat, work, resid, mad_center,
        mad_constant, scale);
    // try IQR if mad is zero
    if (status == ROBSURVEY_ERROR_SCALE_ZERO) {
        if (*verbose)
            PRINT_OUT("\nNote: Scale is computed by IQR because MAD is zero\n");
        status = wiqr(dat, work, resid, scale);
        *used_iqr = 1;
        // make sure that the message is printed only once
        *verbose = 0;
    }
    return(status);
}

/******************************************************************************\
|* euclidean norm                                                             *|
|*                                                                            *|
|* x   array[p]                                                               *|
|* y   array[p]                                                               *|
|* p   dimension                                                              *|
|* NOTE: the algorithm follows S. Hammarling's implementaion of LAPACK:dnorm2 *|
\******************************************************************************/
static inline double norm(const double *x, const double *y, const int p)
{
    double abs, scale = 0.0, ssq = 1.0;
    for (int i = 0; i < p; i++) {
        abs = fabs(x[i] - y[i]);
        if (abs < DBL_EPSILON)
            continue;
        if (scale <= abs) {
            ssq = 1.0 + ssq * _POWER2(scale / abs);
            scale = abs;
        } else {
            ssq += _POWER2(abs / scale);
        }
    }
    return scale * sqrt(ssq);
}

#undef _POWER2
#undef PRINT_OUT

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

#include "robsurvey.h"

// some macros
#define _POWER2(_x) ((_x) * (_x))
#define PRINT_OUT(_f, ...) Rprintf((_f), ##__VA_ARGS__)

// structure: regression data
typedef struct regdata_struct {
    int n;
    int p;
    double *x;
    double *y;
    double *w;
    double *xwgt;
} regdata;

// structure: work arrays
typedef struct workarray_struct {
    int lwork;
    double *work_lapack;
    double *work_x;
    double *work_y;
    double *work_2n;
} workarray;

//#if 0
robsurvey_error_type cov_m_est(workarray *work, double *resid, double *x, double *robwgt,
    double *w, double *k, double *scale, double *scale2, int *n, int *p,
    double (*f_psiprime)(double, const double));
//#endif


// declaration
robsurvey_error_type inverse_qr(workarray*, double*, int*, int*, int);
robsurvey_error_type wmad(regdata*, workarray*, double* restrict, int*,
    double, double*);
robsurvey_error_type rfitwls(regdata*, workarray*, double* restrict,
    double* restrict, double* restrict);
static inline double norm(const double*, const double*, const int)
    __attribute__((always_inline));
static inline void weighting_scheme(regdata*, double (*f_wgt_psi)(double,
    const double), double* restrict, double*, double*, int*, int*,
    double* restrict);
robsurvey_error_type initialize(regdata*, workarray*, double* restrict,
    double* restrict, double*, int*, int*);

//FIXME: head
/******************************************************************************\
|*                                                                            *|
\******************************************************************************/
void rwlslm(double *x, double *y, double *w, double *resid, double *robwgt,
    double *xwgt, int *n, int *p, double *k, double *beta0, double *scale,
    double *tol, int *maxit, int *psi, int *type, int *init, int *mad_center)
{
    // STEP 0: general preparations
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
    double (*f_wgt_psi)(double, double);

    switch (*psi) {
    case 0: // weight of Huber psi-function
        f_wgt_psi = huber_wgt;
        break;
    case 1: // weight of Huber asymmetric psi-function
        f_wgt_psi = huber_wgt_asym;
        break;
    case 2: // weight of Tukey biweight psi-function
        f_wgt_psi = tukey_wgt;
        break;
    default:
        f_wgt_psi = huber_wgt;
    }

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
    status = initialize(dat, work, resid, beta0, scale, init, mad_center);
    if (status != ROBSURVEY_ERROR_OK) {
        PRINT_OUT("Error: %s\n", robsurvey_error(status));
        *maxit = 0;
        goto clean_up;
    }

     // STEP 3: irwls updating
    int iterations = 0, converged = 0;
    while (iterations++ < *maxit) {
        // robustness weights: robwgt
        weighting_scheme(dat, f_wgt_psi, resid, scale, k, psi, type, robwgt);

        // update beta and residuals
        status = rfitwls(dat, work, robwgt, beta1, resid);
        if (status != ROBSURVEY_ERROR_OK) {
            PRINT_OUT("Error: %s\n", robsurvey_error(status));
            *maxit = 0;
            goto clean_up;
        }

        // update estimate of scale
        if (*type == 1) {                       // Mallows GM
            double* restrict dummy_resid = robwgt;
            for (int i = 0; i < *n; i++)
                dummy_resid[i] = resid[i] * sqrt(xwgt[i]);
            status = wmad(dat, work, dummy_resid, mad_center, mad_const,
                scale);

        } else {                                // otherwise
            status = wmad(dat, work, resid, mad_center, mad_const, scale);
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

    // robustness weights
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
\******************************************************************************/
robsurvey_error_type initialize(regdata *dat, workarray *work,
    double* restrict resid, double* restrict beta0, double *scale, int *init,
    int *mad_center)
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
        &double_1, resid, &int_1);

    // compute 'scale' by the weighted mad
    status = wmad(dat, work, resid, mad_center, mad_NORM_CONSTANT, scale);
    return status;
}

/******************************************************************************\
|* Weighting scheme for iteratively reweighted least squares; weighting with  *|
|* respect to sampling, residual and design space outlyingness                *|
|*                                                                            *|
|* dat       typedef struct regata                                            *|
|* f_wgt_psi function pointer to the wgt/psi-function                         *|
|* resid     residuals, array[n]                                              *|
|* k         tuning constant of the psi-function                              *|
|* type      type of estimator                                                *|
|* robwgt    on return: combined robustness (and sampling) weight             *|
\******************************************************************************/
static inline void weighting_scheme(regdata *dat,
    double (*f_wgt_psi)(double, const double), double* restrict resid,
    double *scale, double *k, int *psi, int *type, double* restrict robwgt)
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
            &lwork, &info_dgels);
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
            work_dgels, &lwork, &info_dgels);

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
            &double_1, resid, &int_1);
        return ROBSURVEY_ERROR_OK;
    }
}

/******************************************************************************\
|* Weighted median of the absolute deviation about zero or about the weighted *|
|* median                                                                     *|
|*                                                                            *|
|* dat      typedef struct regdata                                            *|
|* work     typedef struct workarray                                          *|
|* resid    residuals, array[n]                                               *|
|* median   mad centered about the median (1) or about zero (0)               *|
|* constant normalization constant of the mad                                 *|
|* mad      on return: weighted mad                                           *|
\******************************************************************************/
robsurvey_error_type wmad(regdata *dat, workarray *work, double* restrict resid,
    int* median, double constant, double *mad)
{
    int n = dat->n;
    double* restrict w = dat->w;
    double* restrict work_y = work->work_y;
    double* restrict work_2n = work->work_2n;

    // center the absolute deviations about the weighted median or about zero
    double med, prob = 0.5;
    if (*median) {          // centered about the weighted median
        wquantile_noalloc(resid, w, work_2n, &n, &prob, &med);
        for (int i = 0; i < n; i++)
            work_y[i] = fabs(resid[i] - med);
    } else {                // centered about zero
         for (int i = 0; i < n; i++)
            work_y[i] = fabs(resid[i]);
    }

    // compute mad
    wquantile_noalloc(work_y, w, work_2n, &n, &prob, mad);
    *mad *= constant;

    if (*mad < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;
    else
        return ROBSURVEY_ERROR_OK;
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

/******************************************************************************\
|* cov_rwlslm: covariance matrix of the esimated regression coefficients      *|
|*                                                                            *|
|* resid   residuals, array[n]                                                *|
|* x       design matrix, array[n * p]; on return: the p * p cov              *|
|*         matrix  is stored in x[1..(p * p)]                                 *|
|* xwgt    weights in design space, array[n]                                  *|
|* robwgt  robustness weight, array[n]                                        *|
|* w       sampling weight, array[n]                                          *|
|* k       robustness tuning constant                                         *|
|* scale   estimate of scale                                                  *|
|* scale2  on return: estimate of scale (proposal 2)                          *|
|* n, p    dimensions                                                         *|
|* psi     0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight                *|
|* type    0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.              *|
|* ok      1 = ok; 0 = failure                                                *|
\******************************************************************************/
void cov_rwlslm(double *resid, double *x, double *xwgt, double *robwgt,
    double *w, double *k, double *scale, double *scale2, int *n, int *p,
    int *psi, int *type, int *ok)
{
    double tmp, sum_w = 0.0;
    robsurvey_error_type status;

    *ok = 1;

    // initialize and populate structure with work arrays
    double* restrict work_x = (double*) Calloc(*n * *p, double);
    double* restrict work_y = (double*) Calloc(*n, double);
    workarray wwork;
    workarray *work = &wwork;
    work->work_x = work_x;
//    work->work_y = work_y;

    // determine lwork and allocate work_lapack: dgeqrf (used in inverse_qr)
    int lwork = -1, info;
    F77_CALL(dgeqrf)(n, p, x, n, work_x, work_y, &lwork, &info);
    lwork = (int) work_y[0];
    double* restrict work_lapack = (double*) Calloc(lwork, double);
    work->work_lapack = work_lapack;

    // psi-function: function ptrs
    double (*f_psi)(double, double);
    double (*f_psiprime)(double, double);
    switch (*psi) {
    case 0: // Huber psi
        f_psi = huber_psi;
        f_psiprime = huber_psi_prime;
        break;
    case 1: // asymmetric Huber psi
        f_psi = huber_psi_asym;
        f_psiprime = huber_psi_prime_asym;
        break;
    case 2: // Tukey biweight psi
        f_psi = tukey_psi;
        f_psiprime = tukey_psi_prime;
        break;
    default:
        f_psi = huber_psi;
        f_psiprime = huber_psi_prime;
    }

#if 0
    // type of estimator
    switch (*type) {
    case 0: // M-estimator
        status = cov_m_est();
        break;
    case 1: // Mallows GM-estimator
        break;
    case 2: // Schweppe GM-estimator
        break;
    default:
    }

    if (status != ROBSURVEY_ERROR_OK) {
        *ok = 0;
        PRINT_OUT("Error: %s\n", robsurvey_error(status));
        goto clean_up;
    }

#endif

    // Schweppe GM-est. --------------------------------------------------------
    if (*type == 2) {

        for (int i = 0; i < *n; i++) {
            work_y[i] = resid[i] / *scale;
            sum_w += w[i];
        }

        // compute s_1 and s_2
        double tmp2, z;
        for (int i = 0; i < *n; i++) {
            tmp = 0.0; tmp2 = 0.0;

            if (xwgt[i] > DBL_EPSILON) {
            for (int j = 0; j < *n; j++) {
                z = work_y[j] * xwgt[i];
                tmp += w[j] * f_psiprime(z, *k);
                tmp2 += w[j] * _POWER2(f_psi(z, *k) / xwgt[i]);
            }
            tmp /= sum_w;
            tmp2 /= sum_w;

            } else {
                tmp = 1.0;
                tmp2 = 0.0;
            }

            for (int j = 0; j < *p; j++)        // x := sqrt(s_1 * w) o x
                x[*n * j + i] *= sqrt(tmp * w[i]);

            work_x[i] = tmp2 / tmp;             // temporarily store s_2 / s_1
        }

        Memcpy(work_y, work_x, *n);             // temporarily store s_2 / s_1

        // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
        status = inverse_qr(work, x, n, p, 1);
        if (status != ROBSURVEY_ERROR_OK) {
            *ok = 0;
            PRINT_OUT("Error: %s\n", robsurvey_error(status));
            goto clean_up;
        }

        // pre-multiply Q by sqrt(s2 / s1)
        for (int i = 0; i < *n; i++) {
            tmp = sqrt(work_y[i]);
            for (int j = 0; j < *p; j++)
                x[*n * j + i] *= tmp;           // pre-multiply Q
        }

        // B  := Q * R^{-T} (result -> x)
        double done = 1.0, dzero = 0.0;
        F77_CALL(dtrmm)("R", "U", "T", "N", n, p, &done, work_x, p, x, n);

        // compute B^T * B := (x^T * W * W * x)^{-1}
        *scale2 = _POWER2(*scale) / (1.0 - (double)*p / sum_w);
        F77_CALL(dgemm)("T", "N", p, p, n, scale2, x, n, x, n, &dzero, work_x,
            p);
        *scale2 = *scale;

    // M-est. and Mallows GM-est. ----------------------------------------------
    } else {

        //FIXME: this part of the Mallows estimator

        double Epsi_prime = 0.0, Epsi_prime2 = 0.0;
        for (int i = 0; i < *n; i++) {
            tmp = f_psiprime(resid[i] / *scale, *k);
            Epsi_prime += w[i] * tmp;
            Epsi_prime2 += w[i] * _POWER2(tmp);
            sum_w += w[i];
        }

        Epsi_prime /= sum_w;
        Epsi_prime2 /= sum_w;

        // scale estimate
        *scale2 = 0.0;
        for (int i = 0; i < *n; i++)
            *scale2 += w[i] * _POWER2(robwgt[i] * resid[i]);

        *scale2 /= (sum_w - (double)*p) * _POWER2(Epsi_prime);

        // M-est. -----------------------------------------------------
        if (*type == 0) {
            status = cov_m_est(work, resid, x, robwgt, w, k, scale, scale2, n, p, f_psiprime);
            if (status != ROBSURVEY_ERROR_OK) {
                *ok = 0;
                PRINT_OUT("Error: %s\n", robsurvey_error(status));
                goto clean_up;
            }

        // Mallows GM-est. --------------------------------------------
        } else {

            // this is part of the Mallows estimator

            for (int i = 0; i < *n; i++) {
                tmp = sqrt(w[i] * xwgt[i]);
                for (int j = 0; j < *p; j++)
                    x[*n * j + i] *= tmp;           // pre-multiply x
            }

            // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
            status = inverse_qr(work, x, n, p, 1);
            if (status != ROBSURVEY_ERROR_OK) {
                *ok = 0;
                PRINT_OUT("Error: %s\n", robsurvey_error(status));
                goto clean_up;
            }

            // pre-multiply Q by with sqrt(xwgt)
            for (int i = 0; i < *n; i++) {
                tmp = sqrt(xwgt[i]);
                for (int j = 0; j < *p; j++)
                    x[*n * j + i] *= tmp;           // pre-multiply Q
            }

            // B  := Q * R^{-T} (result -> x)
            double done = 1.0, dzero = 0.0;
            F77_CALL(dtrmm)("R", "U", "T", "N", n, p, &done, work_x, p, x, n);

            // compute B^T * B := (x^T * W * W * x)^{-1}
            F77_CALL(dgemm)("T", "N", p, p, n, scale2, x, n, x, n, &dzero,
            work_x, p);
        }
    }
    Memcpy(x, work_x, *p * *p);                     // store in x[1..(p * p)]

clean_up:
    Free(work_lapack); Free(work_x); Free(work_y);
}

/******************************************************************************\
|* Asymptotic covariance matrix of the M-estimator                            *|
|*                                                                            *|
|* data       typedef struct regdata                                          *|
|* work       typedef struct workarray                                        *|
|* resid      residuals, array[n]                                             *|
|* k          robustness tuning constant                                      *|
|* scale      weighted mad                                                    *|
|* scale2     on return: estimate of regression scale                         *|
|* f_psiprime function ptr to the psi-prime function                          *|
\******************************************************************************/
//FIXME: add regata structure
robsurvey_error_type cov_m_est(workarray *work, double *resid,
    double *x, double *robwgt,
    double *w, double *k, double *scale, double *scale2, int *n, int *p,
    double (*f_psiprime)(double, const double))
{
    // E(psi') and E(psi')^2
    double tmp, Epsi_prime = 0.0, Epsi_prime2 = 0.0, sum_w = 0.0;
    for (int i = 0; i < *n; i++) {
        tmp = (*f_psiprime)(resid[i] / *scale, *k);
        Epsi_prime += w[i] * tmp;
        Epsi_prime2 += w[i] * _POWER2(tmp);
        sum_w += w[i];
    }
    Epsi_prime /= sum_w;
    Epsi_prime2 /= sum_w;

    // scale estimate
    *scale2 = 0.0;
    for (int i = 0; i < *n; i++)
        *scale2 += w[i] * _POWER2(robwgt[i] * resid[i]);

    *scale2 /= (sum_w - (double)*p) * _POWER2(Epsi_prime);

    // correction factor (see Huber, 1981, p. 172-174)
    double kappa = 1.0 + (double)*p / sum_w * (Epsi_prime2 /
        _POWER2(Epsi_prime) - 1.0) * (double)*n / (double)(*n - 1);
    *scale2 *= _POWER2(kappa);

    if (*scale2 < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;

    // pre-multiply x by sqrt(weight)
    for (int i = 0; i < *n; i++) {
        tmp = sqrt(w[i]);
        for (int j = 0; j < *p; j++)
            x[*n * j + i] *= tmp;
    }

    // inverse of x^T * x (using QR factorization of)
    robsurvey_error_type status = inverse_qr(work, x, n, p, 0);
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    double* restrict work_x = work->work_x;
    F77_CALL(dtrmm)("R", "U", "T", "N", p, p, scale2, work_x, p, work_x, p);

    return ROBSURVEY_ERROR_OK;
}

/******************************************************************************\
|* Inverse of R matrix and Q matrix of the QR factorization                   *|
|*                                                                            *|
|* x        on return; Q matrix, array[n * p]                                 *|
|* work     typedef struct workarray                                          *|
|* n, p     dimensions                                                        *|
|* qmatrix  toggle whether Q matrix is computed : 0 = no; 1 = yes             *|
|*                                                                            *|
|* NOTE: array x will be overwritten                                          *|
\******************************************************************************/
robsurvey_error_type inverse_qr(workarray *work, double *x, int *n, int *p,
    int qmatrix)
{
    int lwork = work->lwork;
    int info = 1;
    int offset = _POWER2(*p);
    double* restrict R = work->work_x;                  // R matrix of QR
    double* restrict work_dgeqrf = work->work_lapack;

    // QR factorization
    F77_CALL(dgeqrf)(n, p, x, n, R + offset, work_dgeqrf, &lwork, &info);
    if (info != 0)
        return ROBSURVEY_ERROR_QR_DGEQRF;

    for (int i = 0; i < *p * *p; i++)                   // prepare matrix R
        R[i] = 0.0;

    for (int i = 0; i < *p; i++)                        // extract matrix R
        for (int j = 0; j < i + 1; j++)
            R[j + i * *p] = x[j + i * *n];

    F77_CALL(dtrtri)("U", "N", p, R, p, &info);         // inverse of R
    if (info != 0)
        return ROBSURVEY_ERROR_QR_DTRTRI;

    if (qmatrix) {
        F77_CALL(dorgqr)(n, p, p, x, n, R + offset,     // extract matrix Q
            work_dgeqrf, &lwork, &info);
        if (info != 0)
            return ROBSURVEY_ERROR_QR_DORGQR;
    }
    return ROBSURVEY_ERROR_OK;
}
#undef _POWER2
#undef PRINT_OUT

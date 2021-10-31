/* Functions to compute the covariance matrix of weighted (generalized)
   regression M-estimators

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

#include "regression_cov.h"

// some macros
#define _POWER2(_x) ((_x) * (_x))

// variadic arguments in macros are supported by gcc (>=3.0), clang (all
// versions), visual studio (>=2005); the version with ##__VA_ARGS__ (which
// will silently eliminate the trailing comma) is not portable (clang complains)
#define PRINT_OUT(...) Rprintf(__VA_ARGS__)

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

// declaration
robsurvey_error_type variance_est(regdata*, double*restrict, double* restrict,
    double*, double*, double*, double (*)(double, const double));
robsurvey_error_type cov_m_est(regdata*, workarray*, double* restrict,
    double* restrict, double*, double*, double*, double (*)(double,
    const double));
robsurvey_error_type cov_mallows_gm_est(regdata*, workarray*, double*,
    double*, double*, double*, double*, double (*)(double, const double));
robsurvey_error_type cov_schweppe_gm_est(regdata*, workarray*,
    double* restrict, double* restrict, double*, double*, double*,
    double (*)(double, const double), double (*)(double,
    const double));
robsurvey_error_type inverse_qr(workarray*, double* restrict, int*, int*, int);

/******************************************************************************\
|* model-based covariance matrix of the esimated regression coefficients      *|
|*                                                                            *|
|* resid   residuals, array[n]                                                *|
|* x       design matrix, array[n * p]; on return: the p * p cov              *|
|*         matrix  is stored in x[1..(p * p)]                                 *|
|* xwgt    weights in design space, array[n]                                  *|
|* robwgt  robustness weight, array[n]                                        *|
|* w       sampling weight, array[n]                                          *|
|* k       robustness tuning constant                                         *|
|* scale   estimate of scale                                                  *|
|* scale2  on return: estimate of regression variance (proposal 2)            *|
|* n, p    dimensions                                                         *|
|* psi     0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight                *|
|* type    0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.              *|
|* ok      on return: 1 = ok; 0 = failure                                     *|
\******************************************************************************/
void cov_reg_model(double *resid, double *x, double *xwgt, double *robwgt,
    double *w, double *k, double *scale, double *scale2, int *n, int *p,
    int *psi, int *type, int *ok)
{
    // initialize and populate structure with regression-specific data
    regdata data;
    regdata *dat = &data;
    dat->n = *n;
    dat->p = *p;
    dat->x = x;
    dat->w = w;
    dat->xwgt = xwgt;

    // initialize and populate structure with work arrays
    double* restrict work_x = (double*) Calloc(*n * *p, double);
    double* restrict work_y = (double*) Calloc(*n, double);
    workarray wwork;
    workarray *work = &wwork;
    work->work_x = work_x;
    work->work_y = work_y;

    // determine lwork and allocate work_lapack: dgeqrf (used in inverse_qr)
    int lwork = -1, info;
    F77_CALL(dgeqrf)(n, p, x, n, work_x, work_y, &lwork, &info);
    lwork = (int) work_y[0];
    work->lwork = lwork;
    double* restrict work_lapack = (double*) Calloc(lwork, double);
    work->work_lapack = work_lapack;

    // psi and psi-prime functions
    double (*f_psi)(double, double);
    f_psi = get_psi_function(*psi);

    double (*f_psiprime)(double, double);
    f_psiprime = get_psi_prime_function(*psi);

    // type of estimator
    robsurvey_error_type status;

    switch (*type) {
    case 0: // M-estimator
        status = cov_m_est(dat, work, resid, robwgt, k, scale, scale2,
            f_psiprime);
        break;
    case 1: // Mallows GM-estimator
        status = cov_mallows_gm_est(dat, work, resid, robwgt, k, scale,
            scale2, f_psiprime);
        break;
    case 2: // Schweppe GM-estimator
        status = cov_schweppe_gm_est(dat, work, resid, robwgt, k, scale,
            scale2, f_psiprime, f_psi);
        break;
    default: // M-estimator
        status = cov_m_est(dat, work, resid, robwgt, k, scale, scale2,
            f_psiprime);
    }

    if (status != ROBSURVEY_ERROR_OK) {
        *ok = 0;
        PRINT_OUT("Error: %s\n", robsurvey_error(status));
        goto clean_up;
    }

    *ok = 1;

    // copy cov matrix to x[1..(p * p)]
    Memcpy(x, work_x, *p * *p);

clean_up:
    Free(work_lapack); Free(work_x); Free(work_y);
}

/******************************************************************************\
|* Asymptotic covariance matrix of the M-estimator                            *|
|*                                                                            *|
|* data       typedef struct regdata                                          *|
|* work       typedef struct workarray                                        *|
|* resid      residuals, array[n]                                             *|
|* robwgt     robustness weight, array[n]                                     *|
|* k          robustness tuning constant                                      *|
|* scale      weighted mad                                                    *|
|* scale2     on return: estimate of regression variance                      *|
|* f_psiprime function ptr to the psi-prime function                          *|
\******************************************************************************/
robsurvey_error_type cov_m_est(regdata *dat, workarray *work,
    double* restrict resid, double* restrict robwgt, double *k,
    double *scale, double *scale2, double (*f_psiprime)(double, const double))
{
    int n = dat->n, p = dat->p;
    double* restrict x = dat->x;
    double* restrict w = dat->w;
    robsurvey_error_type status;

    // estimate of variance
    status = variance_est(dat, resid, robwgt, scale, scale2, k, f_psiprime);
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    // pre-multiply x by sqrt(weight)
    double tmp;
    for (int i = 0; i < n; i++) {
        tmp = sqrt(w[i]);
        for (int j = 0; j < p; j++)
            x[n * j + i] *= tmp;
    }

    // inverse of x^T * x (using QR factorization)
    status = inverse_qr(work, x, &n, &p, 0);    // inverse of R matrix
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    double* restrict work_x = work->work_x;     // scale2 * R^{-1} * R^{-T}
    F77_CALL(dtrmm)("R", "U", "T", "N", &p, &p, scale2, work_x, &p, work_x,
        &p FCONE FCONE FCONE FCONE);

    return ROBSURVEY_ERROR_OK;
}

/******************************************************************************\
|* Estimate of regression variance                                            *|
|*                                                                            *|
|* dat        typedef struct regdata                                          *|
|* resid      residuals, array[n]                                             *|
|* robwgt     robustness weight, array[n]                                     *|
|* scale      weighted mad                                                    *|
|* scale2     on return: estimate of regression variance                      *|
|* k          robustness tuning constant                                      *|
|* f_psiprime function ptr to the psi-prime function                          *|
\******************************************************************************/
robsurvey_error_type variance_est(regdata *dat, double* restrict resid,
    double* restrict robwgt, double *scale, double *scale2, double *k,
    double (*f_psiprime)(double, const double))
{
    int n = dat->n, p = dat->p;
    double dbl_p = (double)p;
    double* restrict w = dat->w;

    // E(psi') and E(psi')^2
    double tmp, Epsi_prime = 0.0, Epsi_prime2 = 0.0, sum_w = 0.0;
    for (int i = 0; i < n; i++) {
        tmp = (*f_psiprime)(resid[i] / *scale, *k);
        Epsi_prime += w[i] * tmp;
        Epsi_prime2 += w[i] * _POWER2(tmp);
        sum_w += w[i];
    }
    Epsi_prime /= sum_w;
    Epsi_prime2 /= sum_w;

    // variance estimate
    *scale2 = 0.0;
    for (int i = 0; i < n; i++)
        *scale2 += w[i] * _POWER2(robwgt[i] * resid[i]);
    *scale2 /= (sum_w - dbl_p) * _POWER2(Epsi_prime);

    // correction factor (see Huber, 1981, p. 172-174); with modification:
    // numerator is (N-1) not N; viz. MASS::rlm
    double kappa = 1.0 + dbl_p / (sum_w - 1.0) * (Epsi_prime2 /
        _POWER2(Epsi_prime) - 1.0);
    *scale2 *= _POWER2(kappa);

    if (*scale2 < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;
    else
        return ROBSURVEY_ERROR_OK;
}

/******************************************************************************\
|* Asymptotic covariance matrix of the Schweppe GM-estimator                  *|
|*                                                                            *|
|* data       typedef struct regdata                                          *|
|* work       typedef struct workarray                                        *|
|* resid      residuals, array[n]                                             *|
|* robwgt     robustness weight, array[n]                                     *|
|* k          robustness tuning constant                                      *|
|* scale      weighted mad                                                    *|
|* scale2     on return: estimate of regression variance                      *|
|* f_psiprime function ptr to the psi-prime function                          *|
|* f_psi      function ptr to the psi-function                                *|
\******************************************************************************/
robsurvey_error_type cov_schweppe_gm_est(regdata *dat, workarray *work,
    double* restrict resid, double* restrict robwgt, double *k, double *scale,
    double *scale2, double (*f_psiprime)(double, const double),
    double (*f_psi)(double, const double))
{
    int n = dat->n, p = dat->p;
    double* restrict x = dat->x;
    double* restrict w = dat->w;
    double* restrict xwgt = dat->xwgt;
    double* restrict work_x = work->work_x;
    double* restrict work_y = work->work_y;

    double sum_w = 0.0;
    for (int i = 0; i < n; i++) {
        work_y[i] = resid[i] / *scale;
        sum_w += w[i];
    }

    // compute s_1 and s_2
    double tmp, tmp2, z;
    for (int i = 0; i < n; i++) {
        tmp = 0.0; tmp2 = 0.0;

        if (xwgt[i] > DBL_EPSILON) {
            for (int j = 0; j < n; j++) {
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

        // x := sqrt(s_1 * w) o x
        for (int j = 0; j < p; j++)
            x[n * j + i] *= sqrt(tmp * w[i]);

        // temporarily store s_2 / s_1
        work_x[i] = tmp2 / tmp;
    }

    // temporarily store s_2 / s_1
    Memcpy(work_y, work_x, n);

    // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
    robsurvey_error_type status = inverse_qr(work, x, &n, &p, 1);
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    // pre-multiply Q by sqrt(s2 / s1)
    for (int i = 0; i < n; i++) {
        tmp = sqrt(work_y[i]);
        for (int j = 0; j < p; j++)
            x[n * j + i] *= tmp;           // pre-multiply Q
    }

    // B  := Q * R^{-T} (result -> x)
    double done = 1.0, dzero = 0.0;
    F77_CALL(dtrmm)("R", "U", "T", "N", &n, &p, &done, work_x, &p, x,
        &n FCONE FCONE FCONE FCONE);

    // compute B^T * B := (x^T * W * W * x)^{-1}
    *scale2 = _POWER2(*scale) / (1.0 - (double)p / sum_w);
    if (*scale2 < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;

    F77_CALL(dgemm)("T", "N", &p, &p, &n, scale2, x, &n, x, &n, &dzero,
        work_x, &p FCONE FCONE);
    *scale2 = *scale;

    return ROBSURVEY_ERROR_OK;
}

/******************************************************************************\
|* Asymptotic covariance matrix of the Mallows GM-estimator                   *|
|*                                                                            *|
|* data       typedef struct regdata                                          *|
|* work       typedef struct workarray                                        *|
|* resid      residuals, array[n]                                             *|
|* robwgt     robustness weight, array[n]                                     *|
|* k          robustness tuning constant                                      *|
|* scale      weighted mad                                                    *|
|* scale2     on return: estimate of regression variance                      *|
|* f_psiprime function ptr to the psi-prime function                          *|
\******************************************************************************/
robsurvey_error_type cov_mallows_gm_est(regdata *dat, workarray *work,
    double* restrict resid, double* restrict robwgt, double *k,
    double *scale, double *scale2, double (*f_psiprime)(double, const double))
{
    int n = dat->n, p = dat->p;
    double* restrict x = dat->x;
    double* restrict w = dat->w;
    double* restrict xwgt = dat->xwgt;
    double* restrict work_x = work->work_x;
    robsurvey_error_type status;

    // estimate of variance
    status = variance_est(dat, resid, robwgt, scale, scale2, k, f_psiprime);
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    // pre-multiply x by sqrt(xwgt)
    double tmp;
    for (int i = 0; i < n; i++) {
        tmp = sqrt(w[i] * xwgt[i]);
        for (int j = 0; j < p; j++)
            x[n * j + i] *= tmp;
    }

    // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
    status = inverse_qr(work, x, &n, &p, 1);
    if (status != ROBSURVEY_ERROR_OK)
        return status;

    // pre-multiply Q by with sqrt(xwgt)
    for (int i = 0; i < n; i++) {
        tmp = sqrt(xwgt[i]);
        for (int j = 0; j < p; j++)
            x[n * j + i] *= tmp;
    }

    // B := Q * R^{-T} (result -> x)
    double done = 1.0, dzero = 0.0;
    F77_CALL(dtrmm)("R", "U", "T", "N", &n, &p, &done, work_x, &p, x,
        &n FCONE FCONE FCONE FCONE);

    // compute B^T * B := (x^T * W * W * x)^{-1}
    F77_CALL(dgemm)("T", "N", &p, &p, &n, scale2, x, &n, x, &n, &dzero,
        work_x, &p FCONE FCONE);

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
robsurvey_error_type inverse_qr(workarray *work, double* restrict x, int *n,
    int *p, int qmatrix)
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

    F77_CALL(dtrtri)("U", "N", p, R, p, &info FCONE FCONE); // inverse of R
    if (info != 0)
        return ROBSURVEY_ERROR_QR_DTRTRI;

    if (qmatrix) {                                      // extract matrix Q
        F77_CALL(dorgqr)(n, p, p, x, n, R + offset,
            work_dgeqrf, &lwork, &info);
        if (info != 0)
            return ROBSURVEY_ERROR_QR_DORGQR;
    }
    return ROBSURVEY_ERROR_OK;
}

/******************************************************************************\
|* design-based estimate of the regression covariance matrix                  *|
|*                                                                            *|
|* x      model design matrix, array[n, p]                                    *|
|* w      weights, array[n]                                                   *|
|* xwgt   weight in the model's design space                                  *|
|* resid  residual, array[n]                                                  *|
|* scale  estimate of regressions scale                                       *|
|* k      robustness tuning constant                                          *|
|* psi    type of psi-function                                                *|
|* type   0: M-estimator, 1: Mallows GM-, and 2: Schweppe GM-estimator        *|
|* n, p   dimensions                                                          *|
|* ok     on return: 1 = computation is ok; 0 = failure                       *|
|* mat    on return: covariance matrix; on entry: covariance matrix of the    *|
|*        estimated total                                                     *|
|*                                                                            *|
|* NOTE: matrix x will be overwritten                                         *|
\******************************************************************************/
void cov_reg_design(double *x, double *w, double *xwgt, double *resid,
    double *scale, double *k, int *psi, int *type, int *n, int *p, int *ok,
    double *mat)
{
    int info = 1;
    double d_one = 1.0, d_zero = 0.0;
    *ok = 1;
    // allocate memory
    double* M = (double*) Calloc(*p * *p, double);
    double* work_pp = (double*) Calloc(*p * *p, double);
    double* work_np = (double*) Calloc(*n * *p, double);

    // GM-estimators
    if (*type == 1) {                   // Mallows GM-estimator
        for (int i = 0; i < *n; i++)
            w[i] *= xwgt[i];
    }
    if (*type == 2) {                   // Schweppe GM-estimator
        for (int i = 0; i < *n; i++) {
            if (fabs(xwgt[i]) < DBL_EPSILON)
                resid[i] = 0.0;
            else
                resid[i] /= xwgt[i];
        }
    }

    // psi-prime function
    double (*f_psiprime)(double, double);
    f_psiprime = get_psi_prime_function(*psi);

    //--------------------------------------------
    // compute matrix M (i.e., sum w[i] * psi_prime[i] * x[i] * x[i]^T)
    double tmp;
    for (int i = 0; i < *n; i++) {
        tmp = w[i] * (*f_psiprime)(resid[i] / *scale, *k);
        for (int j = 0; j < *p; j++)
            work_np[i + *n * j] = x[i + *n * j] * tmp;
    }
    F77_CALL(dgemm)("T", "N", p, p, n, &d_one, work_np, n, x, n, &d_zero,
        M, p FCONE FCONE);

    // overwrite M with Cholesky factor of M such that M = U * U^T, where
    // U is upper triangular
    F77_CALL(dpotrf)("U", p, M, p, &info FCONE);
    if (info != 0) {
        PRINT_OUT("Error in dpotrf (M matrix)\n");
        *ok = 0;
        goto clean_up;
    }

    // overwrite M with upper triangular matrix of the inverse of M (using
    // Cholesky factor of M)
    F77_CALL(dpotri)("U", p, M, p, &info FCONE);
    if (info != 0) {
        PRINT_OUT("Error in dpotri (M matrix)\n");
        *ok = 0;
        goto clean_up;
    }

    // compute B := M^{-1} * mat
    F77_CALL(dsymm)("L", "U", p, p, &d_one, M, p, mat, p, &d_zero, work_pp,
        p FCONE FCONE);

    // compute B * M^{-1} := covariance
    F77_CALL(dsymm)("R", "U", p, p, &d_one, M, p, work_pp, p, &d_zero, mat,
        p FCONE FCONE);

    //--------------------------------------------
    // covariance matrix

clean_up:
    Free(work_pp); Free(work_np); Free(M);
}
#undef _POWER2
#undef PRINT_OUT

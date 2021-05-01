/* Functions to compute weighted (generalized) regression M-estimators,
   winsorized, and trimmed estimators of location

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

#include "fitwls.h"

/******************************************************************************\
|* Weighted least squares estimate                                            *|
|*                                                                            *|
|*  x       vectorized design matrix, array[n * p]                            *|
|*  work_x  array[n * p], on return: QR factor (dgels)                        *|
|*  y       response vector, array[n]                                         *|
|*  work_y  array[n]                                                          *|
|*  w       weights vector, array[n]                                          *|
|*  resid   on return: residuals vector, array[n]                             *|
|*  beta0   on return: coefficient vector, array[p]                           *|
|*  n, p    array dimensions                                                  *|
|*  work    work array, array[lwork] used for QR factorization                *|
|*  lwork   size of array 'work' (if < 0, then 'dgels' determines the optimal *|
|*          size)                                                             *|
|*  info    on return: info on fitwls (error if != 0)                         *|
\******************************************************************************/
void fitwls(double* restrict x, double* restrict work_x, double* restrict y,
    double* restrict work_y, double* restrict w, double* restrict resid,
    double* restrict beta0, int *n, int *p, double *work, int *lwork,
    int *info)
{
    // define constants for the call of 'dgels'
    const int int_1 = 1;
    int info_dgels = 1;
    *info = 0;

    // STEP 0: determine the optimal size of array 'work'
    if (*lwork < 0) {
        F77_CALL(dgels)("N", n, p, &int_1, x, n, y, n, work, lwork,
            &info_dgels);
        *lwork = (int) work[0];

	// STEP 1: compute least squares fit
	} else {
        // pre-multiply the design matrix and the response vector by sqrt(w)
        double tmp;
        for (int i = 0; i < *n; i++) {
            tmp = sqrt(w[i]);
            work_y[i] = y[i] * tmp;

            for (int j = 0; j < *p; j++)
                work_x[*n * j + i] = x[*n * j + i] * tmp;
        }

        // compute the (weighted) least squares estimate (LAPACK::dgels),
        // solves minimize |B - A*X| for X (using QR factorization)
        F77_CALL(dgels)("N", n, p, &int_1, work_x, n, work_y, n, work, lwork,
            &info_dgels);

        // dgels is not well suited as a rank-revealing procedure; i.e., INFO<0
        // iff a diagonal element of the R matrix is exactly 0. This is not
        // helpful; hence, we check the diagonal elements of R separately and
        // issue and error flag if any(abs(diag(R))) is close to zero
        for (int i = 0; i < *p; i++) {
            if (fabs(work_x[(*n + 1) * i]) < sqrt(DBL_EPSILON)) {
                *info = 1;
                return;
            }
        }

        Memcpy(beta0, work_y, *p);  // retrieve 'betacoefficients'

        // compute the residuals (BLAS::dgemv): y = alpha*A*x + beta*y
        const double double_minus1 = -1.0, double_1 = 1.0;
        Memcpy(resid, y, *n);
        F77_CALL(dgemv)("N", n, p, &double_minus1, x, n, beta0, &int_1,
            &double_1, resid,&int_1);
    }
}

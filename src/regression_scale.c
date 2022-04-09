/* Functions to compute weighted robust scales estimators

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

#include "regression_scale.h"

// variadic arguments in macros are supported by gcc (>=3.0), clang (all
// versions), visual studio (>=2005); the version with ##__VA_ARGS__ (which
// will silently eliminate the trailing comma) is not portable (clang complains)
#define PRINT_OUT(...) Rprintf(__VA_ARGS__)

/******************************************************************************\
|* Weighted normalized IQR                                                    *|
|*                                                                            *|
|* dat          typedef struct regdata                                        *|
|* work         typedef struct workarray                                      *|
|* resid        residuals, array[n]                                           *|
|* iqr          on return: iqr scale estimate                                 *|
\******************************************************************************/
robsurvey_error_type wiqr(regdata *dat, workarray *work, double* restrict resid,
    double *iqr)
{
    int n = dat->n;
    double p25 = 0.25, x25 = 0.0;
    double p75 = 0.75, x75 = 0.0;
    double* restrict w = dat->w;
    double* restrict work_2n = work->work_2n;

    wquantile_noalloc(resid, w, work_2n, &n, &p25, &x25);
    wquantile_noalloc(resid, w, work_2n, &n, &p75, &x75);
    *iqr = (x75 - x25) * iqr_NORM_CONSTANT;

    if (*iqr < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;
    else
        return ROBSURVEY_ERROR_OK;
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
    int *median, double constant, double *mad)
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

    // compute normalized mad
    wquantile_noalloc(work_y, w, work_2n, &n, &prob, mad);
    *mad *= constant;

    if (*mad < DBL_EPSILON)
        return ROBSURVEY_ERROR_SCALE_ZERO;
    else
        return ROBSURVEY_ERROR_OK;
}

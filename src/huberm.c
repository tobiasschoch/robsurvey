/* Weighted Huber proposal 2 estimator of location and scale

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

#include <R.h>
#include "wquantile.h"

// some macros
#define _WGT_HUBER(_x, _k) ((fabs(_x) >= _k) ? _k / fabs(_x) : 1.0)
#define _POWER2(_x) ((_x) * (_x))
#define _MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#define _MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))

// declaration of 'local' functions (inline imperative is GCC specific)
static inline double kappa_huber(const double) __attribute__((always_inline));

/*****************************************************************************\
|* huber mean and scale estimator: proposal 2 estimator (using weights;	     *| 
|* the same tuning constant 'k' used for location and scale; the scale	     *|
|* is normalized to be a Fisher consistent estimator at the Gaussian	     *| 
|* core model)								     *|
|*									     *|
|*    x	       array							     *|
|*    w	       array of weights				 		     *|
|*    robwgt   array of robustness weights (on return)	  		     *|
|*    k	       robustness tuning constant (location and scale)		     *|
|*    loc      robust location (on return)				     *|
|*    scale    robust scale (on return)					     *|
|*    n	       array dimension						     *|
|*    maxit    maximum number of iteration (termination rule); on return:    *|
|*	       effective number of iterations				     *|
|*    tol      numerical convergence tolerance for iterative refinement	     *|
|*	       iterations						     *|
\*****************************************************************************/
void huberm(double *x, double *w, double *robwgt, double *k, double *loc, 
   double *scale, int *n, int *maxit, const double *tol)
{
   int iter;
   double loc0, scale0, tmp, tmp_loc, tmp_scale, kappa, wtotal; 
   double *x_wins, *work_2n;

   if (*n == 1) {
      *loc = *x;
      *scale = 0.0;
   } else {
      // initialize location (weighted median)
      double p50 = 0.5;

      // quantile   
      work_2n = (double*) Calloc(2 * *n, double);
      wquantile_noalloc(x, w, work_2n, n, &p50, &loc0);

      // initialize variable for winsorized x-variable 
      x_wins = (double*) Calloc(*n, double);

      // initialize scale estimate by (weighted) interquartile range (IQR is
      // normalized to be a Fisher consistent estimate of the scale at the 
      // Gaussian distr.)
      double p25 = 0.25, x25 = 0.0;
      double p75 = 0.75, x75 = 0.0; 
      wquantile_noalloc(x, w, work_2n, n, &p25, &x25); 
      wquantile_noalloc(x, w, work_2n, n, &p75, &x75); 
      scale0 = (x75 - x25) / 1.349;
 
      // stop if IQR is zero 
      if (scale0 < DBL_EPSILON) 
	 *loc = x[1];
      else {
	 // compute the weight total
	 wtotal = 0.0;
	 for (int i = 0; i < *n; i++) 
	    wtotal += w[i];

	 // compute consistency correction term
	 kappa = kappa_huber(*k);

	 // loop   
	 for (iter = 0; iter < *maxit; iter ++) {

	    // update location and scale 
	    tmp_loc = 0.0; tmp_scale = 0.0;
	    tmp = *k * scale0;

	    for (int i = 0; i < *n; i++) {
	       // winsorized x-variable
	       x_wins[i] = _MIN(_MAX(loc0 - tmp, x[i]), loc0 + tmp);
	       // location step
	       tmp_loc += w[i] * x_wins[i];
	       // scale
	       tmp_scale += w[i] * _POWER2(x_wins[i] - loc0);
	    } 

	    *loc = tmp_loc / wtotal;
	    *scale = tmp_scale / wtotal;
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

	 Free(x_wins);
      }

      Free(work_2n);
   }
}

/*****************************************************************************\
|*  This function computes the multiplicative consistency correction term    *|
|*  for the Huber 'Proposal 2' type of scale estimator.                      *|
\*****************************************************************************/
static inline double kappa_huber(const double k) 
{
   double pdf_k, cdf_k; 
   if (k < 10.0) {
      pdf_k = exp(-_POWER2(k) / 2.0) / 2.5066282746;

      // lower tail of cdf; i.e., 1 - Phi(.)
      cdf_k = 0.5 - 0.5 * erf(k / 1.4142135623);
      return 1.0 - 2.0 * (k * pdf_k + (1.0 - k * k) * cdf_k);
   } else 
      return 1.0;
} 

#undef _WGT_HUBER
#undef _POWER2
#undef _MIN
#undef _MAX 

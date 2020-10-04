/*****************************************************************************\
|* PROJECT  robsurvey							     *|
|* SUBEJCT  basic functions						     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include "robsurvey.h"

/* some macros */
#define _WGT_HUBER(_x, _k) ((fabs(_x) >= _k) ? _k / fabs(_x) : 1.0)
#define _POWER2(_x) ((_x) * (_x))
#define _MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#define _MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))

/*declaration of 'local' functions (inline imperative is GCC specific)*/
static inline void fitwls(double*, double*, double*, double*, double*, double*, 
   double*, int*, int*, double*, int*, int*) __attribute__((always_inline));
static inline double kappa_huber(const double) __attribute__((always_inline));
static inline double euclidean_norm(const double*, const double*, int)
    __attribute__((always_inline));
static inline double wmad(double*, double*, int, double)
   __attribute__((always_inline));
static inline void robweight(double*, double*, double*, double*, double*, int*, 
   int*, int*) __attribute__((always_inline));

/*****************************************************************************\
|*  rwlslm: iteratively reweighted least squares			     *|
|*                                                                           *|
|*    x		  vectorized design matrix, array[n * p]		     *|
|*    y		  response vector, array[n]	   			     *| 
|*    w		  weights, array[n]					     *|
|*    resid	  on return: residuals vector, array[n]			     *|
|*    robwgt	  on return: robustness weights, array[n]		     *|
|*    xwgt	  weight in design space, GM-estimator, array[n]	     *|
|*    n		  (array dimensions)					     *|
|*    p		  (array dimensions)					     *|
|*    k		  robustness tuning constant				     *|
|*    beta0	  on return: coefficient vector, array[p]		     *|
|*    scale	  on return: scale estimate				     *|
|*    maxit	  max iterations (on input); iterations (on return)	     *|
|*    tol	  numerical tolerance criterion (stoping rule in rwls	     *|
|*		  updating rule)					     *|
|*    psi	  0 = Huber, 1 = asymmetric Huber, 2 = Tukey	     	     *|
|*    type	  0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.	     *|
|*    psi2	  expected value of psi^2, i.e.: sum(w * psi^2)		     *|
|*    psiprime expected value of psi', i.e..: sum(w * psi')		     *|
\*****************************************************************************/
void rwlslm(double *x, double *y, double *w, double *resid, double *robwgt, 
   double *xwgt, int *n, int *p, double *k, double *beta0, double *scale, 
   int *maxit, double *tol, int *psi, int *type, double *psi2, double *psiprime)
{
   int info = 0, iterations = 0, converged = 0;
   double mad_const = 1.4826;	 
   double *wx, *wy, *work, *beta_new;


   // STEP 0 preparations 
   beta_new = (double*) Calloc(*p, double);
   wx = (double*) Calloc(*n * *p, double);	// used in fitwls
   wy = (double*) Calloc(*n, double);		// used in fitwls

   // determine optimal size of array 'work' and allocat it 
   int lwork = -1; 
   fitwls(x, wx, y, wy, w, resid, beta0, n, p, wx, &lwork, &info);
   work = (double*) Calloc(lwork, double);

   // preparations for GM-estimators
   if (*type == 1) {	   // Mallows: xwgt combined with sampling weight
      for (int i = 0; i < *n; i++)
	 w[i] *= xwgt[i];
      
      // Mallows type has special normalization constant (see 'mallows.c')
      mad_const = mallows_mad_normalization(xwgt, n);
   }
   if (*type == 2) {	   // Schweppe: xwgt turned into multiplicative weight 
      for (int i = 0; i < *n; i++) {
   	 if (fabs(xwgt[i]) < DBL_EPSILON)
	    xwgt[i] = 0.0;	
	 else 
	    xwgt[i] = 1.0 / xwgt[i];
      }
   }

   // STEP 1: initialize beta by weighted least squares
   fitwls(x, wx, y, wy, w, resid, beta0, n, p, work, &lwork, &info);
   if (info > 0){
      error("Error: the design matrix is rank deficient (or nearly so)\n");
      *maxit = 0; 
      return;
   }

   // STEP 2: initialize scale by weighted MAD (ignore that Mallows is special)
   *scale = wmad(resid, w, *n, mad_const);
   if (*scale < DBL_EPSILON) {
      error("Error: the estimate of scale is zero (or nearly so)\n");
      *maxit = 0;
      return;
   }

   // STEP 3: irls updating
   while (!converged && ++iterations < *maxit) {
      // compute irwls weights 
      Memcpy(robwgt, w, *n);
      robweight(resid, robwgt, xwgt, k, scale, n, psi, type);

      // update beta and residuals
      fitwls(x, wx, y, wy, robwgt, resid, beta_new, n, p, work, &lwork, &info);

      // update scale
      if (*type == 1) {				      // Mallows
	 for (int i = 0; i < *n; i++)
	    resid[i] *= sqrt(xwgt[i]);
	 *scale = wmad(resid, w, *n, mad_const); 
      } else					      // otherwise
	 *scale = wmad(resid, w, *n, mad_const); 

      // check for convergence
      converged = (euclidean_norm(beta0, beta_new, *p) < *tol * *scale) ? 1: 0;

      // prepare the next while run 
      Memcpy(beta0, beta_new, *p); 
   }
   *maxit = (converged) ? iterations : 0;

   // compute robustness weights (using MAD estimate of scale and beta) 
   for (int i = 0; i < *n; i++)
      robwgt[i] = 1.0;
   robweight(resid, robwgt, xwgt, k, scale, n, psi, type);

//FIXME

   // compute 'Proposal 2' estimate of scale 
   double tmp = 0.0;
   for (int i = 0; i < *n; i++)  
      tmp += w[i] * _POWER2(robwgt[i] * resid[i]); 

   *scale = sqrt(tmp / (*n * kappa_huber(*k)));

   // empirical estimate of E(psi^2) 
   *psi2 = 0.0;

   // empirical estimate of E(psi')
   *psiprime = 0.0;

   Free(beta_new); Free(wx); Free(wy); Free(work);
}

/*****************************************************************************\
|*  Weighted least squares estimate			                     *|
|*                                                                           *|
|*    x           vectorized design matrix, array[n * p]                     *|
|*    wx          empty array[n * p]					     *|
|*    y           response vector, array[n]                                  *|
|*    wy          empty array[n]					     *|
|*    w           weights vector, array[n]                                   *|
|*    resid       on return: residuals vector, array[n]                      *|
|*    beta0       on return: coefficient vector, array[p]                    *|
|*    n, p        (array dimensions)                                         *|
|*    work        work array, array[lwork]                                   *|
|*    lwork	  size of array 'work' (if < 0, then 'dgels' determines      *|
|*                optimal size)                                              *|
|*    info	  on return: info on fitwls (error if != 0)		     *|
\*****************************************************************************/
static inline void fitwls(double *x, double *wx, double *y, double *wy, 
   double *w, double *resid, double *beta0, int *n, int *p, double *work, 
   int *lwork, int *info)
{
   // define constants for the call of 'dgels'
   const int int_1 = 1;
   int info_dgels = 1;
   *info = 0;

   // STEP 0: determine the optimal size of array 'work'
   if (*lwork < 0) {
      F77_CALL(dgels)("N", n, p, &int_1, x, n, y, n, work, lwork, &info_dgels);	 
      *lwork = (int) work[0]; 

   // STEP 1: compute least squares fit
   } else {
      // pre-multiply the design matrix and the response vector by sqrt(w) 
      double tmp;
      for (int i = 0; i < *n; i++) { 
	 tmp = sqrt(w[i]);
	 wy[i] = y[i] * tmp;

	 for (int j = 0; j < *p; j++) 
	    wx[*n * j + i] = x[*n * j + i] * tmp;
      }

      // compute the (weighted) least squares estimate (LAPACK::dgels),
      // solves minimize |B - A*X| for X (using QR factorization)
      F77_CALL(dgels)("N", // TRANS 
	 n,		   // M, dimension
	 p,		   // N, dimension
	 &int_1,	   // NRHS, no. of. columns of array B
	 wx,		   // on entry: A, array[M, N]; on return: QR factor
	 n,		   // LDA, dimension
	 wy,		   // on entry: B, array[LDB, NRHS], on return: least 
			   // squares coefficients (rows 1:N) 
	 n,		   // LDB, dimension
	 work,		   // WORK, array[LWORK]
	 lwork,		   // LWORK, dimension
	 &info_dgels);	   // INFO

      // dgels is not well suited as a rank-revealing procedure; i.e., INFO < 0
      // iff a diagonal element of the R matrix is exactly 0. This is not
      // helpful; hence, we check the diagonal elements of R separately and 
      // issue and error flag if any(abs(diag(R))) is close to zero  
      for (int i = 0; i < *p; i++) {
	 if (fabs(wx[(*n + 1) * i]) < sqrt(DBL_EPSILON)) {
	    *info = 1;
	    return;
	 }
      }

      // retrieve 'betacoefficients' (from array 'wy')
      Memcpy(beta0, wy, *p); 

      // compute the residuals (BLAS::dgemv): y = alpha*A*x + beta*y
      const double double_minus1 = -1.0, double_1 = 1.0; 
      Memcpy(resid, y, *n); 
      F77_CALL(dgemv)("N", // TRANS
	 n,		   // M, dimension
	 p,		   // N, dimension
	 &double_minus1,   // ALPHA, scalar
	 x,		   // A, array [LDA = M, N]
	 n,		   // LDA, dimension
	 beta0,		   // X, array
	 &int_1,	   // INCX 
	 &double_1,	   // BETA, scalar
	 resid,		   // Y, array (on return: modified Y)
	 &int_1);	   // INCY
   }
} 

/*****************************************************************************\
|*  psi functions (in fact, wgt-functions)				     *|
\*****************************************************************************/
static inline void robweight(double *resid, double *robwgt, double *xwgt, 
   double *k, double *scale, int *n, int *psi, int *type)
{
   // pre-treatment for Schweppe type GM
   if (*type == 2) {   
      for (int i = 0; i < *n; i++) 
	 resid[i] *= xwgt[i];	       // note: xwgt = 1 / xwgt (in rwlslm)
   }

   // M- and GM-estimators
   double z; 
   switch (*psi) { 
      case 0: // Huber weight 
	 for (int i = 0; i < *n; i++) { 
	    z = fabs(resid[i] / *scale); 
	    robwgt[i] *= z >= *k ? *k / z : 1.0; 
	 } 
	 break; 
      case 1: // asymmetric Huber weight 
	 for (int i = 0; i < *n; i++) { 
	    z = resid[i] / *scale; 
	    robwgt[i] *= fabs(z) >= *k ? *k / z : 1.0; 
	 } 
	 break; 
      case 2: // Tukey biweight 
	 for (int i = 0; i < *n; i++) {
	    z = resid[i] / *scale; 
	    robwgt[i] *= fabs(z) >= *k ? 0.0 : _POWER2( 1.0 - _POWER2(z / *k)); 
	 } 
	 break; 
   } 
}

/*****************************************************************************\
|*  weighted trimmed mean (scalar)					     *| 
|*                                                                           *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    lo	  lower bound [0, 1]					     *|
|*    hi	  upper bound [0, 1] s.t. lo < hi			     *|
|*    mean	  on return: weighted trimmed mean			     *|
|*    n		  dimension						     *|
\*****************************************************************************/
void wtrimmedmean(double *x, double *w, double *lo, double *hi, double *mean, 
   int *n)
{
   double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;   

   // quantiles   
   wquantile(x, w, n, lo, &quantile_lo); 
   wquantile(x, w, n, hi, &quantile_hi); 

   // trimmed mean 
   for (int i = 0; i < *n; i++) {
      if (quantile_lo <= x[i] && x[i] <= quantile_hi) {
	 sum_x += x[i] * w[i];
	 sum_w += w[i];
      }
   }

   if (sum_w > DBL_EPSILON) 
      *mean = sum_x / sum_w; 
   else {
      *mean = 0.0;
      error("Error: trimmed mean: division by zero\n");
   }
}

/*****************************************************************************\
|*  weighted winsorized mean (scalar)					     *| 
|*                                                                           *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    lo	  lower bound [0, 1]					     *|
|*    hi	  upper bound [0, 1] s.t. lo < hi			     *|
|*    mean	  on return: weighted trimmed mean			     *|
|*    n		  dimension						     *|
\*****************************************************************************/
void wwinsorizedmean(double *x, double *w, double *lo, double *hi, 
   double *mean, int *n)
{
   double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;   

   // quantiles   
   wquantile(x, w, n, lo, &quantile_lo); 
   wquantile(x, w, n, hi, &quantile_hi); 

   // winsorized mean 
   for (int i = 0; i < *n; i++) {
      if (x[i] < quantile_lo) 
	 sum_x += quantile_lo * w[i];
      else {
	 if (x[i] < quantile_hi)
	    sum_x += x[i] * w[i];
	 else
	    sum_x += quantile_hi * w[i];
      } 
      sum_w += w[i];
   }
   
   *mean = sum_x / sum_w;
}

/*****************************************************************************\
|*  weighted k one-sided winsorized mean (scalar)			     *| 
|*                                                                           *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    k		  k-th largest element (zero index)			     *|
|*    mean	  on return: weighted trimmed mean			     *|
|*    n		  dimension						     *|
|*    prob	  on return: estimated probability			     *|
\*****************************************************************************/
void wkwinsorizedmean(double *x, double *w, int *k, double *mean, 
   int *n, double *prob)
{
   double sum_xw = 0.0, below_sum_w = 0.0, above_sum_w = 0.0;

   // determine k-th largest element (it puts element k into its final 
   // position in the sorted array) 
   wselect0(x, w, 0, *n - 1, *k);

   // partial sums of x[0..k] * w[0..k] and w[0..k]  
   for (int i = 0; i < *k + 1; i++) {
      sum_xw += w[i] * x[i];
      below_sum_w += w[i];
   }

   // partial sum of w[(k+1)..(n-1)]
   for (int i = *k + 1; i < *n; i++)
      above_sum_w += w[i];

   sum_xw += above_sum_w * x[*k]; 
   
   *mean = sum_xw / (below_sum_w + above_sum_w);   // winsorized mean

   *prob = below_sum_w / (below_sum_w + above_sum_w); // estimate of prob
}

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
   double *x_wins;

   if (*n == 1) {
      *loc = *x;
      *scale = 0.0;
   } else {
      // initialize location (weighted median)
      double p50 = 0.5;
      wquantile(x, w, n, &p50, &loc0);

      // initialize variable for winsorized x-variable 
      x_wins = (double*) Calloc(*n, double);

      // initialize scale estimate by (weighted) interquartile range (IQR is
      // normalized to be a Fisher consistent estimate of the scale at the 
      // Gaussian distr.)
      double p25 = 0.25, x25 = 0.0;
      double p75 = 0.75, x75 = 0.0; 
      wquantile(x, w, n, &p25, &x25); 
      wquantile(x, w, n, &p75, &x75); 
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

/*****************************************************************************\
|*  euclidean norm							     *|
\*****************************************************************************/
static inline double euclidean_norm(const double *x, const double *y,
   const int p)
{
   double s = 0.0;
   for (int i = 0; i < p; i++) 
      s += _POWER2(x[i] - y[i]);

   return sqrt(s);
}

/*****************************************************************************\
|*  weighted median of the absolute deviations from the weighted median; the *|
|*  mad is normalized to be a consistent estimator of scale at the Gaussian  *|
|*  core model (i.e., times the constant 1.4826)	    		     *|
\*****************************************************************************/
static inline double wmad(double *x, double *w, int n, double constant)
{
   double med, mad, prob = 0.5;
   double *absdev;
   absdev = (double*) Calloc(n, double);

   wquantile(x, w, &n, &prob, &med); 

   // compute absolute deviation from the weighted median
   for(int i = 0; i < n; i++)
      absdev[i] = fabs(x[i] - med);

   wquantile(absdev, w, &n, &prob, &mad); 

   Free(absdev);
   return constant * mad;
}

#undef _WGT_HUBER
#undef _POWER2
#undef _MIN
#undef _MAX 

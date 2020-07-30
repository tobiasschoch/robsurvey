/*****************************************************************************\
|* sctbase								     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sct library							     *|
|* SUBEJCT  basic statistics functions					     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), January 19, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include "robsurvey.h"

/* some macros */
#define _WGT_HUBER(_x, _k) ((fabs(_x) >= _k) ? _k / fabs(_x) : 1.0)
#define _WGT_HUBERasym(_x, _k) ((fabs(_x) >= _k) ? _k / _x : 1.0)
#define _PSI_HUBER(_x, _k) ((_x <= -_k) ? -_k : ((_x < _k) ? _x : _k))
#define _PSI_PRIME_HUBER(_x, _k) ((fabs(_x) <= _k) ? 1.0 : 0.0) 
#define _POWER2(_x) ((_x) * (_x))
#define _MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#define _MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))

/*declaration of 'local' functions (inline imperative is GCC specific)*/
static inline void fitwls(double*, double*, double*, double*, double*, int*, 
   int*, int*)  __attribute__((always_inline));
static inline double kappa_huber(const double) __attribute__((always_inline));
static inline double euclidean_norm(const double*, const double*, int)
    __attribute__((always_inline));
static inline double wmad(double*, double*, int)
   __attribute__((always_inline));

/*****************************************************************************\
|*  Weighted least squares estimate			                     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    X           vectorized design matrix, array[n * p]                     *|
|*    y           response vector, array[n]                                  *|
|*    w           weights vector, array[n]                                   *|
|*    resid       on return: residuals vector, array[n]                      *|
|*    beta0       on return: coefficient vector, array[p]                    *|
|*    n, p        (array dimensions)                                         *|
|*    lwork	  size of array 'work' (if < 0, then 'dgels' determines      *|
|*                optimal size)                                              *|
|*                                                                           *|
|*  DEPENDENCIES                                                             *|
|*    dgels      LAPACK                                                      *|
|*    dgemv      LAPACK                                                      *|
|*                                                                           *|
|*   NOTE: if lwork < 0, rsaefit_wls determines the optimal size of array    *|
|* 	   work and returns it (lwork)	                                     *|
\*****************************************************************************/
static inline void fitwls(double *x, double *y, double *w, double *resid, 
			  double *beta0, int *ptrn, int *ptrp, int *ptrlwork)
{
   // define constants for the call of 'dgels'
   const int int_1 = 1;
   int info = 1;
   double *work, *wy, *wx, tmp;

   // STEP 0: compute optimal size of array work
   if (*ptrlwork < 0) {
      work = (double*) Calloc(*ptrn, double);
   } else {
      work = (double*) Calloc(*ptrlwork, double);
   }
   wy = (double*) Calloc(*ptrn, double);
   wx = (double*) Calloc(*ptrn * *ptrp, double);

   // pre-multiply the design matrix and the response vector by sqrt(w) 
   for (int i = 0; i < *ptrn; i++) { 
      tmp = sqrt(w[i]);
      wy[i] = y[i] * tmp;
      for (int j = 0; j < *ptrp; j++) {
	 wx[*ptrn * j + i] = x[*ptrn * j + i] * tmp;
      }
   }

   // compute the (weighted) least squares estimate (LAPACK::dgels),
   // solves minimize |B - A*X| for X (using QR factorization)
   F77_CALL(dgels)("N",	// TRANS 
      ptrn,		// M, dimension
      ptrp,		// N, dimension
      &int_1,		// NRHS, no. of. columns of array B
      wx,		// on entry: A, array[M, N]; on return: QR factor
      ptrn,		// LDA, dimension
      wy,		// on entry: B, array[LDB, NRHS], on return: least 
			// squares coefficients (rows 1:N) 
      ptrn,		// LDB, dimension
      work,		// WORK, array[LWORK]
      ptrlwork,		// LWORK, dimension
      &info);		// INFO

   // STEP 1: compute weighted ls
   if( *ptrlwork < 0 ) {
      *ptrlwork = (int) work[0]; // optimal value of 'lwork'
   } else {
      const double double_minus1 = -1.0, double_1 = 1.0; 

      // through an error if info is not zero  
      if (info > 0){
	 Rprintf("wls: dgels: not full rank\n");
      }
      if (info < 0){
	 Rprintf("wls: dgels: illegal value of argument %d\n", info);
      }

      // STEP 2: obtain 'betacoefficients' (from array 'wy')
      Memcpy(beta0, wy, *ptrp); 

      // STEP 3: compute the residuals (BLAS::dgemv): y = alpha*A*x + beta*y
      Memcpy(resid, y, *ptrn); 
      F77_CALL(dgemv)("N", // TRANS
	 ptrn,		   // M, dimension
	 ptrp,		   // N, dimension
	 &double_minus1,   // ALPHA, scalar
	 x,		   // A, array [LDA = M, N]
	 ptrn,		   // LDA, dimension
	 beta0,		   // X, array
	 &int_1,	   // INCX 
	 &double_1,	   // BETA, scalar
	 resid,		   // Y, array (on return: modified Y)
	 &int_1);	   // INCY
   }
 
   Free(work); Free(wy); Free(wx);
} 

/*****************************************************************************\
|*  weighted trimmed mean (scalar)					     *| 
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    ptrlo	  lower bound [0, 1]					     *|
|*    ptrhi	  upper bound [0, 1] s.t. ptrlo < ptrhi			     *|
|*    ptrmean	  on return: weighted trimmed mean			     *|
|*    ptrn	  dimension						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    wquantile		                                                     *|
\*****************************************************************************/
void wtrimmedmean(double *x, double *w, double *ptrlo, double *ptrhi, 
		  double *ptrmean, int *ptrn)
{
   double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;   

   // quantiles   
   wquantile(x, w, ptrn, ptrlo, &quantile_lo); 
   wquantile(x, w, ptrn, ptrhi, &quantile_hi); 

   // trimmed mean 
   for (int i = 0; i < *ptrn; i++){
      if (quantile_lo <= x[i] && x[i] <= quantile_hi){
	 sum_x += x[i] * w[i];
	 sum_w += w[i];
      }
   }

   if (sum_w > 0.0){
      *ptrmean = sum_x / sum_w; 
   }else{
      *ptrmean = 0.0;
      Rprintf("Trimmed mean: Division by zero\n");
   }
}

/*****************************************************************************\
|*  weighted winsorized mean (scalar)					     *| 
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    ptrlo	  lower bound [0, 1]					     *|
|*    ptrhi	  upper bound [0, 1] s.t. ptrlo < ptrhi			     *|
|*    ptrmean	  on return: weighted trimmed mean			     *|
|*    ptrn	  dimension						     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    wquantile		                                                     *|
\*****************************************************************************/
void wwinsorizedmean(double *x, double *w, double *ptrlo, double *ptrhi, 
		     double *ptrmean, int *ptrn)
{
   double quantile_lo, quantile_hi, sum_w = 0.0, sum_x = 0.0;   

   // quantiles   
   wquantile(x, w, ptrn, ptrlo, &quantile_lo); 
   wquantile(x, w, ptrn, ptrhi, &quantile_hi); 

   // winsorized mean 
   for (int i = 0; i < *ptrn; i++){
      if (x[i] < quantile_lo){
	 sum_x += quantile_lo * w[i];
      }else{
	 if (x[i] < quantile_hi){
	    sum_x += x[i] * w[i];
	 } else{
	    sum_x += quantile_hi * w[i];
	 }
      } 
      sum_w += w[i];
   }
   
   *ptrmean = sum_x / sum_w;
}

/*****************************************************************************\
|*  weighted k one-sided winsorized mean (scalar)			     *| 
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    x		  data, array[n]					     *|
|*    w		  weights, array[n]					     *|
|*    k		  k-th largest element (zero index)			     *|
|*    ptrmean	  on return: weighted trimmed mean			     *|
|*    ptrn	  dimension						     *|
|*    prob	  on return: estimated probability			     *|
|*                                                                           *|
|*  DEPENDENCIES							     *|
|*    wselect0		                                                     *|
\*****************************************************************************/
void wkwinsorizedmean(double *x, double *w, int *ptrk, double *ptrmean, 
   int *ptrn, double *prob)
{
   double sum_xw = 0.0, below_sum_w = 0.0, above_sum_w = 0.0;

   // determine k-th largest element (it puts element k into its final 
   // position in the sorted array) 
   wselect0(x, w, 0, *ptrn - 1, *ptrk);

   // partial sums of x[0..k] * w[0..k] and w[0..k]  
   for (int i = 0; i < *ptrk + 1; i++){
      sum_xw += w[i] * x[i];
      below_sum_w += w[i];
   }

   // partial sum of w[(k+1)..(n-1)]
   for (int i = *ptrk + 1; i < *ptrn; i++){
      above_sum_w += w[i];
   }

   sum_xw += above_sum_w * x[*ptrk]; 
   
   *ptrmean = sum_xw / (below_sum_w + above_sum_w);   // winsorized mean

   *prob = below_sum_w / (below_sum_w + above_sum_w); // estimate of prob
}

/*****************************************************************************\
|* huber mean and scale estimator: proposal 2 estimator (using weights;	     *| 
|* the same tuning constant 'k' used for location and scale; the scale	     *|
|* is normalized to be a Fisher consistent estimator at the Gaussian	     *| 
|* core model)								     *|
|*									     *|
|* PARAMETERS								     *|
|*    x	       array							     *|
|*    w	       array of weights				 		     *|
|*    robwgt   array of robustness weights (on return)	  		     *|
|*    ptrk     robustness tuning constant (location and scale)		     *|
|*    ptrloc   robust location (on return)				     *|
|*    ptrscale robust scale (on return)					     *|
|*    ptrn     array dimension						     *|
|*    ptrmaxit maximum number of iteration (termination rule); on return:    *|
|*	       effective number of iterations				     *|
|*    ptrtol   numerical convergence tolerance for iterative refinement	     *|
|*	       iterations						     *|
|*									     *|
|* DEPENDENCIES								     *|
|*    wquantile	  							     *|
\*****************************************************************************/
void huberm(double *x, double *w, double *robwgt, double *ptrk, 
	    double *ptrloc, double *ptrscale, int *ptrn, int *ptrmaxit, 
	    const double *ptrtol)
{
   int iter;
   double loc0, scale0, tmp, tmp_loc, tmp_scale, kappa, wtotal; 
   double *x_wins;
   if(*ptrn == 1){
      *ptrloc = *x;
      *ptrscale = 0.0;
   }else{
      // initialize location (weighted median)
      double p50 = 0.5;
      wquantile(x, w, ptrn, &p50, &loc0);

      // initialize variable for winsorized x-variable 
      x_wins = (double*) Calloc(*ptrn, double);

      // initialize scale estimate by (weighted) interquartile range (IQR is
      // normalized to be a Fisher consistent estimate of the scale at the 
      // Gaussian distr.)
      double p25 = 0.25, x25 = 0.0;
      double p75 = 0.75, x75 = 0.0; 
      wquantile(x, w, ptrn, &p25, &x25); 
      wquantile(x, w, ptrn, &p75, &x75); 
      scale0 = (x75 - x25) / 1.349;
 
      // stop if IQR is zero or negative
      if (scale0 <= 0.0){
	 *ptrloc = x[1];
      }else{
	 // compute the weight total
	 wtotal = 0.0;
	 for (int i = 0; i < *ptrn; i++) wtotal += w[i];

	 // compute consistency correction term
	 kappa = kappa_huber(*ptrk);

	 // loop   
	 for (iter = 0; iter < *ptrmaxit; iter ++){

	    // update location and scale 
	    tmp_loc = 0.0; tmp_scale = 0.0;
	    tmp = *ptrk * scale0;
	    for (int i = 0; i < *ptrn; i++){
	       // winsorized x-variable
	       x_wins[i] = _MIN(_MAX(loc0 - tmp, x[i]), loc0 + tmp);
	       // location step
	       tmp_loc += w[i] * x_wins[i];
	       // scale
	       tmp_scale += w[i] * _POWER2(x_wins[i] - loc0);
	    } 
	    *ptrloc = tmp_loc / wtotal;
	    *ptrscale = tmp_scale / wtotal;
	    // normalize the variance/ scale  
	    *ptrscale = sqrt(*ptrscale / kappa);

	    // termination rule
	    if (fabs(*ptrloc - loc0) < *ptrtol * scale0 &&	
	       fabs(*ptrscale / scale0 - 1) < *ptrtol){
	       break;
	    }else{
	       loc0 = *ptrloc;
	       scale0 = *ptrscale;
	    }
	 }
	 // on return: maxit = effective number of iterations
	 *ptrmaxit = iter;

	 // compute robustness weights 
	 for (int i = 0; i < *ptrn; i++){ 
	    robwgt[i] = _WGT_HUBER((x[i] - *ptrloc) / *ptrscale, *ptrk); 
	 }

	 Free(x_wins);
      }
   }
}

/*****************************************************************************\
|*                                                                           *|
|*  This function computes the multiplicative consistency correction term    *|
|*  for the Huber 'Proposal 2' type of scale estimator.                      *|
|*                                                                           *|
\*****************************************************************************/
static inline double kappa_huber(const double k) 
{
   double pdf_k, cdf_k; 
   if (k < 10.0) {
      pdf_k = exp(-_POWER2(k) / 2.0) / 2.5066282746;
      /*lower tail of cdf; i.e., 1 - Phi(.)*/
      cdf_k = 0.5 - 0.5 * erf(k / 1.4142135623);
      return 1.0 - 2.0 * (k * pdf_k + (1.0 - k * k) * cdf_k);
   } else {
      return 1.0;
   }
}

/*****************************************************************************\
|*  rwlslm: iteratively reweighted least squares			     *|
|*                                                                           *|
|*  PARAMETERS								     *|
|*    x		  vectorized design matrix, array[n * p]		     *|
|*    y		  response vector, array[n]	   			     *| 
|*    w		  weights, array[n]					     *|
|*    resid	  on return: residuals vector, array[n]			     *|
|*    robwgt	  on return: robustness weights, array[n]		     *|
|*    ptrn	  (array dimensions)					     *|
|*    ptrp	  (array dimensions)					     *|
|*    ptrk	  robustness tuning constant				     *|
|*    beta0	  on return: coefficient vector, array[p]		     *|
|*    ptrscale	  on return: scale estimate				     *|
|*    ptrmaxit	  max iterations (on input); iterations (on return)	     *|
|*    ptrtol	  numerical tolerance criterion (stoping rule in rwls	     *|
|*		  updating rule)					     *|
|*    ptrpsi	  psi-function (0 = Huber, 1 = asymmetric Huber	     	     *|
|*    ptrpsi2	  expected value of psi^2, i.e.: sum(w * psi^2)		     *|
|*    ptrpsiprime expected value of psi', i.e..: sum(w * psi')  	     *|
|*                                                                           *|
|*  DEPENDENCIES                                                             *|
|*    sctbase: fitwls                                                        *|
|*    wmad                                                                   *|
\*****************************************************************************/
void rwlslm(double *x, double *y, double *w, double *resid, double *robwgt, 
	    int *ptrn, int *ptrp, double *ptrk, double *beta0, 
	    double *ptrscale, int *ptrmaxit, double *ptrtol, int *ptrpsi, 
	    double *ptrpsi2, double *ptrpsiprime)
{
   int iterations = 0, converged = 0;
   double scale, tmp;
   double *beta_new, *ui_weight;

   // STEP 1: initialize beta by weighted least squares
   beta_new = (double* ) Calloc(*ptrp, double);
   ui_weight = (double *) Calloc(*ptrn, double);

   // determine optimal size of array 'work' in 'fitwls'
   int lwork = -1;
   fitwls(x, y, w, resid, beta0, ptrn, ptrp, &lwork);

   // compute initial wls fit 
   // first, we compute a modified weight 
   fitwls(x, y, w, resid, beta0, ptrn, ptrp, &lwork);

   // STEP 2: initialize scale estimate by weighted MAD
   scale = wmad(resid, w, *ptrn);

   // STEP 3: irls updating
   while (!converged && ++iterations < *ptrmaxit){

      // STEP 3.1: update beta
      if (*ptrpsi == 0){ // Huber psi
	 for (int i = 0; i < *ptrn; i++){
	    ui_weight[i] = w[i] * _WGT_HUBER(resid[i] / scale, *ptrk);
	 }
      } else{ // asymmetric Huber psi
	 for (int i = 0; i < *ptrn; i++){
	    ui_weight[i] = w[i] * _WGT_HUBERasym(resid[i] / scale, *ptrk);
	 }
      }
      // STEP 3.1: update beta and residuals
      fitwls(x, y, ui_weight, resid, beta_new, ptrn, ptrp, &lwork);
 
      // STEP 3.2: update scale
      scale = wmad(resid, w, *ptrn);

      // check for convergence
      converged = (euclidean_norm(beta0, beta_new, *ptrp) < *ptrtol * 
	 scale) ? 1: 0;

      // prepare the next while 'run'
      Memcpy(beta0, beta_new, *ptrp);
   }
   *ptrmaxit = (converged) ? iterations : 0;

   // compute robustness weights (using MAD estimate of scale and beta) 
   if (*ptrpsi == 0){ // Huber psi
      for (int i = 0; i < *ptrn; i++){
	 robwgt[i] = _WGT_HUBER(resid[i] / scale, *ptrk);

//FIXME: psi^2 and psi'

      }
   } else { // asymmetric Huber psi
      for (int i = 0; i < *ptrn; i++){
	 robwgt[i] = _WGT_HUBERasym(resid[i] / scale, *ptrk);
      }
   }

   // compute 'Proposal 2' estimate of scale 
   tmp = 0.0;
   for (int i = 0; i < *ptrn; i++){  
      tmp += w[i] * _POWER2(robwgt[i] * resid[i]); 
   }
   *ptrscale = sqrt(tmp / (*ptrn * kappa_huber(*ptrk)));

   // empirical estimate of E(psi^2) 
   *ptrpsi2 = 0.0;

   // empirical estimate of E(psi')
   *ptrpsiprime = 0.0;

   Free(beta_new); Free(ui_weight);
}

/*****************************************************************************\
|*  euclidean norm							     *|
|*                                                                           *|
|*  PARAMETERS								     *|
|*    x	 array[n]							     *|
|*    y	 array[n]							     *|
|*    n	 dimension							     *|
\*****************************************************************************/
static inline double euclidean_norm(const double *x, const double *y,
				   const int n)
{
   double s = 0.0;
   for (int i = 0; i < n; i++) {
      s += _POWER2(x[i] - y[i]);
   }
   return(sqrt(s));
}

/*****************************************************************************\
|*  weighted median of the absolute deviations from the weighted median; the *|
|*  mad is normalized to be a consistent estimator of scale at the Gaussian  *|
|*  core model (i.e., times the constant 1.4826)	    		     *|
|*                                                                           *|
|*  PARAMETERS								     *|
|*    x	 data, array[n]							     *|
|*    w	 weights, array[n]						     *|
|*    n	 dimension							     *|
\*****************************************************************************/
static inline double wmad(double *x, double *w, int n)
{
   double med, mad, prob = 0.5;
   double *absdev;
   absdev = (double*) Calloc(n, double);

   wquantile(x, w, &n, &prob, &med); 

   // compute absolute deviation from the weighted median
   for(int i = 0; i < n; i++){
      absdev[i] = fabs(x[i] - med);
   }

   wquantile(absdev, w, &n, &prob, &mad); 

   Free(absdev);
   return(1.4826 * mad);
}

#undef _WGT_HUBER
#undef _POWER2
#undef _MIN
#undef _MAX
#undef _PSI_HUBER
#undef _PSI_PRIME_HUBER



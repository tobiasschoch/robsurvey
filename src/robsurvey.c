/* Functions to compute weighted (generalized) regression M-estimators,  
   winsorized, and trimmed estimators of location

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

#include "robsurvey.h"

// some macros
#define _POWER2(_x) ((_x) * (_x))

// declaration of 'local' functions (inline imperative is GCC specific)
static inline void fitwls(double*, double*, double*, double*, double*, double*, 
   double*, int*, int*, double*, int*, int*) __attribute__((always_inline));
static inline double euclidean_norm(const double*, const double*, int)
    __attribute__((always_inline));
static inline double wmad(double*, double*, double*, int, double)
   __attribute__((always_inline));
static inline void robweight(double*, double*, double*, double*, double*, 
   double*, double*, int*, int*, int*) __attribute__((always_inline));

double huber_psi(double, double); 
double huber_psi_prime(double, double); 
static inline double huber_wgt(double, double) __attribute__((always_inline));

double huber_psi_asym(double, double); 
double huber_psi_prime_asym(double, double); 
static inline double huber_wgt_asym(double, double) 
   __attribute__((always_inline));

double tukey_psi(double, double);
double tukey_psi_prime(double, double);
static inline double tukey_wgt(double, double) __attribute__((always_inline));

static inline void inverse_qr(double*, double*, double*, int*, int*, int*, int)
   __attribute__((always_inline));

/*****************************************************************************\
|*  rwlslm: iteratively reweighted least squares			     *|
|*                                                                           *|
|*    x		  vectorized design matrix, array[n * p]		     *|
|*    y		  response vector, array[n]	   			     *| 
|*    w		  weights, array[n] (on return: a vector of ones)	     *|
|*    resid	  on return: residuals vector, array[n]			     *|
|*    robwgt	  on return: robustness weights, array[n]		     *|
|*    xwgt	  weight in design space, GM-estimator, array[n]	     *|
|*    n		  (array dimensions)					     *|
|*    p		  (array dimensions)					     *|
|*    k		  robustness tuning constant				     *|
|*    beta0	  on return: coefficient vector, array[p]		     *|
|*    scale	  on return: normalized MAD				     *|
|*    tol	  numerical tolerance criterion (stoping rule in rwls	     *|
|*		  updating rule)					     *|
|*    maxit	  max iterations (on input); iterations (on return)	     *|
|*    psi	  0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight  	     *|
|*    type	  0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.	     *|
\*****************************************************************************/
void rwlslm(double *x, double *y, double *w, double *resid, double *robwgt, 
   double *xwgt, int *n, int *p, double *k, double *beta0, double *scale, 
   double *tol, int *maxit, int *psi, int *type)
{
   int info = 0, iterations = 0, converged = 0;
   double mad_const = 1.482602;	 
   double *work_x, *work_y, *work, *beta_new, *w_mallows;

   // STEP 0 preparations 
   beta_new = (double*) Calloc(*p, double);
   work_x = (double*) Calloc(*n * *p, double);	// work array
   work_y = (double*) Calloc(*n, double);	// work array 
   w_mallows = (double*) Calloc(*n, double);	// only used for Mallows GM 

   // determine optimal size of array 'work' and allocate it 
   int lwork = -1; 
   fitwls(x, work_x, y, work_y, w, resid, beta0, n, p, work_x, &lwork, &info);
   work = (double*) Calloc(lwork, double);

   // STEP 1: initialize beta by weighted least squares
   fitwls(x, work_x, y, work_y, w, resid, beta0, n, p, work, &lwork, &info);
   if (info > 0){
      error("The design matrix is rank deficient (or nearly so)\n");
      *maxit = 0; 
      return;
   }

   // preparations for GM-estimators
   if (*type == 1) {					       // Mallows 
      Memcpy(w_mallows, w, *n);
      for (int i = 0; i < *n; i++)
	 w[i] *= xwgt[i];  // modify w to include xwgt 
     
      // Mallows type has special normalization constant (see 'mallows.c')
      mad_const = mallows_mad_normalization(xwgt, n);
   }
   if (*type == 2) {					       // Schweppe  
      for (int i = 0; i < *n; i++) {
   	 if (fabs(xwgt[i]) < DBL_EPSILON)
	    xwgt[i] = 0.0;	
	 else 
	    xwgt[i] = 1.0 / xwgt[i]; // xwgt turned into multiplicative weight
      }
   }

   // STEP 2: initialize scale by weighted MAD (ignore that Mallows is special)
   *scale = wmad(resid, w, work_y, *n, mad_const);
   if (*scale < DBL_EPSILON) {
      error("The estimate of scale is zero (or nearly so)\n");
      *maxit = 0;
      return;
   }

   // STEP 3: irls updating
   while (!converged && ++iterations < *maxit) {
      // compute robwgt (i.e. total wgt = sampling wgt * xwgt * robwgt)
      robweight(resid, robwgt, xwgt, w, work_y, k, scale, n, psi, type);

      // update beta and residuals
      fitwls(x, work_x, y, work_y, robwgt, resid, beta_new, n, p, work, &lwork, 
	 &info);

      // update scale 
      if (*type == 1) {					       // Mallows GM 

	 // we work on work_y, which is a temporary copy of resid 
	 Memcpy(work_y, resid, *n);	
	 for (int i = 0; i < *n; i++)
	    work_y[i] *= sqrt(xwgt[i]);

	 *scale = wmad(work_y, w_mallows, work_x, *n, mad_const);
 
      } else						       // otherwise
	 *scale = wmad(resid, w, work_y, *n, mad_const); 

      // check for convergence
      converged = (euclidean_norm(beta0, beta_new, *p) < *tol * *scale) ? 1: 0;

      // prepare the next while run 
      Memcpy(beta0, beta_new, *p); 
   }
   *maxit = (converged) ? iterations : 0;

   // compute robwgt (without sampling weights) 
   for (int i = 0; i < *n; i++) 
      robwgt[i] /= w[i];

   Free(beta_new); Free(work_x); Free(work_y); Free(work); Free(w_mallows);
}

/*****************************************************************************\
|* cov_rwlslm: covariance matrix of the esimated regression coefficients     *|
|*									     *|
|*    resid	     residuals, array[n]				     *|
|*    x		     design matrix, array[n * p]; on return: the p * p cov   *|
|*		     matrix is stored in x[1..(p * p)]			     *|
|*    xwgt	     weights in design space, array[n]			     *|
|*    robwgt	     robustness weight, array[n]			     *|
|*    w		     sampling weight, array[n]				     *|
|*    k		     robustness tuning constant				     *|
|*    scale	     estimate of scale					     *|
|*    scale2	     on return: estimate of scale (proposal 2)		     *|
|*    n, p	     dimensions						     *|
|*    psi	     0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight     *|
|*    type	     0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.   *|
\*****************************************************************************/
void cov_rwlslm(double *resid, double *x, double *xwgt, double *robwgt, 
   double *w, double *k, double *scale, double *scale2, int *n, int *p, 
   int *psi, int *type)
{
   double tmp, sum_w = 0.0;
   double *work, *work_x, *work_y;
  
   work_x = (double*) Calloc(*n * *p, double);
   work_y = (double*) Calloc(*n, double);

   // determine lwork and allocate work: dgeqrf (used in inverse_qr)
   int lwork = -1, info;
   F77_CALL(dgeqrf)(n, p, x, n, work_x, work_y, &lwork, &info); 	
   lwork = (int) work_y[0]; 
   work = (double*) Calloc(lwork, double);

   // function ptrs (initialized, otherwise we get an "uninitialized" warning)
   double (*f_psi)(double, double) = huber_psi;		    
   double (*f_psiprime)(double, double) = huber_psi_prime; 

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
   }

   // Schweppe GM-est. ---------------------------------------------------------
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

	 for (int j = 0; j < *p; j++)		// x := sqrt(s_1 * w) o x
	    x[*n * j + i] *= sqrt(tmp * w[i]);    

	 work_x[i] = tmp2 / tmp;		// temporarily store s_2 / s_1
      }

      Memcpy(work_y, work_x, *n);		// temporarily store s_2 / s_1

      // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
      inverse_qr(x, work_x, work, n, p, &lwork, 1);	

      // pre-multiply Q by sqrt(s2 / s1) 
      for (int i = 0; i < *n; i++) {	  
	 tmp = sqrt(work_y[i]);    
	 for (int j = 0; j < *p; j++) 
	    x[*n * j + i] *= tmp;   // pre-multiply Q 
      }

      // B  := Q * R^{-T} (result -> x)
      double done = 1.0, dzero = 0.0; 
      F77_CALL(dtrmm)("R", "U", "T", "N", n, p, &done, work_x, p, x, n);

      // compute B^T * B := (x^T * W * W * x)^{-1}
      *scale2 = _POWER2(*scale) / (1.0 - (double)*p / sum_w);
      F77_CALL(dgemm)("T", "N", p, p, n, scale2, x, n, x, n, &dzero, work_x, p);

      *scale2 = *scale;

   // M-est. and Mallows GM-est. -----------------------------------------------
   } else {					

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

	 // correction factor (see Huber, 1981, p. 172-174) 
	 double kappa = 1.0 + (double)*p / sum_w * (Epsi_prime2 / 
	    _POWER2(Epsi_prime) - 1.0) * (double)*n / (double)(*n - 1);
	 *scale2 *= _POWER2(kappa);

	 // QR factorization of x (goal: R^{-1} * R^{-T} =: inverse of x^T * x)
	 for (int i = 0; i < *n; i++) {    
	    tmp = sqrt(w[i]); // pre-multiply x by sampling weight  		  
	    for (int j = 0; j < *p; j++) 
	       x[*n * j + i] *= tmp;      
	 }

	 inverse_qr(x, work_x, work, n, p, &lwork, 0);	

	 F77_CALL(dtrmm)("R", "U", "T", "N", p, p, scale2, work_x, p, work_x, p);

      // Mallows GM-est. --------------------------------------------
      } else {		

	 for (int i = 0; i < *n; i++) {    
	    tmp = sqrt(w[i] * xwgt[i]);   
	    for (int j = 0; j < *p; j++) 
	       x[*n * j + i] *= tmp;   // pre-multiply x 
	 }

	 // QR factorization: Q -> x; R^{-1} -> work_x[1..(p * p)]
	 inverse_qr(x, work_x, work, n, p, &lwork, 1);	

	 // pre-multiply Q by with sqrt(xwgt) 
	 for (int i = 0; i < *n; i++) {	  
	    tmp = sqrt(xwgt[i]);    
	    for (int j = 0; j < *p; j++) 
	       x[*n * j + i] *= tmp;   // pre-multiply Q 
	 }

	 // B  := Q * R^{-T} (result -> x)
	 double done = 1.0, dzero = 0.0; 
	 F77_CALL(dtrmm)("R", "U", "T", "N", n, p, &done, work_x, p, x, n);
    
	 // compute B^T * B := (x^T * W * W * x)^{-1}
	 F77_CALL(dgemm)("T", "N", p, p, n, scale2, x, n, x, n, &dzero, work_x, p);
      }
   }
   Memcpy(x, work_x, *p * *p);	       // put result in x[1..(p * p)]
   Free(work); Free(work_x); Free(work_y);
}

/*****************************************************************************\
|*    Inverse of R matrix and Q matrix of the QR factorization		     *|
|*                                                                           *|
|*    x		  on return; Q matrix, array[n * p]			     *|
|*    work_x	  on return: inv. R is on work_x[1..(p * p)], array[n * p]   *|
|*    work	  work array used for QR factorization, array[lwork]	     *|
|*    n, p, lwork dimensions						     *|
|*    qmatrix	  toogle whether Q matrix is computed : 0 = no; 1 = yes	     *|
|*									     *|
|*    NOTE: array x will be overwritten					     *|
\*****************************************************************************/
static inline void inverse_qr(double *x, double *work_x, double *work, 
   int *n, int *p, int *lwork, int qmatrix)
{
   int info = 1;					 // QR factoriz. of x
   int offset = _POWER2(*p);
   F77_CALL(dgeqrf)(n, p, x, n, work_x + offset, work, lwork, &info); 	
   if (info != 0)
      error("dgeqrf failed\n");

   for (int i = 0; i < *p * *p; i++)			 // prepare matrix R
      work_x[i] = 0.0;

   for (int i = 0; i < *p; i++)				 // extract matrix R
      for (int j = 0; j < i + 1; j++) 
	 work_x[j + i * *p] = x[j + i * *n];
      
   F77_CALL(dtrtri)("U", "N", p, work_x, p, &info);	 // inverse of R 
   if (info != 0)
      error("dtrtri failed\n");

   if (qmatrix) {
      F77_CALL(dorgqr)(n, p, p, x, n, work_x + offset,	 // extract matrix Q
         work, lwork, &info);
      if (info != 0)
         error("dorgqr failed\n");
   }
}

/*****************************************************************************\
|*  Weighted least squares estimate			                     *|
|*                                                                           *|
|*    x           vectorized design matrix, array[n * p]                     *|
|*    work_x	  array[n * p], on return: QR factor (dgels)		     *|
|*    y           response vector, array[n]                                  *|
|*    work_y	  array[n]						     *|
|*    w           weights vector, array[n]                                   *|
|*    resid       on return: residuals vector, array[n]                      *|
|*    beta0       on return: coefficient vector, array[p]                    *|
|*    n, p        (array dimensions)                                         *|
|*    work        work array, array[lwork] used for QR factorization         *|
|*    lwork	  size of array 'work' (if < 0, then 'dgels' determines      *|
|*                optimal size)                                              *|
|*    info	  on return: info on fitwls (error if != 0)		     *|
\*****************************************************************************/
static inline void fitwls(double *x, double *work_x, double *y, double *work_y, 
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
	 work_y[i] = y[i] * tmp;

	 for (int j = 0; j < *p; j++) 
	    work_x[*n * j + i] = x[*n * j + i] * tmp;
      }

      // compute the (weighted) least squares estimate (LAPACK::dgels),
      // solves minimize |B - A*X| for X (using QR factorization)
      F77_CALL(dgels)("N", n, p, &int_1, work_x, n, work_y, n, work, lwork, 
	 &info_dgels);	  

      // dgels is not well suited as a rank-revealing procedure; i.e., INFO < 0
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
      F77_CALL(dgemv)("N", n, p, &double_minus1, x, n, beta0, &int_1, &double_1,
	 resid,&int_1);	 
   }
} 

/*****************************************************************************\
|* Huber psi-, psi-prime- and wgt-functions				     *|
\*****************************************************************************/
double huber_psi(double x, double k)
{
   return (x <= -k) ? -k : ((x < k) ? x : k);
}

double huber_psi_prime(double x, double k)
{
   return (fabs(x) <= k) ? 1.0 : 0.0; 
}

static inline double huber_wgt(double x, double k)
{
   double z = fabs(x);
   if (z > DBL_EPSILON)
      return (z >= k) ? k / z : 1.0;
   else 
      return 0.0;
}

/*****************************************************************************\
|* Huber asymmetric psi-, psi-prime- and wgt-functions			     *|
\*****************************************************************************/
double huber_psi_asym(double x, double k)
{
   return (x <= k) ? x : k;
}

double huber_psi_prime_asym(double x, double k)
{
   return (x <= k) ? 1.0 : 0.0;
}

static inline double huber_wgt_asym(double x, double k)
{
   double z = fabs(x);
   if (z > DBL_EPSILON)
      return (x <= k) ? 1.0 : k / x;
   else
      return 0.0;
}

/*****************************************************************************\
|* Tukey biweigt psi-, psi-prime- and wgt-functions			     *|
\*****************************************************************************/
double tukey_psi(double x, double k)
{
   if (fabs(x) > k)
      return 0.0;
   else {
      double z = x / k;
      double u = 1.0 - _POWER2(z);
      return x * _POWER2(u);
   }
}

double tukey_psi_prime(double x, double k)
{
   if (fabs(x) > k)
	return 0.0;
    else {
      x /= k;
      double z = _POWER2(x);
      return (1.0 - z) * (1.0 - 5.0 * z);
    }
}

static inline double tukey_wgt(double x, double k)
{
   if (fabs(x) > k)
      return 0.0;
   else {
      double z = x / k;
      z = (1.0 - z) * (1.0 + z);
      return _POWER2(z);
    }
}

/*****************************************************************************\
|* Robustness wgt-functions						     *|
|*                                                                           *|
|*    resid	  residuals, array[n]					     *|
|*    robwgt	  on return: robustness weight, array[n]		     *|
|*    xwgt	  weight on the design space, array[n]			     *|
|*    w		  sampling weight, array[n]				     *|
|*    work_y	  work array, array[n]					     *|
|*    k		  robustness tuning constant				     *|
|*    scale	  estimate of scale					     *|
|*    n		  dimension						     *|
|*    psi	  0 = Huber, 1 = asymmetric Huber, 2 = Tukey biweight	     *|
|*    type	  0 = M-est., 1 = Mallows GM-est., 2 = Schweppe GM-est.	     *|
\*****************************************************************************/
static inline void robweight(double *resid, double *robwgt, double *xwgt, 
   double *w, double *work_y, double *k, double *scale, int *n, int *psi, 
   int *type)
{
   double *ptr;

   if (*type == 2) {		       // pre-treatment for Schweppe type GM
      Memcpy(work_y, resid, *n);

      for (int i = 0; i < *n; i++) 
	 work_y[i] *= xwgt[i];	       // note: xwgt := 1.0 / xwgt (in rwlslm)

      ptr = resid;
      resid = work_y;		       // resid -> work_y
   } 
 
   switch (*psi) {		       // M- and GM-estimators
      case 0: // Huber weight 
	 for (int i = 0; i < *n; i++) 
	    robwgt[i] = w[i] * huber_wgt(resid[i] / *scale, *k);
	 break; 

      case 1: // asymmetric Huber weight
	 for (int i = 0; i < *n; i++) 
	    robwgt[i] = w[i] * huber_wgt_asym(resid[i] / *scale, *k);
	 break; 

      case 2: // Tukey biweight weight 
	 for (int i = 0; i < *n; i++) 
	    robwgt[i] = w[i] * tukey_wgt(resid[i] / *scale, *k);
	 break; 
   }
   
   if (*type == 2)		       // undo resid -> work_w
      resid = ptr;      
}

/*****************************************************************************\
|*  weighted median of the absolute deviations from the weighted median; the *|
|*  mad is normalized to be a consistent estimator of scale at the Gaussian  *|
|*  core model (i.e., times the constant 1.4826)	    		     *|
|*									     *|
|*    x	       data, array[n]						     *|
|*    w	       weights, array[n]					     *|
|*    work     work array[n]						     *|
|*    n	       dimension						     *|
|*    constant normalization constant					     *|
\*****************************************************************************/
static inline double wmad(double *x, double *w, double *work, int n, 
   double constant)
{
   double med, mad, prob = 0.5;

   wquantile(x, w, &n, &prob, &med); 

   // compute absolute deviation from the weighted median
   for (int i = 0; i < n; i++)
      work[i] = fabs(x[i] - med);

   wquantile(work, w, &n, &prob, &mad); 
 
   return constant * mad;
}

/*****************************************************************************\
|*  euclidean norm							     *|
|*									     *|
|*    x	    array[p]							     *|
|*    y	    array[p]							     *|
|*    p	    dimension							     *|
\*****************************************************************************/
static inline double euclidean_norm(const double *x, const double *y,
   const int p)
{
   double s = 0.0;
   for (int i = 0; i < p; i++) 
      s += _POWER2(x[i] - y[i]);

   return sqrt(s);
}

#undef _POWER2 

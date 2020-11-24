/* Implementation of the weighted BACON algorithm of Billor et al. (2000), 
   respectively, Béguin and Hulliger (2008) 

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

   Billor N, Hadi AS, Vellemann PF (2000). BACON: Blocked Adaptative 
      Computationally efficient Outlier Nominators. Computational Statistics 
      and Data Analysis 34, pp. 279-298.
   Béguin C, Hulliger B (2008). The BACON-EEM Algorithm for Multivariate 
      Outlier Detection in Incomplete Survey Data. Survey Methodology 34, 
      pp. 91-103.
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "wquantile.h" 

// macros
#define _RANK_TOLERANCE 1.0e-8
#define _POWER2(_x) ((_x) * (_x))

// prototype of local function
static inline double cutoffval(double, int, int, int) 
   __attribute__((always_inline));
static inline void initialsubset(double*, double*, double*, int*, int*, int*, 
   int*, int*) __attribute__((always_inline));
static inline void weightedscatter(double*, double*, double*, double*, double*, 
   int*, int*) __attribute__((always_inline));
static inline void weightedmean(double*, double*, double*, int*, int*) 
   __attribute__((always_inline));
static inline void mahalanobis(double*, double*, double*, double*, double*, 
   int*, int*) __attribute__((always_inline));

/*****************************************************************************\
|*  Chi-squared cutoff value (p-dimensional, sample size n	  	     *|
|*    alpha	  scalar (0 <= alpha <= 1)		                     *|
|*    n		  scalar						     *|
|*    p		  scalar	 					     *|
|*                                                                           *|
|*    dependency: qchisq <R_ext/Utils.h>                                     *|
\*****************************************************************************/
static inline double cutoffval(double alpha, int n, int k, int p)
{
   double h = ((double)n + (double)p + 1.0) / 2.0;
   double chr = fmax(0.0, (h - (double)k) / (h + (double)k));
   double cnp = 1.0 + ((double)p + 1.0) / ((double)n - (double)p) + 
      2.0 / ((double)n - 1.0 - 3.0 * (double)p);
   double cnpr = cnp + chr;
   double cutoff = cnpr * sqrt(qchisq(alpha / (double)n, p, 0, 0));

   return cutoff;
}

/*****************************************************************************\
|* initialsubset	   		     				     *|
|*    x		  array[n, p]				                     *|
|*    w		  on exit: elements not in subset have w = 0, array[n]       *|
|*    dist	  array[n]						     *|
|*    subset	  on exit: initial subset		                     *|
|*    subsetsize  on exit: size initial subset		                     *|
|*    n, p	  dimensions		  		                     *|
|*    verbose	  0 = quiet; 1 =  verbose				     *|
\*****************************************************************************/
static inline void initialsubset(double *x, double *w, double *dist, 
   int *subset, int *subsetsize, int *n, int *p, int *verbose)
{
   int *iarray;
   double *dist_sorted, *x_sorted, *w_sorted, *w_sortedcpy, *center, *scatter, 
      *work;

   center = (double*) Calloc(*p, double);
   scatter = (double*) Calloc(*p * *p, double);
   work = (double*) Calloc(*n * *p, double);

   // STEP 0: sort Mahalanobis distance from small to large (make copies 
   // because the arrays are sorted in place)
   dist_sorted = (double*) Calloc(*n, double);
   x_sorted = (double*) Calloc(*n * *p, double);
   w_sorted = (double*) Calloc(*n, double);
   Memcpy(dist_sorted, dist, *n);
   Memcpy(x_sorted, x, *n * *p);
   Memcpy(w_sorted, w, *n);

   // generate and populate 'iarray'
   iarray = (int*) Calloc(*n, int);

   for (int i = 0; i < *n; i++)
      iarray[i] = i;

   // sort dist (iarray is sorted along with dist_sorted)
   rsort_with_index(dist_sorted, iarray, *n);

   // use iarray to sort x (x_sorted) and w (w_sorted)
   for (int i = 0; i < *n; i++) {
      w_sorted[i] = w[iarray[i]];

      for (int j = 0; j < *p; j++)
	 x_sorted[i + *n * j] = x[iarray[i] + *n * j]; 
   }

   int m = (int)fmin(4.0 * (double)*p, (double)*n * 0.5); // subset size

   w_sortedcpy = (double*) Calloc(*n, double);
   Memcpy(w_sortedcpy, w_sorted, *n);

   // STEP 1 check if scatter matrix of the first 1:m observations has full rank 
   int rank;
   while (m < *n) {
      // set weights of observations (m+1):n to zero (note: we use int i = m,
      // because of C's zero indexing
      for (int i = m; i < *n; i++)
	 w_sortedcpy[i] = 0.0;	  

      weightedscatter(x_sorted, work, w_sortedcpy, center, scatter, n, p); 
      
      // Cholesky decomposition of scatter matrix 
      int info;
      F77_CALL(dpotrf)("L", p, scatter, p, &info);		
      if (info != 0) error("Initial subset: DPOTRF failed\n");

      // count diagonal elements of Cholesky factor > tol to determine the rank
      // note: rank revealing by Cholesky is cheap, but not numerically stable
      rank = 0;
      for (int i = 0; i < *p; i++)
	 rank += scatter[i * (*p + 1)] > _RANK_TOLERANCE ? 1 : 0;

      if (rank == *p) 
	 break;
      else 
	 error("Scatter matrix is rank defficient: subset is enlarged\n");	 

      m++;	
      Memcpy(w_sortedcpy, w_sorted, *n);
   }

   // STEP 2: generate initial subset    
   for (int i = 0; i < m; i++)	    
      subset[iarray[i]] = 1; 
   
   for (int i = m; i < *n; i++)	 // set weight = 0 when not in subset
      w[iarray[i]] = 0.0;  

   *subsetsize = m;

   Free(iarray); Free(x_sorted); Free(w_sorted); Free(w_sortedcpy); 
   Free(center); Free(dist_sorted); Free(scatter); Free(work);
}

/*****************************************************************************\
|*  weighted BACON							     *|
|*    x		  array[n, p]                                                *|
|*    w		  array[n]                                                   *|
|*    center	  on return: array[p]                                        *|
|*    scatter	  on return: array[p, p]                                     *|
|*    dist	  on return: array[n]                                        *|
|*    n, p	  dimensions                                                 *|
|*    alpha	  scalar                                                     *|
|*    subset	  array[n]                                                   *|
|*    cutoff	  on return: scalar                                          *|
|*    maxiter	  maximal number of iterations                               *|
|*    verbose	  0 = quiet; 1 =  verbose				     *|
\*****************************************************************************/
void wbacon(double *x, double *w, double *center, double *scatter, double *dist, 
   int *n, int *p, double *alpha, int *subset, double *cutoff, 
   int *maxiter, int *verbose)
{
   int *subset0;
   double *w_cpy, *work; 

   subset0 = (int*) Calloc(*n, int); 
   w_cpy = (double*) Calloc(*n, double); 
   work = (double*) Calloc(*n * *p, double); 

   // STEP 0: establish initial subset 
   double dhalf = 0.5;
   for (int j = 0; j < *p; j++)	 // center: coordinate-wise median
      wquantile(x + *n * j, w, n, &dhalf, &center[j]);   

   weightedscatter(x, work, w, center, scatter, n, p);	 // covariance 
   mahalanobis(x, work, center, scatter, dist, n, p);	 // Mahalan. dist. 

   Memcpy(w_cpy, w, *n);
   int subsetsize; 
   initialsubset(x, w_cpy, dist, subset, &subsetsize, n, p, verbose);

   // STEP 1: update iteratively
   int iter = 1;
   for (;;) {

      if (*verbose) {
	 double percentage = 100.0 * (double)subsetsize / (double)*n;
	 Rprintf("Subset %d: n = %d (%.1f%%)\n", iter, subsetsize, percentage);
      }

      weightedmean(x, w_cpy, center, n, p);		       // mean on subset
      weightedscatter(x, work, w_cpy, center, scatter, n, p);  // cov on subset
      mahalanobis(x, work, center, scatter, dist, n, p);       // Mahalan. dist.

      // check for convergence (XOR current with last subset)
      int converged = 0;
      for (int i = 0; i < *n; i++) {
	 converged += subset0[i] ^ subset[i];    
	 if (converged) break;
      }
      if (converged == 0) { 
	 *maxiter = iter; 
	 break;
      } 

      *cutoff = cutoffval(*alpha, *n, subsetsize, *p);   

      // generate new subset
      Memcpy(subset0, subset, *n); 
      Memcpy(w_cpy, w, *n); 

      subsetsize = 0;
      for (int i = 0; i < *n; i++) {
	 if (dist[i] < *cutoff) {
	    subset[i] = 1;  
	    subsetsize += 1;
	 } else
	    w_cpy[i] = 0.0;
      }

      iter++;
      if (iter > *maxiter) break;
   }
 
   Free(subset0); Free(w_cpy); Free(work);
}

/*****************************************************************************\
|*  weighted mean (vector values)					     *|
|*    x		  array[n, p]				                     *|
|*    w		  array[n], weights			                     *|
|*    center	  on return: array[p]		      	                     *|
|*    n, p	  dimensions						     *|
\*****************************************************************************/
static inline void weightedmean(double *x, double *w, double *center, int *n, 
   int *p)
{
   for (int i = 0; i < *p; i++)
      center[i] = 0.0;

   double sum_w = 0.0;
   for (int i = 0; i < *n; i++) {
      sum_w += w[i]; 
      for (int j = 0; j < *p; j++) {
	 center[j] += x[*n * j + i] * w[i]; 
      }
   }

   for (int j = 0; j < *p; j++)
      center[j] /= sum_w; 
}

/*****************************************************************************\
|*  weighted covariance/ scatter matrix					     *|
|*    x		  array[n, p]				                     *|
|*    work	  array[n, p]						     *|
|*    w		  array[n], weights			                     *|
|*    center	  array[p]		  	      	                     *|
|*    scatter	  on return: array[p, p]	      	                     *|
|*    n, p	  dimensions						     *|
\*****************************************************************************/
static inline void weightedscatter(double *x, double *work, double *w, 
   double *center, double *scatter, int *n, int *p)
{
   Memcpy(work, x, *n * *p); 

   // center the data and multiply by sqrt(w) 
   double sum_w = 0.0;
   for (int i = 0; i < *n; i++) {
      sum_w += w[i];
 
      for (int j = 0; j < *p; j++) {
	 work[*n * j + i] -= center[j];
	 work[*n * j + i] *= sqrt(w[i]);
      }
   }

   // cross product: scatter matrix = t(work) %*% work; 
   const double done = 1.0, dzero = 0.0;
   F77_CALL(dgemm)("T",	"N", p, p, n, &done, work, n, work, n, &dzero, scatter, 
      p);	

   for (int i = 0; i < (*p * *p); i++)
      scatter[i] /= (sum_w - 1); 
}

/*****************************************************************************\
|*  Mahalanobis distance						     *|
|*    x		  array[n, p]				                     *|
|*    work	  array[n, p]						     *|
|*    center	  array[p]	       		      	                     *|
|*    scatter	  array[p, p]		     	      	                     *|
|*    dist	  on return: array[n]	     	      	                     *|
|*    n, p	  dimensions						     *|
\*****************************************************************************/
static inline void mahalanobis(double *x, double *work, double *center, 
   double *scatter, double *dist, int *n, int *p)
{
   Memcpy(work, x, *n * *p);		  // copy of 'x' 

   for (int i = 0; i < *n; i++)       
      for (int j = 0; j < *p; j++) 
	 work[*n * j + i] -= center[j];	  // center the data

   Memcpy(dist, scatter, *p * *p); // store 'scatter' temporarily in 'dist'

   // Cholesky decomposition of scatter matrix (stored in 'dist') 
   int info;
   F77_CALL(dpotrf)("L", p, dist, p, &info);
   if (info != 0) error("mahalanobis: DPOTRF failed\n");

   // Solve for y in A * y = B by forward subsitution (A = Cholesky factor) 
   double done = 1.0;
   F77_CALL(dtrsm)("R",	"L", "T", "N", n, p, &done, dist, p, work, n);		
 
   for (int i = 0; i < *n; i++) {   // Mahalanobis distance
      dist[i] = 0.0;

      for (int j = 0; j < *p; j++) 
	 dist[i] += _POWER2(work[*n * j + i]);

      dist[i] = sqrt(dist[i]);
   }
}

#undef _POWER2 

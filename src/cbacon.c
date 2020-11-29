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

#include "cbacon.h" 

// macros
#define R_PACKAGE 1		   // 1: *.dll/*.so for R; 0: standalone binary 
#define DEBUG_MODE 0		   // 1: debug mode on; 0: off
#define _RANK_TOLERANCE 1.0e-8	   // criterion to detect rank deficiency
#define _POWER2(_x) ((_x) * (_x))

#if R_PACKAGE			    
   #define PRINT_TO_CONSOLE(_f, ...) Rprintf((_f), ##__VA_ARGS__)
#else
   #define PRINT_TO_CONSOLE(_f, ...) printf((_f), ##__VA_ARGS__)
#endif

// prototype of local function
static inline double cutoffval(double, int, int, int) 
   __attribute__((always_inline));
static inline void initialsubset(double*, double*, double*, int*, int*, int*, 
   int*, int*)  __attribute__((always_inline));
static inline void weightedscatter(double*, double*, double*, double*, double*, 
   int*, int*) __attribute__((always_inline));
static inline void weightedmean(double*, double*, double*, int*, int*) 
   __attribute__((always_inline));
static inline void mahalanobis(double*, double*, double*, double*, double*, 
   double*, int*, int*) __attribute__((always_inline));
static inline double qchisq2(double, double) __attribute__((always_inline));

#if DEBUG_MODE
void print_cov(double*, int*);
#endif

/******************************************************************************\
|* Quantile of chi-square distr. (approximation by Lin, 1994)		      *|	
|*    p	    probability							      *|	
|*    n	    degrees of freedom						      *|	
|*									      *|	
|* Lin, J.-T. (1994). New approximations for the precentage points of the     *|	
|*    chi-square distribution, Probab. Eng. Inf. Sci. 8, pp. 135-146.	      *|	
\******************************************************************************/
static inline double qchisq2(double p, double n) 
{
   return _POWER2(-0.655 + 0.975 * sqrt(n + 0.6) + 1.839 * sqrt(-log10(p))) - 3.0;
}

/*****************************************************************************\
|* Chi-squared cutoff value (p-dimensional, sample size n	  	     *|
|*    alpha	  scalar (0 <= alpha <= 1)		                     *|
|*    k		  size of subset					     *|
|*    n, p	  dimensions	 					     *|
\*****************************************************************************/
static inline double cutoffval(double alpha, int n, int k, int p)
{
   double dn = (double)n, dp = (double)p, dk = (double)k;
   double h = (dn + dp + 1.0) / 2.0;
   double chr = fmax(0.0, (h - dk) / (h + dk));
   double cnp = 1.0 + (dp + 1.0) / (dn - dp) + 2.0 / (dn - 1.0 - 3.0 * dp);
   return (cnp + chr) * sqrt(qchisq2(alpha / dn, dp));
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
   double *w_cpy, *work_np, *work_pp, *work_2n; 

   subset0 = (int*) Calloc(*n, int);	     
   w_cpy = (double*) Calloc(*n, double); 
   work_np = (double*) Calloc(*n * *p, double); 
   work_pp = (double*) Calloc(*p * *p, double); 
   work_2n = (double*) Calloc(2 * *n, double); 


   // STEP 0: establish initial subset 
   double dhalf = 0.5;
   for (int j = 0; j < *p; j++)	 // center: coordinate-wise median
      wquantile_noalloc(x + *n * j, w, work_2n, n, &dhalf, &center[j]);   

   // scatter and Mahalanobis distances
   weightedscatter(x, work_np, w, center, scatter, n, p);   
   mahalanobis(x, work_np, work_pp, center, scatter, dist, n, p);  

   Memcpy(w_cpy, w, *n);  // copy w (w_cpy will be set 0 if not in subset)
   int subsetsize; 
   initialsubset(x, w_cpy, dist, subset, &subsetsize, n, p, verbose);

   // STEP 1: update iteratively
   int iter = 1, is_different;
   for (;;) {

      if (*verbose) {
	 double percentage = 100.0 * (double)subsetsize / (double)*n;
	 if (iter > 1)
	    PRINT_TO_CONSOLE("Subset %d: n = %d (%.1f%%); cutoff: %.2f\n", iter, 
	       subsetsize, percentage, *cutoff);
	 else 
	    PRINT_TO_CONSOLE("Subset %d: n = %d (%.1f%%)\n", iter, subsetsize, 
	       percentage);
      }

      // location, scatter and the Mahalanobis distances
      weightedmean(x, w_cpy, center, n, p);			  
      weightedscatter(x, work_np, w_cpy, center, scatter, n, p);  
#if DEBUG_MODE
      print_cov(work_np, p);
#endif
      mahalanobis(x, work_np, work_pp, center, scatter, dist, n, p);       

      // check whether the subsets differ (XOR current with previous subset)
      is_different = 0;
      for (int i = 0; i < *n; i++) {
	 if (subset0[i] ^ subset[i]) {   
	    is_different = 1;
	    break;
	 }
      }
      if (is_different == 0) { 
	 *maxiter = iter; 
	 break;
      } 

      // chi-square cutoff value (quantile)
      *cutoff = cutoffval(*alpha, *n, subsetsize, *p);   
#if DEBUG_MODE
      PRINT_TO_CONSOLE("cutoff: %.4f\n", *cutoff);
#endif

      // generate new subset (based on updated Mahalanobis dist.)
      Memcpy(subset0, subset, *n); 
      Memcpy(w_cpy, w, *n); 

      subsetsize = 0;
      for (int i = 0; i < *n; i++) {
	 if (dist[i] < *cutoff) {      
	    subset[i] = 1;	       
	    subsetsize += 1;
	 } else {		       
	    subset[i] = 0;
	    w_cpy[i] = 0.0;	    // weight = 0 (if obs. is not in subset)
	 }
      }

      iter++;
      if (iter > *maxiter) break;
   }

   Free(subset0); Free(w_cpy); Free(work_np); Free(work_pp); Free(work_2n); 
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

   if (sum_w > DOUBLE_EPS) 
      for (int j = 0; j < *p; j++)
	 center[j] /= sum_w; 
   else 
      error("Weighted mean: division by zero\n");
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

   if (sum_w > DOUBLE_EPS) 
      for (int i = 0; i < (*p * *p); i++)
	 scatter[i] /= (sum_w - 1); 
   else 
      error("Weighted scatter: division by zero\n");
}

/*****************************************************************************\
|*  Mahalanobis distance						     *|
|*    x		  array[n, p]				                     *|
|*    work_np	  array[n, p]						     *|
|*    work_pp	  array[p, p]						     *|
|*    center	  array[p]	       		      	                     *|
|*    scatter	  array[p, p]		     	      	                     *|
|*    dist	  on return: array[n]	     	      	                     *|
|*    n, p	  dimensions						     *|
\*****************************************************************************/
inline void mahalanobis(double *x, double *work_np, double *work_pp, 
   double *center, double *scatter, double *dist, int *n, int *p)
{
   Memcpy(work_np, x, *n * *p);		     // copy of 'x' 

   for (int i = 0; i < *n; i++)       
      for (int j = 0; j < *p; j++) 
	 work_np[*n * j + i] -= center[j];   // center the data

   // Cholesky decomposition of scatter matrix 
   Memcpy(work_pp, scatter, *p * *p); 
   int info;
   F77_CALL(dpotrf)("L", p, work_pp, p, &info);
   if (info != 0) error("Mahalanobis distance: DPOTRF failed\n");

   // Solve for y in A * y = B by forward subsitution (A = Cholesky factor) 
   double done = 1.0;
   F77_CALL(dtrsm)("R",	"L", "T", "N", n, p, &done, work_pp, p, work_np, n);		

   for (int i = 0; i < *n; i++) {   // Mahalanobis distance
      dist[i] = 0.0;

      for (int j = 0; j < *p; j++) 
	 dist[i] += _POWER2(work_np[*n * j + i]);

      dist[i] = sqrt(dist[i]);
   }
}

/*****************************************************************************\
|* Print square matrix							      *|
\*****************************************************************************/
#if DEBUG_MODE
void print_cov(double *x, int *p)
{
   for (int i = 0; i < *p; i++) {
      for (int j = 0; j < *p; j++) {
	 PRINT_TO_CONSOLE("%.4f\t", x[j * *p + i]);
      }
      PRINT_TO_CONSOLE("\n");
   }
   PRINT_TO_CONSOLE("\n");
}
#endif
#undef _POWER2 

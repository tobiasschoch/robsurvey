/* weighted quantile and selection of k-th largest element 

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

   Reference:  Bentley, J.L. and D.M. McIlroy (1993). Engineering a 
	       Sort Function, Software - Practice and Experience 23, 
	       pp. 1249-1265    
   Note:       The extension of the method to weighted problems is ours.	
*/

# define _medium_array	40 // switch from insertion sort to quickselect; 
			   // using https://github.com/google/benchmark on
			   // Intel i7-8550U (mobile processor, 2017)

# define _large_array	50 // pivotal element determined by ninther
# define DEBUG_MODE	0  // debug mode (0 = off; 1 = activated)

#include "wquantile.h"

static inline void swap2(double*, double*, int, int) 
   __attribute__((always_inline));
static inline double med3(double*, int, int, int) 
   __attribute__((always_inline));
static inline int min(int, int) __attribute__((always_inline));
static inline int choose_pivot(double*, int, int)
   __attribute__((always_inline));
static inline int is_equal(double, double) __attribute__((always_inline));
static inline void partition_3way(double*, double*, int, int, int*, int*)
   __attribute__((always_inline));
static inline double insertionselect(double*, double*, int, int, double) 
   __attribute__((always_inline));

void wquant0(double*, double*, double, int, int, double, double*);

// debugging tools
#if DEBUG_MODE 
#include <stdio.h>
#include <string.h>	   
void debug_print_data(double*, double*, int, int, char*);
void debug_print_state(int, int);
#endif

/******************************************************************************\
|* wquantile  : weighted quantile					      *|	
|*									      *| 
|* array       array[lo..hi]						      *|	
|* weights     array[lo..hi]						      *|
|* n	       dimension						      *|	
|* prob	       probability defining quantile (0 <= prob <= 1)		      *|	
|* result      on return: weighted quantile				      *|	
\******************************************************************************/
void wquantile(double *array, double *weights, int *n, double *prob, 
   double *result)
{
   if (is_equal(*prob, 0.0)) {			// prob = 0.0
      wselect0(array, weights, 0, *n - 1, 0);
      *result = array[0];
   } else if (is_equal(*prob, 1.0)) {		// prob = 1.0
      wselect0(array, weights, 0, *n - 1, *n - 1);
      *result = array[*n - 1];
   } else {				
      // make copy 'array' and 'weights' because 'wquant0' is destructive
      double *array2;
      // array2 = [array[0..(n-1)], weights[0..(n-1)]], i.e. 'weights' is 
      // appended to 'array' s.t. we have one contiguous chunk of memory 
      // => cache optimized
      array2 = (double*) Calloc(2 * *n, double);
      Memcpy(array2, array, *n); 
      Memcpy(array2 + *n, weights, *n); 
      wquant0(array2, array2 + *n, 0.0, 0, *n - 1, *prob, result);      
      Free(array2); 
   }
}

/******************************************************************************\
|* wquant0  : weighted quantile (recursive function; for internal use)	      *|	
|*									      *| 
|* array    array[lo..hi]						      *|	
|* weights  array[lo..hi]						      *|
|* sum_w    total sum of weights: initialized with 0.0			      *|
|* lo	    dimension (usually: 0)					      *|	
|* hi	    dimension (usually: n - 1)					      *|	
|* prob	    probability defining quantile (0 <= prob <= 1)		      *|	
|* result   on return: weighted median					      *|	
|*									      *| 
|* NOTE:    wquant0 uses a weighted quickselect based on the the 3-way	      *|
|*	    partitioning of Bentley & McIlroy's (1993)			      *| 
\******************************************************************************/
void wquant0(double *array, double *weights, double sum_w, int lo, int hi, 
   double prob, double *result)
{ 
   #if DEBUG_MODE 
   debug_print_data(array, weights, lo, hi, "init");
   #endif

   if (hi <= lo) return;   // case: n = 1

   if (hi - lo == 1) {	   // case: n = 2
      double one_minus = 1.0 - prob;
      if (is_equal(one_minus * weights[lo], prob * weights[hi])) 
	 *result = (array[lo] + array[hi]) / 2.0;  	 
      else if (one_minus * weights[lo] > prob * weights[hi]) 
	 *result = array[lo];
      else 
	 *result = array[hi];

      return; 
   }

   if (sum_w < DBL_EPSILON) { // sum_w is only computed at initialization 
      for (int k = lo; k <= hi; ++k) 
	 sum_w += weights[k]; 
   }

   // case: n <= _medium_array
   if (hi - lo + 1 <= _medium_array) {
      *result = insertionselect(array, weights, lo, hi, prob);
      return;
   }

   // case: n > _medium_array: weighted quickselect
   // Bentley-McIlroy's 3-way partitioning (weighted): the positions of the 
   // sentinels 'i' and 'j' are returned 
   int i, j; 
   partition_3way(array, weights, lo, hi, &i, &j);

   // sum of weights of the elements smaller and larger than the pivot 
   // (determined with the help of the sentinels' positions) 
   double sum_w_lo = 0.0, sum_w_hi = 0.0; 
   for (int k = lo; k <= j; ++k) 
      sum_w_lo += weights[k]; 
   for (int k = i; k <= hi; ++k) 
      sum_w_hi += weights[k];

   #if DEBUG_MODE 
   debug_print_data(array, weights, lo, hi, "");
   debug_print_state(i, j);
   #endif

   // termination criterion: sum of weights on both sides are smaller than 0.5
   if (sum_w_lo < prob * sum_w && sum_w_hi < (1.0 - prob) * sum_w) {
      *result = array[j + 1]; 
   } else {  
      // tail recursion only on the partitioning with larger sum of weights
      if ((1 - prob) * sum_w_lo > prob * sum_w_hi) { 
	 weights[j + 1] = sum_w - sum_w_lo;  // weight of ignored part is dumped 
	 wquant0(array, weights, sum_w, lo, j + 1, prob, result);
      }
      else { 
	 weights[i - 1] = sum_w - sum_w_hi;  // weight is dumped
	 wquant0(array, weights, sum_w, i - 1, hi, prob, result);
      }
   }
}

/******************************************************************************\
|* wselect0 select the k-th largest elemens (with weights; for internal use)  *|	
|*									      *| 
|* array    array[lo..hi]						      *|	
|* weights  array[lo..hi]						      *|
|* lo	    dimension (usually: 0)					      *|	
|* hi	    dimension (usually: n - 1)					      *|	
|* k	    integer in 0:(n - 1)					      *|	
|*									      *| 
|* NOTE	    wselect0 sorts 'array' partially, such that element 'k' is in its  *| 
|*	    final (sorted) position => array[k] gives the k-th largest element*| 
\******************************************************************************/
void wselect0(double *array, double *weights, int lo, int hi, int k)
{
   if (hi <= lo) return;   // case: n = 1

   // Bentley-McIlroy's 3-way partitioning (weighted): the positions of the 
   // sentinels 'i' and 'j' are returned 
   int i, j; 
   partition_3way(array, weights, lo, hi, &i, &j);

   // tail recursion only on the partitioning where element 'k' lies  
   if (k <= j) { 
      wselect0(array, weights, lo, j, k);
   }
   else if (k >= i) { 
      wselect0(array, weights, i, hi, k);
   }
}

/******************************************************************************\
|* partition_3way: Bentley and McIlroy's (1993) 3-way partitioning scheme     *| 
\******************************************************************************/
static inline void partition_3way(double *array, double *weights, int lo, 
   int hi, int *i, int *j)
{
   // determine pivot and swap it into position 'lo' (i.e., position 0)  
   swap2(array, weights, choose_pivot(array, lo, hi), lo);
   double pivot = array[lo];

   // Bentley-McIlroy's 3-way partitioning weighted with sentinels i and j, 
   // respectively, scanning up and down until they cross; elements equal to 
   // the pivot are swapped to the far left and right, 

   int p = lo, q = hi + 1;
   *i = lo; *j = hi + 1;		  // initialize the sentinels

   for (;;) {
      while (array[++(*i)] < pivot)
	 if (*i == hi) break;

      while (pivot < array[--(*j)])
	 if (*j == lo) break;

      if (*i == *j && is_equal(array[*i], pivot))
	 swap2(array, weights, ++p, *i);

      if (*i >= *j) break;		  // check if sentinels cross 
      swap2(array, weights, *i, *j);	       

      // swap equal elements to the far left and right, respectively
      if (is_equal(array[*i], pivot)) swap2(array, weights, ++p, *i);
      if (is_equal(array[*j], pivot)) swap2(array, weights, --q, *j);
   }

   // swap equal elements from the borders back to the center
   *i = *j + 1;
   for (int k = lo; k <= p; k++)
      swap2(array, weights, k, (*j)--);

   for (int k = hi; k >= q; k--)
      swap2(array, weights, k, (*i)++);
}

/******************************************************************************\
|* choose_pivot: choose pivotal element (methods depend on the array size     *| 
\******************************************************************************/
static inline int choose_pivot(double *array, int lo, int hi)
{
   int n = hi - lo + 1; 
   int mid = lo + n / 2;	 // small array: median of three
   if (n > _large_array) {	 // large array: Tukey's ninther
      int eps = n / 8;
      lo = med3(array, lo, lo + eps, lo + eps + eps);
      mid = med3(array, mid - eps, mid, mid + eps);
      hi = med3(array, hi - eps - eps, hi - eps, hi); 
   }
   return med3(array, lo, mid, hi);
}

/******************************************************************************\
|* swap: swap two elements in double array				      *| 
\******************************************************************************/
static inline void swap2(double *array, double *weights, int i, int j)
{  
   double tmp = array[i]; array[i] = array[j]; array[j] = tmp;
   // swap weights
   tmp = weights[i]; weights[i] = weights[j]; weights[j] = tmp;
}

/******************************************************************************\
|* min: minimum of two integer values					      *|
\******************************************************************************/
static inline int min(int a, int b)
{ 
   return a < b ? a : b;
}

/******************************************************************************\
|* is_equal: check whether two doubles are equal using Knuth's notion of      *|
|* essential equality							      *| 
\******************************************************************************/
static inline int is_equal(double a, double b)
{
   return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : 
      fabs(a)) * DBL_EPSILON);
}

/******************************************************************************\
|* med3: median-of-three (but without swaps)		 		      *| 
\******************************************************************************/
static inline double med3(double *array, int i, int j, int k)
{
   return array[i] < array[j] ? 
         (array[j] < array[k] ? j : array[i] < array[k] ? k : i) 
       : (array[j] > array[k] ? j : array[i] > array[k] ? k : i);
}

/******************************************************************************\
|* insertionselect (i.e., insertion sort => select with weights)	      *| 
\******************************************************************************/
static inline double insertionselect(double *array, double *weights, int lo, 
   int hi, double prob) 
{
   // part: sort
   int exch = 0;		       
   for (int i = hi; i > lo; i--) {	  // smallest element as sentinel
      if (array[i] < array[i - 1]) {
	 swap2(array, weights, i, i - 1);
	 exch++;
      }
   }
   if (exch != 0){			  // insertion sort with half-exchanges
      for (int i = lo + 2; i <= hi; i++) {
	 double pivot = array[i];
	 double pivot_weight = weights[i];
	 int j = i;
	 while (pivot < array[j - 1]) {
	    array[j] = array[j - 1];
	    weights[j] = weights[j - 1];
	    j--;
	 }
	 array[j] = pivot;
	 weights[j] = pivot_weight;
      }
   }

   // part: select 
   double sum_w = 0.0;
   for (int k = lo; k <= hi; k++)	  // total sum of weight
      sum_w += weights[k]; 

   int k;			        
   double cumsum = 0.0; 
   for (k = lo; k <= hi; k++) {		  // cumulative sum of weight
      cumsum += weights[k]; 
      if (cumsum > prob * sum_w) 
	 break; 
   }

   if (k == lo)
      return array[k];
   else {
      cumsum -= weights[k];
      if (is_equal((1 - prob) * cumsum, prob * (sum_w - cumsum)))
	 return (array[k - 1] + array[k]) / 2.0;  	 
      else
	 return array[k];
   }
}

/******************************************************************************\
|* DEBUGGING TOOLS							      *| 
\******************************************************************************/
#if DEBUG_MODE 
void debug_print_data(double *array, double *weights, int lo, int hi, 
   char *message)
{ 
   if (strlen(message) > 0) {
      printf("------\n");
      printf("%s x\t", message);
   } else 
      printf("x \t");

   for (int i = lo; i <= hi; ++i) 
      printf("%.2f\t", array[i]);

   printf("\n");
   //
   if (strlen(message) > 0) 
      printf("%s w\t", message);
   else 
      printf("w \t");

   for (int i = lo; i <= hi; ++i) 
      printf("%.2f\t", weights[i]);

   printf("\n");
}

void debug_print_state(int i, int j)
{ 
   printf("final:\ti = %d\tj = %d\n", i, j);
}
#endif 

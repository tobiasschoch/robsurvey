/*****************************************************************************\
|* wquantile: weighted quantiles 					     *|
|* ------------------------------------------------------------------------- *|
|* PROJECT  sctbase 							     *|
|* SUBEJCT  weighted quantile						     *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), February 10, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include "wquantile.h"

#define DEBUG_MODE 0	// debug mode (0 = off; 1 = activated)

// declaration of 'local' functions (inline imperative is GCC specific)
static inline int partition(double*, int, int, int) 
   __attribute__((always_inline));
static inline void swap(double*, int, int, int)
   __attribute__((always_inline));
static inline void median_of_three(double*, int, int, int)
   __attribute__((always_inline));

#if DEBUG_MODE 
// print function in debug mode
static inline void print(double*, int, int, int, int, int)
   __attribute__((always_inline));
static inline void print(double *xw, int n, int left, int right, int i, 
			 int init)
{
   if (init){
      for (int k = 0; k < n; k++){
	 Rprintf("%d\t", k);
      }
      Rprintf("\n");
   } else{
      for (int k = 0; k < n; k++){
	 Rprintf("%.2f\t", xw[k]);
      }
      Rprintf("left = %d, right = %d, i = %d\n", left, right, i);
   }
}
#endif

/*****************************************************************************\
|* wquantile: weighted quantiles					     *|
|*									     *| 
|* PARAMETERS						    		     *|
|*    x	       data, array[n]						     *|
|*    w	       weights, array[n]					     *|
|*    ptrn     dimension n						     *|
|*    ptrq     prob q in [0,1]						     *|
|*    ptresult on return, q-th quantile					     *|
\*****************************************************************************/
void wquantile(double *x, double *w, int *ptrn, double *ptrq, 
	       double *ptrresult)
{ 
   int i, left = 0, right = *ptrn - 1;
   double sum_w = 0.0, sum_w_lo = 0.0, sum_w_hi = 0.0;
   double *xw;
  
   // weight total 
   for (int k = 0; k < *ptrn; k++) sum_w += w[k];

   // generate array xw = (x, w)  
   xw = (double*) Calloc(2 * *ptrn, double); 
   Memcpy(xw, x, *ptrn); 
   Memcpy(xw + *ptrn, w, *ptrn); 

   // normalized weight, s.t. sum(w_i) = 1
   for (int k = 0; k < *ptrn; k++) xw[k + *ptrn] /= sum_w;

#if DEBUG_MODE 
Rprintf("initial\n");
print(xw, *ptrn, left, right, 0, 1);
print(xw, *ptrn, left, right, 0, 0);
Rprintf("---\n");
#endif

   while (right > left){
      
      // choosing the pivotal element by the median-of-three method
      median_of_three(xw, *ptrn, left, right); 

#if DEBUG_MODE 
Rprintf("median swaps\n");
print(xw, *ptrn, left, right, i, 0);
#endif

      i = partition(xw, *ptrn, left, right);

#if DEBUG_MODE 
Rprintf("partitioning\n");
print(xw, *ptrn, left, right, i, 0);
#endif

      // here, we decide which partition to keep active (i.e., partition that
      // contains the k-th element); first compute the weights below and above 
      // of the pivotal element
      sum_w_lo = 0.0;
      for (int k = 0; k < i; k++) sum_w_lo += xw[k + *ptrn];
      sum_w_hi = 0.0;
      for (int k = i + 1; k < *ptrn; k++) sum_w_hi += xw[k + *ptrn];

//FIXME: compute weights lo and hi only for range [left, right]; in the 
//	 decision (below), we put the weights of the inactive partition onto 
//	 the pivot

      if (sum_w_lo >= *ptrq){
	 right = i - 1;
      } else if (sum_w_hi <= 1.0 - *ptrq){ 
	 break;
      } else{
	 left = i + 1;
      }
   }

#if DEBUG_MODE 
Rprintf("on return\n");
print(xw, *ptrn, left, right, i, 0);
#endif

   // modification
   if (left == right) i = left;

   *ptrresult = xw[i];
   Free(xw); 
}

// old implementation of the decision (only median)
      /* sum_w_lo = 0.0; */
      /* for (int k = 0; k < i; k++) sum_w_lo += xw[k + *ptrn]; */
      /* sum_w_hi = 0.0; */
      /* for (int k = i + 1; k < *ptrn; k++) sum_w_hi += xw[k + *ptrn]; */
      /*  */
      /* if (sum_w_lo >= sum_w_hi + xw[i + *ptrn]){ */
	/*  right = i - 1; */
      /* } else if (sum_w_lo + xw[i + *ptrn] >= sum_w_hi){  */
	/*  break; */
      /* } else{ */
	/*  left = i + 1; */
      /* } */



/*****************************************************************************\
|* A.C.R. Hoare's partitioning scheme (part of algorithm 'FIND'; i.e. quick- *|
|* select); see also Knuth (1998): The Art of Computer Programming, vol. 3,  *|
|* 2nd ed., Addison-Wesley, p. 113-123 and Exercise 31 on p. 131	     *|
\*****************************************************************************/
static inline int partition(double *a, int n, int left, int right)
{ 
   int i = left - 1;
   int j = right; 
   double pivot = a[right]; // the pivot's position = rightmost element of a 

   // Hoare's partitioning scheme (i.e. two sentinels or pointers, i and j, 
   // that scan up (i) and scan down (j) until they cross (i.e., i >= j)
   for (;;) { 
      // scan up to find values larger than pivot (note that we start at 
      // i = left -1; also, the scan is automatically stopped because the 
      // rightmost element in a is the pivot) 
      while (a[++i] < pivot);
      // scan up to find values smaller than pivot (the scan must be stopped
      // at the leftmost element )
      while (pivot < a[--j]) if (j == left) break;
      // when the pointers cross, we terminate
      if (i >= j) break;
      // otherwise we swap the elements  (i.e., we found two elements s.t. 
      // x[i] > pivot and x[j] < pivot; hence, we swap them)
      swap(a, n, i, j);
   }
   // insert pivotal element
   swap(a, n, i, right);
   return i;
}

/*****************************************************************************\
|* swap two elements of array a = (x, w), which is of size 2 * n; elements   *|
|* in x and w are swapped in 'parallel'					     *|
\*****************************************************************************/
static inline void swap(double *a, int n, int i, int j)
{  
   double tmp;
   // swap data values (first n observations in array a)
   tmp = a[i]; a[i] = a[j]; a[j] = tmp;
   // swap weights (second n observations in array a)
   tmp = a[i + n]; a[i + n] = a[j + n]; a[j + n] = tmp;
}

/*****************************************************************************\
|* median_of_three: we sort x s.t. x[left] < x[mid] < x[right] and then we   *|
|* take x[mid] as pivotal or partitioning element rather than x[right]       *|
|* See Sedgewick (1983): Algorithms, Addison-Wesley, p. 112-113	  	     *|
\*****************************************************************************/
static inline void median_of_three(double *x, int n, int left, int right)
{
   int mid = (left + right) / 2;
   if (x[right] < x[left]){
      swap(x, n, left, right);        
   }
   if (x[mid] < x[left]){ 
      swap(x, n, mid, left);
   }
   if (x[right] < x[mid]){
      swap(x, n, right, mid);
   }
   // move 'mid' to 'right' (right = pivot's position)
   swap(x, n, mid, right);
}



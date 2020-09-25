/*****************************************************************************\
|* PROJECT  robsurvey							     *|
|* SUBEJCT  Normalization constant of weighted MAD for Mallows GM-estimator  *|
|* AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), September 25, 2020	     *|
|* LICENSE  GPL >= 2							     *|
|* COMMENT  [none]							     *|
\*****************************************************************************/
#include <R.h>
#include <Rmath.h>
#include "zeroin.h" 

// parameters of R_zeroin2 (see 'zeroin.c') 
#define zeroin_MAXIT 30
#define zeroin_TOL 1e-5
#define zeroin_LEFTBOUND 0.1
#define zeroin_RIGHTBOUND 10

struct zeroin_data_t {
   int n;
   double totalweight;
   double data[];
};

double zeroin_mallows_mad(double, void*);

/*****************************************************************************\
|* Compute normalization constant (wrapper function for Richard Brent's root  *|
|* finding algorithm 'zeroin')						      *|
|*									      *|
|*    xwgt  weights in the design space, array[n]			      *|
|*    n	    dimension							      *|
\*****************************************************************************/
double mallows_mad_normalization(double *xwgt, int *n)
{
   // allocate structure
   struct zeroin_data_t *info;
   info = (struct zeroin_data_t*) malloc(sizeof(struct zeroin_data_t) + 
      sizeof(double) * *n);

   // intialize structure 
   info->n = *n;
   Memcpy(info->data, xwgt, *n);

   // compute the root with Richard Brent's zeroin algorithm ('zeroin.c')
   int Maxit = zeroin_MAXIT;
   double ax = zeroin_LEFTBOUND, bx = zeroin_RIGHTBOUND, Tol = zeroin_TOL; 
   double fa = zeroin_mallows_mad(ax, info); 
   double fb = zeroin_mallows_mad(bx, info); 
   double result = R_zeroin2(ax, bx, fa, fb, zeroin_mallows_mad, info, &Tol, 
      &Maxit);	
   free(info);

   if (Maxit > 0)
      return (result > DBL_EPSILON) ? 1.0 / result : 1.4826;
   else {
      warning("Normalization for weighted MAD: not converged\n");
      return 1.4826;	// failure of convergence: the 'standard' is returned 
   }  
}

/*****************************************************************************\
|* Estimating equation whose root is the normalization constant of the        *|
|* weighted MAD for Mallows GM-estimator of regression.			      *|
|*									      *|
|*    beta     argument (root)						      *|
|*    info     void pointer (called with an instance of struct zeroin_data_t) *|
\*****************************************************************************/
double zeroin_mallows_mad(double beta, void *info)
{
   int n = ((struct zeroin_data_t*)info)->n; // cast void pointer to struct 
   double tmp = 0.0;

   for (int i = 0; i < n; i++) {
      tmp += pnorm(beta, 0.0, sqrt(((struct zeroin_data_t*)info)->data[i]), 
	 1, 0);		// function 'pnorm' is defined in 'R/include/Rmath.h'
   }
   
   tmp /= (double)n;
   return tmp - 0.75; 
}

#undef zeroin_MAXIT 
#undef zeroin_TOL 
#undef zeroin_LEFTBOUND
#undef zeroin_RIGHTBOUND 

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

// consistency correction constant of the median absolute deviations about
// the median at the Gaussian core model (see regression.c)
#define mad_NORM_CONSTANT 1.482602
// consistency correction constant interquartile range at the Gaussian core
// model (see huber2.c and regresion.c)
#define iqr_NORM_CONSTANT 0.741301

// parameterization of Brent's root finding algorithm R_zeroin2 ('zeroin.c')
#define zeroin_MAXIT 30             // maximum number of iterations
#define zeroin_TOL 1e-5             // numerical tolerance (comparison)
#define zeroin_LEFTBOUND 0.1        // lower bound of the search interval
#define zeroin_RIGHTBOUND 10        // upper bound of the search interval

#endif

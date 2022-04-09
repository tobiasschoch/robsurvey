#ifndef _REGRESSION_DATA_H
#define _REGRESSION_DATA_H

// structure: regression data
typedef struct regdata_struct {
    int n;
    int p;
    double *x;
    double *y;
    double *w;
    double *xwgt;
} regdata;

// structure: work arrays
typedef struct workarray_struct {
    int lwork;
    double *work_lapack;
    double *work_x;
    double *work_y;
    double *work_2n;
} workarray;

#endif

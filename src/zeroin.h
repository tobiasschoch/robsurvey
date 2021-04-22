#include <R.h>
#ifndef _ZEROIIN_H
#define _ZEROIIN_H
double R_zeroin2(double, double, double, double, double (*f)(double, void*),
    void*, double*, int*);
#endif

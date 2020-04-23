// RegisteringDynamic Symbols

#include <R_ext/Rdynload.h>
#include "robsurvey.h"

// create arrays describing each C routine 
static const R_CMethodDef cMethods[]  = {
   {"rwlslm", (DL_FUNC) &rwlslm, 14},   
   {"wmad", (DL_FUNC) &wmad, 4},   
   {"wquantile", (DL_FUNC) &wquantile, 5},   
   {"wmeantrimmed", (DL_FUNC) &wmeantrimmed, 6},   
   {"wmeanwinsorized", (DL_FUNC) &wmeanwinsorized, 6},   
   {NULL, NULL, 0}
};

// register the C routines to R
void R_init_rht(DllInfo* info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

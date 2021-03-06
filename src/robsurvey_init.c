// RegisteringDynamic Symbols
#include <R_ext/Rdynload.h>
#include "robsurvey.h"
#include "trimmedwinsorized.h"
#include "huberm.h"

// create arrays describing each C routine
static const R_CMethodDef cMethods[]  = {
    {"wtrimmedmean", (DL_FUNC) &wtrimmedmean, 7},
    {"wwinsorizedmean", (DL_FUNC) &wwinsorizedmean, 6},
    {"huberm", (DL_FUNC) &huberm, 11},
    {"rwlslm", (DL_FUNC) &rwlslm, 17},
    {"wquantile", (DL_FUNC) &wquantile, 5},
    {"wkwinsorizedmean", (DL_FUNC) &wkwinsorizedmean, 6},
    {"cov_reg_model", (DL_FUNC) &cov_reg_model, 13},
    {"cov_reg_design", (DL_FUNC) &cov_reg_design, 12},
    {NULL, NULL, 0}
};

// register the C routines to R
void R_init_robsurvey(DllInfo* info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

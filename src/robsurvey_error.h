#include <stddef.h>

#ifndef _ROBSURVEY_ERROR_H
#define _ROBSURVEY_ERROR_H

// error handling
typedef enum robsurvey_error_enum {
    ROBSURVEY_ERROR_OK = 0,                 // no error
    ROBSURVEY_ERROR_SCALE_ZERO,             // scale estimate is zero
    ROBSURVEY_ERROR_RANK_DEFICIENT,         // design matrix is rank deficient
    ROBSURVEY_ERROR_MALLOWS_NOT_CONVERGED,  // Mallows normalization const.
    ROBSURVEY_ERROR_COUNT,                  // [not an actual error type]
} robsurvey_error_type;

// declaration
const char* robsurvey_error(robsurvey_error_type);
#endif

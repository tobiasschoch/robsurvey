/* Error handling function and array with error messages

   Copyright (C) 2020-21 Tobias Schoch (e-mail: tobias.schoch@gmail.com)

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
*/

#include <stddef.h>
#include "robsurvey_error.h"

// human readable errors
const char* const ROBSURVEY_ERROR_STRINGS[] = {
    "no errors",
    "scale estimate is zero (or nearly so)",
    "design matrix is rank deficient (or nearly so)",
    "Mallows normalization constant: Algorithm did not converge",
    "QR factorization: dgeqrf failed",
    "QR factorization: dtrtri failed",
    "QR factorization: dorgqr failed"
};

// obtain a human readable error message
const char* robsurvey_error(robsurvey_error_type err)
{
    if (err >= ROBSURVEY_ERROR_COUNT)
        return NULL;
    else
        return ROBSURVEY_ERROR_STRINGS[err];
}

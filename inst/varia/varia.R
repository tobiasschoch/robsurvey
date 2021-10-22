# check whether the vector of (heteroscedastic) variances is in the column
# space of the design matrix
colspace <- function(x, var, tol = sqrt(.Machine$double.eps))
{
    R <- qr.R(qr(cbind(x, var)))
    any(abs(R[, NCOL(x) + 1]) < tol)
}



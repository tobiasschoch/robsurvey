#==============================================================================
# huberm: proposal 2 Huber M-estimator of location and scale
# -----------------------------------------------------------------------------
# PROJECT  sctbase
# SUBEJCT  R function and test cases
# AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), February 10, 2020
# LICENSE  GPL >= 2
# COMMENT  [none]
#==============================================================================
setwd("C:/My/C/sctbase")
dyn.load("sctbase.dll")


hub <- function(x, w = NULL, k, tol = 1e-6, maxit = 50){
    stopifnot(k > 0)
    n <- length(x)
    if (is.null(w))
        w <- rep(1, n)
    tmp <- .C("huberm", x = as.double(x), w = as.double(w),
        robwgt = as.double(numeric(n)), k = as.double(k),
        loc = as.double(numeric(1)), scale = as.double(numeric(1)),
        n = as.integer(n), maxit = as.integer(maxit),
        tol = as.double(tol))
    tmp$convergend <- TRUE
    if (tmp$maxit == maxit)
        tmp$convergend <- FALSE
    tmp
}

x <- c(200,0.55,-0.34,-0.91,0.61,0.45,-1.07,0.3,0.46,-1.22,-0.04,-1.19,-0.17,
    -0.31,0.99,0.35,-0.38,0.9,0.31,-1.35)

hub(x, k = 2) # 10

hub(x, k = 1.6) # 9

hub(x, k = 1.3) # 15

hub(x, k = 1.2) # 20

hub(x, k = 1.1) # 30



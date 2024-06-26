# compute first-order inclusion probabilities
pps_probabilities <- function(size, n)
{
    stopifnot(is.numeric(size), all(size >= 0), is.numeric(n))
    N <- length(size)
    if (n >= N)
        stop("Sample size must be smaller than population size\n")

    res <- .C(C_pps_prob, size = as.double(size), N = as.integer(N),
              n = as.integer(n), pik = as.double(numeric(N)))$pik
    attr(res, "n") <- n
    attr(res, "N") <- N
    class(res) <- "prob_pps"
    res
}

# S3 print method
print.prob_pps <- function(x, ...)
{
    n <- attr(x, "n"); N <- attr(x, "N")
    cat(paste0("Probabilities for pps sampling (w/o replacement; n = ", n,
        ", N = ", N,")\n"))
    n_wp1 <- sum(x == 1)
    if (n_wp1 == 1)
        cat("One element has prob. equal to 1\n")
    if (n_wp1 > 1)
        cat(n_wp1, "elements have prob. equal to 1\n")
}

# pps sampling
pps_draw <- function(x, method = "brewer", sort = TRUE)
{
    if (!inherits(x, "prob_pps"))
        stop("Argument 'x' must be an object generated by 'pps_probabilities()'\n")
    method <- match.arg(tolower(method), c("brewer"))
    n <- attr(x, "n"); N <- attr(x, "N")
    tmp <- .C(C_pps_brewer, probs = as.double(x), N = as.integer(N),
              n = as.integer(n), sample = as.integer(rep(0, n)))
    if (sort)
        sort(tmp$sample)
    else
        tmp$sample
}

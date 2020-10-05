weighted_mean <- function(x, w, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   return(sum(dat$x * dat$w) / sum(dat$w))
} 

weighted_total <- function(x, w, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   return(sum(dat$w * dat$x))
} 

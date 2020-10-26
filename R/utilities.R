summary.formula <- function(object, design, na.rm = FALSE, ...)
{
   mf <- stats::model.frame(object, design$variables, na.action = 
      stats::na.pass)
   n <- nrow(mf)
   xname <- names(mf)
   y <- mf[[1]]
   w <- as.numeric(1 / design$prob)

   if (is.factor(y)) {	# case: factor
      dat <- data.frame(y = y, w = w)
      cc <- stats::complete.cases(dat)
      if (sum(cc) != n) {
	 if (na.rm) {
	    dat <- dat[cc, ]
	 } else{
	    return(NULL)
	 }
      }
      res <- rbind(table(dat$y), sapply(split(dat$w, dat$y), sum)) 
      rownames(res) <- c("n", "N") 
   } else {		# case: numeric 
      tmp <- weighted_quantile(y, w, probs = c(0.25, 0.5, 0.75), na.rm)
      m <- weighted_mean(y, w, na.rm) 
      res <- c(min(y), tmp[1:2], m[[1]], tmp[3], max(y)) 
      names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
      res <- rbind(res, summary(y))
      rownames(res) <- c("weighted", "classical")
   }
   res
}

# weight associated with the Huber psi-function
huberWgt <- function(x, k = 1.345)
{
   pmin.int(1, k / abs(x))
} 

# weight associated with the Tukey biweight psi-function
tukeyWgt <- function(x, k = 4.685)
{
   (1 - (x / k)^2)^2 * (abs(x) <= k)
} 

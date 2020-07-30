#' Weighted five-number summary of a variable 
#'
#' Weighted five-number summary of a variable (similar to \code{base::summary} for \code{[numeric vectors]})
#'
#' A weighted five-number summary (numeric variable) or a frequency table (factor variable).  
#'
#' @param object \code{[character]} name of the variable for which a summary is desired. 
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param ... additional arguments. 
#' @return A numerical summary 
#' @method summary formula
#' @export
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

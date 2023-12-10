# summary method for survey.design object
svysummary <- function(object, design, na.rm = FALSE, ...)
{
    mf <- model.frame(object, design$variables, na.action = na.pass)
    n <- nrow(mf)
    y <- mf[[1]]
    w <- as.numeric(1 / design$prob)

    if (is.factor(y)) {
        dat <- data.frame(y = y, w = w)
        cc <- complete.cases(dat)
        if (sum(cc) != n) {
	        if (na.rm)
	            dat <- dat[cc, ]
	        else
	            return(NULL)
        }
        res <- rbind(table(dat$y), sapply(split(dat$w, dat$y), sum))
        rownames(res) <- c("n", "N")
    } else {
        tmp <- weighted_quantile(y, w, probs = c(0.25, 0.5, 0.75), na.rm)
        m <- weighted_mean(y, w, na.rm)
        res <- c(min(y), tmp[1:2], m[[1]], tmp[3], max(y))
        names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
        res <- rbind(res, summary(y))
        rownames(res) <- c("weighted", "classical")
    }
    res
}

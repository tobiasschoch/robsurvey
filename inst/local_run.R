#library(robsurvey); library(survey)
library(survey)
ROOT <- "C:/My/code/robsurvey"
load(paste0(ROOT, "/data/workplace.RData"))
dyn.load(paste0(ROOT, "/src/robsurvey.dll"))
f <- c("robline.R","robsurvey-internal.R","robsvyreg.R","svymean_huber.R",
   "svymean_reg.R","svymean_reg_huber.R", 
   "svymean_winsorized.R", "svyreg.R", "svyreg_huber.R","weighted_mad.R",
   "weighted_mean.R", "weighted_mean_huber.R", "weighted_mean_trimmed.R",
   "weighted_mean_winsorized.R", "weighted_quantile.R", "utilities.R",
   "weighted_mean_dalen.R", "svymean_trimmed.R", "svymean_k_winsorized.R")       
for (i in 1:length(f)){
   source(paste0(ROOT, "/R/", f[i]))
}

set.seed(1)
x <- rnorm(2, 100)
z <- 2 * x + rnorm(100)
w <- rep(1, length(x))
dn <- svydesign(ids = ~1, weights = ~w , data = data.frame(x = x, w = w, z = z))  


svymean_winsorized(~x, dn, 0, 1, simple_var = TRUE)



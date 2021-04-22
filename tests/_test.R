library(data.table)
library(survey)
ROOT <- "C:/My/code/robsurvey"
dyn.load(paste0(ROOT, "/src/robsurvey.dll"))
setwd(paste0(ROOT, "/R"))
files <- dir(pattern = ".R$")
files <- files[-which(files == "_test.R")]
sapply(files, source, .GlobalEnv)



# MU284 population
library(sampling)
data(MU284)
MU284 <- as.data.table(MU284)
# add region size
setkey(MU284, REG)
MU284 <- MU284[MU284[, .N, keyby = REG]]
setkey(MU284, LABEL)
# stratified srswor by region (sample size per stratum = 9)
dat <- as.data.table(sampling:::strata(MU284, "REG", rep(9, 8),
    method = "srswor"))
dat[, REG := NULL]
setkey(dat, ID_unit)
dat <- MU284[dat]
dat[, weight := 1 / Prob]
d_strat_MU284 <- svydesign(id = ~1, strata = ~ Stratum, weights = ~weight,
    fpc = ~N, data = dat)


w <- weights(d_strat_MU284)
x <- d_strat_MU284$variables$ME84

weighted_mean_huber(x, w)

svymean_huber(~ME84, d_strat_MU284)


robsvyreg(rep(1, length(x)), x, w, k = 100, 1, na.rm = TRUE)

svyreg_huber(ME84~P85, d_strat_MU284, k = 2)


weighted_mean_trimmed(x, w, LB = 0.05, UB = 0.95)

weighted_mean_winsorized(x, w, LB = 0.05, UB = 0.95)






library(sampling)
data(MU284)
#-------------------------------------------------------------------------------
# Stratified by region + take all stratum of large citites

# MU281: population of the 281 smallest municipalities (based on 1975 population)
MU281 <- MU284[MU284$P75 < 200, ]
MU281_N_strata <- sapply(split(MU281, MU281$REG), NROW)
N <- 281
n <- 57
MU281_n_strata <- unname(round(n / N * MU281_N_strata))
# stratified simple random sample from MU281
set.seed(1)
s <- strata(MU281, "REG", size = MU281_n_strata, method = "srswor")
MU281_strat <- getdata(MU281, s)
MU281_strat$ID_unit <- NULL
# three largest municipalities (Stockholm, Göteborg, and Malmö)
MU3 <- MU284[MU284$P75 > 200, ]
MU3$Prob <- 1
MU3$Stratum <- 9
# population stratum size
MU284_N_strata <- c(MU281_N_strata, "9" = 3)
# sample
MU284strat <- rbind(MU281_strat, MU3[, names(MU281_strat)])
MU284strat$weights <- 1 / MU284strat$Prob
MU284strat$Prob <- NULL
# types
double_var <- c("P85", "P75", "RMT85", "CS82", "SS82", "S82", "ME84", "REV84",
    "weights")
MU284strat[, double_var] <- apply(MU284strat[, double_var], 2, as.double)
int_var <- c("LABEL", "Stratum", "CL", "REG")
MU284strat[, int_var] <- apply(MU284strat[, int_var], 2, as.integer)
# add fpc
MU284strat$fpc <- NA_real_
for (i in 1:length(MU284_N_strata)) {
    MU284strat[
        MU284strat$Stratum == as.numeric(names(MU284_N_strata[i])), "fpc"] <-
            MU284_N_strata[i]
}
# save
save(MU284strat, file = "MU284strat.RData")

# design
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
    weights = ~weights, data = MU284strat)

#-------------------------------------------------------------------------------
# Sample with probabilites proportional to population size in 1975 (P75)
n <- 50
pik <- inclusionprobabilities(MU284$P75, n)
set.seed(1)
s_index <- UPbrewer(pik)
MU284pps <- MU284[s_index == 1, ]
MU284pps$pi <- pik[s_index == 1]
# calibrate the weights that are not equal to 1.0
MU284pps$intercept <- 1
MU284pps <- MU284pps[order(1 / MU284pps$pi), ]
at <- 4:50
totals <- c(284 - 3, sum(MU284$P75) - sum(MU284pps$P75[1:3]))
MU284pps$g[at] <- calib(MU284pps[at, c("intercept", "P75")],
    d = 1 / MU284pps$pi[at], total = totals, method = "truncated",
    bounds = c(0.8, 10))
MU284pps$weights[at] <- MU284pps$g[at] / MU284pps$pi[at]
MU284pps$weights[1:3] <- 1
MU284pps[, c("pi", "g", "intercept")] <- NULL
MU284pps$pi <- 1 / MU284pps$weights

# with(MU284pps, sum(weights))
# with(MU284pps, sum(weights * P75))
# plot(weights ~ P75, MU284pps)

save(MU284pps, file = "MU284pps.RData")

# design
dn <- svydesign(ids = ~LABEL, fpc = ~pi, data = MU284pps, pps = "brewer")
svytotal(~P75, dn)



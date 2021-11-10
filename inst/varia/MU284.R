library(sampling)
data(MU284)

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

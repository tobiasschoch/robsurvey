library(testthat)
library(MASS)
library(robsurvey, quietly = TRUE)
attach(workplace)

# make a copy of the function MASS::hubers
hubers_mod <- hubers
# replace the mad by the (scaled) IQR as initial scale estimator
body(hubers_mod)[[7]][[3]][[2]] <- substitute(s0 <- IQR(y, type = 2) * 0.7413)

TOLERANCE <- 1e-8
expect_equal(huber2(payroll, rep(1, length(payroll)), tol = TOLERANCE),
    hubers_mod(payroll, tol = TOLERANCE)$mu)

#-------------------------------------------------------------------------------
# VARIA
#-------------------------------------------------------------------------------
# my <- numeric(9)
# mass <- numeric(9)
# for (i in 1:9) {
#     my[i] <- huber2(payroll, rep(1, 142), tol = 10^-i)
#     mass[i] <- MASS::hubers(payroll, tol = 10^-i)$mu
# }
#
# plot(1:9, mass, ylim = range(mass, my), type = "n", ylab = "estimate",
#     xlab = expression(10^-x), main = "Estimate vs. Precision (tol)")
# points(1:9, my, type = "b", col = 2)
# points(1:9, mass, type = "b", col = 4)
# legend("topright", col = c(2, 4), lty = c(1, 1),
#     legend = c("robsurvey", "MASS"))
#
# reference function: for documentation purposes
# hh <- function(x, w, k = 1.5, tol = 1e-5, maxit = 50, df_cor = TRUE)
# {
#     loc0 <- weighted_median(x, w)
#     scale0 <- weighted_IQR(x, w)
#     w_total <- sum(w)
#     th <- 2 * pnorm(k) - 1
#     kappa <- th + k^2 * (1 - th) - 2 * k * dnorm(k)
#     for (i in 1:maxit) {
#         x_wins <- pmin(pmax(loc0 - k * scale0, x), loc0 + k * scale0)
#         loc1 <- sum(w * x_wins) / w_total
#         ssq <- sum(w * (x_wins - loc1)^2)
#         ssq <- if (df_cor) ssq / (w_total - 1) else ssq / w_total
#         scale1 <- sqrt(ssq / kappa)
#         if (abs(loc0 - loc1) < tol * scale0 && abs(scale1 / scale0 - 1) < tol) {
#             break
#         } else {
#             loc0 <- loc1
#             scale0 <- scale1
#         }
#     }
#     loc1
# }

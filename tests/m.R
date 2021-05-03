library(robsurvey, quietly = TRUE)
attach(workplace)
weighted_mean_huber(payroll, weight)

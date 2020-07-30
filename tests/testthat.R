library(testthat)

formula <- x
design <- dn

expect_equal(coef(svymean(~formula, design)), 
             coef(svymean_huber(~formula, design, k = Inf)))
expect_equal(SE(svytotal(~formula, design)), 
             SE(svytotal(~formula, design, k = Inf)))



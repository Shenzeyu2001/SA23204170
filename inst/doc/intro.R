## -----------------------------------------------------------------------------
library(Rcpp)
library(SA23204170)
error <- generateARMA(3,4,1000,c(.5,-.2),c(.3))
l_hat <-BACM(error,20,30)
l_hat


rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/mpsBART.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
n_ <- 50
set.seed(42)

# Simulation 1
fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = 0.01)
x <- fried_sim$x[,1:5,drop = FALSE]
x_new <- x
y <- fried_sim$y

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)


# Testing the mpsBART
bart_test <- rbart(x_train = x,y = unlist(c(y)),x_test = x_test,
                   n_tree = 1,n_mcmc = 2500,alpha = 0.95,dif_order = 0,
                   beta = 2,nIknots = 100,delta = 1,
                   n_burn = 500,scale_bool = TRUE)


# Running BART
bartmod <- dbarts::bart(x.train = x,y.train = unlist(c(y)),ntree = 200,x.test = x_test)

# Convergence plots
par(mfrow = c(1,2))
plot(bart_test$tau_post,type = "l", main = expression(tau),ylab=  "")
plot(bartmod$sigma^-2, type = "l", main = paste0("BART: ",expression(tau)),ylab=  "")

par(mfrow = c(1,2))
plot(bart_test$y_hat %>% rowMeans(),y, main = 'mpsBART', xlab = "mpsBART pred", ylab = "y")
plot(bartmod$yhat.train.mean,y, main = "BART", xlab = "BART pred", ylab = "y")

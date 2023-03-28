rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/mpsBART.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
n_ <- 100
set.seed(42)

# Simulation 1
fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = 1)
x <- fried_sim$x[,1:5,drop = FALSE]
x_new <- x
y <- fried_sim$y

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)

# Testing the GP-BART
bart_test <- rbart(x_train = x,y = unlist(c(y)),x_test = x_test,
                   n_tree = 10,n_mcmc = 2500,alpha = 0.95,dif_order = 2,
                   beta = 2,nIknots = 20,delta = 1,
                   n_burn = 500,scale_bool = TRUE)

# Convergence plots
par(mfrow = c(3,1))
plot(bart_test$tau_post[-2500],type = "l", main = expression(tau),ylab=  "")
plot(bart_test$tau_b_post[-2500],type = "l", main = expression(tau[b]),ylab=  "")
plot(bart_test$tau_b_post_intercept[-2500],type = "l", main = expression(tau[b[0]]),ylab=  "")

bartmod <- dbarts::bart(x.train = x,y.train = unlist(c(y)),ntree = 200,x.test = x_test)


# Getting the splines model
library(mgcv)
splinemod <- gam(y ~ s(x, bs = "tp"), data = data.frame(y = y, x = x))
pred_spline <- predict(splinemod,newdata = data.frame(x = x_test))

# library(MOTRbart)
# motrbartmod <- MOTRbart::motr_bart(x = cbind(1,x),y = c(y))
# pred_motrbart = predict_motr_bart(motrbartmod,cbind(1,x_test),type = "all")

# par(mfrow=c(2,1))
# plot(y$x,bart_test$y_hat %>% rowMeans())
# plot(y$x,bartmod$yhat.train.mean)

# plot(x$x,bart_test$y_hat %>% rowMeans(), col = "blue")
# plot(x$x,bartmod$yhat.train.mean, col = "red")
# plot(x$x,pred_spline, col = "darkgreen")

# All ll trees prediction plot
all_tree_posterior_mean <- Reduce("+",bart_test$all_tree_post)/length(bart_test$all_tree_post)
if(!(ncol(all_tree_posterior_mean) %>% is_null())){
     colnames(all_tree_posterior_mean) <- paste0("tree.",1:ncol(all_tree_posterior_mean))
} else {
     all_tree_posterior_mean <- all_tree_posterior_mean %>% as.matrix()
     colnames(all_tree_posterior_mean) <- "tree.1"
}
all_tree_posterior_mean <- all_tree_posterior_mean %>% as.data.frame %>% add_column(x) %>% pivot_longer(starts_with("tree"))


# Getting quantiles for ribon
pi_ <- bart_test$y_hat_test %>% apply(1,function(x){quantile(x,probs = c(0.025,0.975))}) %>% t %>% cbind(x_new)
colnames(pi_) <- c("lower","upper","x")
pi_ <- pi_ %>% as.data.frame()
# Replicating the same plot on ggplot
ggplot()+
     geom_ribbon(data = pi_,mapping = aes(x = x,ymin = lower, ymax = upper), col = NA, alpha = 0.1, lty = "dashed",fill = "blue")+
     geom_point(data = data.frame(x = x, y = y ), mapping = aes(x = x, y =y ),alpha = 0.1)+
     # geom_point(data = data.frame(x = x, y = apply(bart_test[[1]],1,mean)), mapping = aes(x = x,  y = y), col = "darkblue", alpha = 0.7, pch= 3)+
     geom_line(data = data.frame(x = x_new, y = apply(bart_test[[2]],1,mean)), mapping = aes(x = x, y = y) , col = "blue") +
     geom_line(data = data.frame(x = x_new, y = bartmod$yhat.test.mean), mapping = aes(x = x, y = y), col ="red")+
     # geom_line(data = data.frame(x = x_new, y = pred_motrbart %>% colMeans()), mapping = aes(x = x, y = y), col ="darkgreen")+
     # geom_line(data = data.frame(x = x_new, y = unlist(pred_spline) ), mapping = aes(x = x, y = y), col ="darkgreen")+
     geom_line(data = all_tree_posterior_mean,
               mapping = aes(x = x, y = value, col = name), alpha = 0.5,show.legend = FALSE)+
     # ylim(y = range(y)*2)+
     theme_bw()

# plot(bart_test$tau_post, type = "l")
# microbenchmark::microbenchmark(bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 5000,
#                                     n_burn = 0,tau = 1,mu = 1,
#                                     tau_mu = 4*4*200,naive_sigma = 1,alpha = 0.95,
#                                     beta = 2,a_tau = 1,d_tau = 1,nsigma = 1),
#                                dbarts::bart(x.train = x,y.train = y,x.test = ,ntree = 200,ndpost = 5000,nskip = 0),times = 10)
# plot(bart_test$tau_b_post[-length(bart_test$tau_b_post)], type = "l")


# Traceplots
# par(mfrow=c(1,2))
# plot(bart_test$tau_b_post,type = "l", main = expression(tau[b]))
# plot(bart_test$tau_post,type = "l", main = expression(tau))


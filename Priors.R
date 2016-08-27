library(boot)
library(ggplot2)

alpha_mu<-rnorm(1e5,0,1/sqrt(0.01))
beta_mu<-rnorm(1e5,0,1/sqrt(0.1))

qplot(alpha_mu)
qplot(beta_mu)

#cut the tail off
fprior<-exp(0 + beta_mu * sample(seq(0,60,10),1e5,replace=T))

#just look at reasonable range, can't be less than 50 visits

fprior
qplot(fprior)

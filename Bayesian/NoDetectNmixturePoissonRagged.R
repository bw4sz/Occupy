
sink("Bayesian/NoDetectNmixturePoissonRagged.jags")

cat("
    model {
    for (x in 1:Nobs){
    
    # Covariates for true state   
    log(lambda[Bird[x],Plant[x],Time[x]]) <- alpha[Bird[x]] + beta[Bird[x]] * traitmatch[x] 
    
    #True State
    N[x] ~ dpois(lambda[Bird[x],Plant[x],Time[x]] )    
    
    #Observation Process
    Yobs[x] ~ dbin(detect[Bird[x]],N[x])    
    
    #Assess Model Fit
    
    #Fit discrepancy statistics
    eval[x]<-detect[Bird[x]]*N[x]
    E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    ynew[x]~dbin(detect[Bird[x]],N[x])
    E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    
    for (i in 1:Birds){
    detect[i] ~ dnorm(dprior,tau_detect) 
    alpha[i] ~ dnorm(intercept,tau_alpha)
    beta[i] ~ dnorm(gamma,tau_beta)    
    }
    
    #Hyperpriors
    #Slope grouping
    gamma~dnorm(0,0.0001)
    
    #Intercept grouping
    intercept~dnorm(0,0.0001)
    
    #detection prior
    dprior~dunif(.99,1)

    # Group intercept variance
    tau_alpha ~ dgamma(0.0001,0.0001)
    sigma_int<-pow(1/tau_alpha,0.5) #Derived Quantity
    
    #Slope variance, turning precision to sd
    tau_beta ~ dgamma(0.0001,0.0001)
    sigma_slope<-pow(1/tau_beta,0.5)
    
    #detect variance, turning precision to sd
    tau_detect ~ dgamma(0.01,0.01)
    sigma_detect<-pow(1/tau_detect,0.5)

    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

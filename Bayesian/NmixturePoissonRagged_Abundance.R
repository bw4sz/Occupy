
sink("Bayesian/NmixturePoissonRagged_Abundance.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Times){
    
    #Process Model with log normal overdispersion
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * resources[i,j,k] + epsilon[i,j,k]
    
    #variance
    epsilon[i,j,k] ~ dnorm(0,tauE)

    #For each Time - there is a latent count, log transformed.
    N[i,j,k] ~ dpois(lambda[i,j,k])
    }
    }
    }
    
    
    #Observed counts for each day of sampling at that Time
    for (x in 1:Nobs){
    
    #Observation Process
    Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Time[x]])    
    
    #Assess Model Fit
    
    #Fit discrepancy statistics
    eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Time[x]]
    E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Time[x]])
    E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    
    #Priors
    #Observation model
    #Detect priors, logit transformed - Following lunn 2012 p85
    
    for(x in 1:Birds){
    #For Cameras
    logit(detect[x])<-dcam[x]
    dcam[x]~dnorm(omega,omega_tau)
    }
    

    #Process Model
    #Species level priors
    for (i in 1:Birds){
    
    #Intercept
    alpha[i] ~ dnorm(alpha_mu,alpha_tau)
    
    #Traits slope 
    beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
    }
    
    #Group Detection Prior
    omega<-dnorm(0,0.386)
    omega_tau ~ dt(0,1,1)I(0,)

    #Group process priors
    
    #Intercept 
    alpha_mu ~ dnorm(0,0.01)
    alpha_tau ~ dt(0,1,1)I(0,)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #SLope
    beta1_mu~dnorm(0,0.01)
    beta1_tau ~ dt(0,1,1)I(0,)
    beta1_sigma<-pow(1/beta1_tau,0.5)

    #Overdispersion
    tauSigma ~ dunif(0,0.5)
    tauE <- pow(1/tauSigma,2)

    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

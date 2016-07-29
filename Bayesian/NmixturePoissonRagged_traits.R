
sink("Bayesian/NmixturePoissonRagged_traits.jags")

cat("
    model {
    #Compute true state for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Times){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]
    
    #True State
    N[i,j,k] ~ dpois(lambda[i,j,k])
    }
    }
    }
    
    #Observation Model
    for (x in 1:Nobs){
      
      #Observation Process for cameras
      Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Time[x]])    
      
      #     #Assess Model Fit - Posterior predictive check
      eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Time[x]]
      E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
      #     
      ynew[x]~dbin(detect[Bird[x]], N[Bird[x],Plant[x],Time[x]])
      E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    }
    
    #Priors
    #Observation model
    #Detect priors, logit transformed - Following lunn 2012 p85
    
    for(x in 1:Birds){
    #For Cameras
      logit(detect[x])<-dcam[x]
      dcam[x] ~ dnorm(dprior,tau_dcam)
    }
    
    #Detection group prior
    dprior ~ dnorm(0,0.386)
    
    #Group effect detect camera
    tau_dcam ~ dunif(0,1000)
    sigma_dcam<-pow(1/tau_dcam,.5)
    
    #Process Model
    #Species level priors
    for (i in 1:Birds){
    
    #Intercept
    alpha[i] ~ dnorm(alpha_mu,alpha_tau)
    
    #Traits slope 
    beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
    }
    
    #Group process priors
    
    #Intercept 
    alpha_mu ~ dnorm(0,0.386)
    alpha_tau ~ dunif(0,1000)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #Trait
    beta1_mu~dnorm(0,0.386)
    beta1_tau ~ dunif(0,1000)
    beta1_sigma<-pow(1/beta1_tau,0.5)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

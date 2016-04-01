
sink("Bayesian/NmixturePoissonRagged.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
      for (j in 1:Plants){
    
    #Process Model
      log(lambda[i,j])<-alpha[i] + beta[i] * Traitmatch[i,j]
      }
    }

    #For each camera - there is a latent count
    for(x in 1:Birds){
      for (y in 1:Plants){
        for (z in 1:Cameras){
          # true latent count
          N[x,y,z] ~ dpois(lambda[x,y])
        }
      }
    }

    #Observed counts for each day of sampling at that camera
    for (x in 1:Nobs){

    #Observation Process
    Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])    
    
    #Assess Model Fit
    
    #Fit discrepancy statistics
    eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Camera[x]]
    E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])
    E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    
    for (i in 1:Birds){
    logit(detect[i]) <- detect_logit[i]
    detect_logit[i] ~ dnorm(dprior,tau_detect)
    alpha[i] ~ dnorm(intercept,tau_alpha)
    beta[i] ~ dnorm(gamma,tau_beta)    
    }
    
    #Hyperpriors
    #Slope grouping
    gamma~dnorm(0,0.0001)
    
    #Intercept grouping
    intercept~dnorm(0,0.0001)

  # Group intercept variance
    sigma_int ~ dt(0,1,1)I(0,)
    tau_alpha <- pow(sigma_int,-2)

    #Detect grouping
    dprior ~ dnorm(0,.5)

  # Detect variance
    tau_detect ~ dunif(0,10)
    sigma_detect<-pow(1/tau_detect,.5) 
    
    #Derived Quantity
    #Slope variance, turning precision to sd
    sigma_slope ~ dt(0,1,1)I(0,)
    tau_beta <- pow(sigma_slope,-2)

    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

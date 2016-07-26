
sink("Bayesian/NmixturePoissonRagged_traits.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Cameras){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]
    
    #For each camera - there is a latent count
    N[i,j,k] ~ dpois(lambda[i,j,k])
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
    logit(detect[i]) <- dtrans[i]
    dtrans[i] ~ dnorm(dprior,tau_detect)
    alpha[i] ~ dnorm(alpha_mu,alpha_tau)
    beta1[i] ~ dnorm(beta1_mu,beta_tau)  
    }
    
    #Hyperpriors
    
    #Intercept grouping
    alpha_mu~dnorm(0,0.0001)

    # Group intercept variance
    alpha_tau ~ dgamma(0.0001,0.0001)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #Detect grouping
    dprior ~ dnorm(0,0.5)
    
    # Detect variance
    tau_detect ~ dunif(0,5)
    sigma_detect<-pow(1/tau_detect,0.5) 
    
    #Trait Slope
    beta1_mu~dnorm(0,0.0001)
    
    #Slope variance, turning precision to sd
    beta1_tau ~ dgamma(0.0001,0.0001)
    beta1_sigma<-pow(1/beta1_tau,0.5)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()


sink("Bayesian/NmixturePoissonRagged_Abundance.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Cameras){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Resources[i,j,k]
    
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
    alpha[i] ~ dnorm(intercept,tau_alpha)
    beta1[i] ~ dnorm(gamma1,tau_beta1)  
    }
    
    #Hyperpriors
    
    #Intercept grouping
    intercept~dnorm(0,0.0001)
    
    #Group intercept variance
    sigma_alpha ~ dt(0,10,1)I(0,)
    tau_alpha <- pow(sigma_alpha,-2)
    
    #Detect grouping
    dprior ~ dnorm(0,0.5)
    
    # Detect variance
    tau_detect ~ dunif(0,10)
    sigma_detect<-pow(1/tau_detect,0.5) 
    
    #Trait Slope
    
    #Mean
    gamma1 ~ dnorm(0,0.0001)
    
    #Variance
    sigma_beta1 ~ dt(0,10,1)I(0,)
    tau_beta1 <- pow(sigma_beta1,-2)
    
    #derived posterior check
    
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

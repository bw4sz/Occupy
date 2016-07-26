
sink("Bayesian/NoDetectNmixturePoissonRagged.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Cameras){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]
    gamma[i,j,k] <- beta2[i] * resources[i,j,k]
    }
    }
    }
    
    for (x in 1:Nobs){
    
      # Observed State   
      Yobs[x] ~ dpois(lambda[Bird[x],Plant[x],Camera[x]] * gamma[Bird[x],Plant[x],Camera[x]]) 
      
      #Assess Model Fit
      
      #Fit discrepancy statistics
      eval[x]<-lambda[Bird[x],Plant[x],Camera[x]]
      E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
      
      ynew[x]~dpois(lambda[Bird[x],Plant[x],Camera[x]])
      E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
      
      }
      
      for (i in 1:Birds){
      alpha[i] ~ dnorm(alpha_mu,alpha_tau)
      beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
      beta2[i] ~ dnorm(beta2_mu,beta2_tau)    
      }
      
    #Hyperpriors
    
    #Intercept grouping
    alpha_mu~dnorm(0,0.0001)
    
    # Group intercept variance
    alpha_sigma ~ dt(0,1,1)I(0,)
    alpha_tau <- pow(alpha_sigma,-2)
    
    #Detect grouping
    dprior ~ dnorm(0,.5)
    
    # Detect variance
    tau_detect ~ dunif(0,5)
    sigma_detect<-pow(1/tau_detect,.5) 
    
    #Trait Slope
    
    #Mean
    beta1_mu~dnorm(0,0.0001)
    
    #Variance
    beta1_sigma ~ dt(0,1,1)I(0,)
    beta1_tau <- pow(beta1_sigma,-2)
    
    #Abundance slope
    
    #Mean
    beta2_mu~dnorm(0,0.0001)
    
    beta2_sigma ~ dt(0,1,1)I(0,)
    beta2_tau <- pow(beta2_sigma,-2)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

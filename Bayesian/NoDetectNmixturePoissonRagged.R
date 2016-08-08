
sink("Bayesian/NoDetectNmixturePoissonRagged.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Cameras){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]
    }
    }
    }
    
    for (x in 1:Nobs){
      
      # Observed State   
      Yobs[x] ~ dpois(lambda[Bird[x],Plant[x],Camera[x]]) 
      
      #Assess Model Fit
      
      #Fit discrepancy statistics
      eval[x]<-lambda[Bird[x],Plant[x],Camera[x]]
      E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
      
      ynew[x]~dpois(lambda[Bird[x],Plant[x],Camera[x]])
      E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
      
      }
      
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
    alpha_mu ~ dnorm(0,0.0001)
    alpha_tau ~ dt(0,1,1)I(0,)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #Trait
    beta1_mu~dnorm(0,0.0001)
    beta1_tau ~ dt(0,1,1)I(0,)
    beta1_sigma<-pow(1/beta1_tau,0.5)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    
    }
    ",fill=TRUE)

sink()

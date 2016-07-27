
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
      
      #expected intensity
      N[x]<-lambda[Bird[x],Plant[x],Camera[x]] * resources[Bird[x],Plant[x],Camera[x]]
      
      # Observed State   
      Yobs[x] ~ dpois(N[x]+0.00000001) 
      
      #Assess Model Fit
      
      #Fit discrepancy statistics
      eval[x]<-lambda[Bird[x],Plant[x],Camera[x]] * resources[Bird[x],Plant[x],Camera[x]]
      E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
      
      ynew[x]~dpois(lambda[Bird[x],Plant[x],Camera[x]])
      E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
      
      }
      
      for (i in 1:Birds){
      alpha[i] ~ dnorm(alpha_mu,alpha_tau)
      beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
      }
      
    #Hyperpriors
    
    #Intercept grouping
    alpha_mu~dnorm(0,0.0001)
    
    # Group intercept variance
    alpha_sigma ~ dt(0,1,1)I(0,)
    alpha_tau <- pow(alpha_sigma,-2)
    
    #Trait Slope
    #Mean
    beta1_mu~dnorm(0,0.0001)
    
    #Variance
    beta1_sigma ~ dt(0,1,1)I(0,)
    beta1_tau <- pow(beta1_sigma,-2)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

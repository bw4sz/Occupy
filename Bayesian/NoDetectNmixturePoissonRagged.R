
sink("Bayesian/NoDetectNmixturePoissonRagged.jags")

cat("
    model {
      #Compute intensity for each pair of birds and plants
      for (i in 1:Birds){
        for (j in 1:Plants){
        
        #Process Model
        log(lambda[i,j])<-alpha[i] + beta[i] * Traitmatch[i,j]
      }
    }
    
    for (x in 1:Nobs){
    
      # Covariates for true state   
      Yobs[x] ~ dpois(lambda[Bird[x],Plant[x]])    
      
      #Assess Model Fit
      
      #Fit discrepancy statistics
      eval[x]<-lambda[Bird[x],Plant[x]]
      E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
      
      ynew[x]~dpois(lambda[Bird[x],Plant[x]])
      E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
      
      }
      
      for (i in 1:Birds){
      alpha[i] ~ dnorm(intercept,tau_alpha)
      beta[i] ~ dnorm(gamma,tau_beta)    
      }
      
      #Hyperpriors
      #Slope grouping
      gamma~dnorm(0,0.0001)
      
      #Intercept grouping
      intercept~dnorm(0,0.0001)
      dprior~dnorm(0,0.5)
      
      # Group intercept variance
      tau_alpha ~ dgamma(0.0001,0.0001)
      sigma_int<-pow(1/tau_alpha,0.5) 
      
      #Derived Quantity
      
      #Slope variance, turning precision to sd
      tau_beta ~ dgamma(0.0001,0.0001)
      sigma_slope<-pow(1/tau_beta,0.5)
      
      #derived posterior check
      fit<-sum(E[]) #Discrepancy for the observed data
      fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

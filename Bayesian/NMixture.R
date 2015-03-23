
sink("Bayesian/NMixture.jags")

cat("
    model {
    for (i in 1:Birds){
      for (j in 1:Plants){
      
      # True state model for the only partially observed true state    
      log(lambda[i,j])<-alpha[i] + beta[i] * traitmatch[i,j] + polybeta[i] * pow(traitmatch[i,j],2)
      N[i,j] ~ dpois(lambda[i,j])
      
      for (k in 1:Months) {   
      # Observation model for the actual observations
      Y[i,j,k] ~ dbin(detect[i],N[i,j]) 
          }
      }
    }
    
    for (i in 1:Birds){
    detect[i] ~ dunif(0,1) # Detection for each bird species
    alpha[i] ~ dnorm(intercept,tau_alpha)
    beta[i] ~ dnorm(gamma,tau_beta)    
    polybeta[i] ~ dnorm(polygamma,tau_polybeta)    
    }
    
    #Hyperpriors
    gamma~dnorm(0.001,0.001)
    polygamma~dnorm(0.001,0.001)
    intercept~dnorm(0.001,0.001)
    
    tau_alpha ~ dgamma(0.001,0.001)
    sigma_int<-pow(1/tau_alpha,0.5) #Derived Quantity
    tau_beta ~ dgamma(0.001,0.001)
    sigma_slope<-pow(1/tau_beta,0.5)
    
    tau_polybeta ~ dgamma(0.001,0.001)
    sigma_polyslope<-pow(1/tau_polybeta,0.5)
    
    }
    ",fill=TRUE)

sink()

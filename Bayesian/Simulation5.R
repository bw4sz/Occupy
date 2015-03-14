setwd("C:/Users/Ben/Documents/Occupy/Bayesian")

sink("Simulation5.jags")

cat("
    model {
    
    for (i in 1:Birds){
      for (j in 1:Plants){
    
    # True state model for the only partially observed true state    
    logit(occ[i,j])<-alpha + beta * traitmatch[i,j]
    present[i,j] ~ dbern(occ[i,j])
    
    for (k in 1:Months) {   
    # Observation model for the actual observations
    sightp[i,j,k] <- present[i,j] * detect[i]
    Y[i,j,k] ~ dbern(sightp[i,j,k]) 
    }
    }
    }
    
    for (i in 1:Birds){
    detect[i] ~ dunif(0,1) # Detection for each bird species
    
    #Derived
    #transform
    }

    alpha ~ dnorm(0.001,0.001)
    beta ~ dnorm(0.001,0.001)
    
    }
    ",fill=TRUE)

sink()

setwd("C:/Users/Ben/Documents/Occupy/Bayesian")

sink("Simulation2.jags")

cat("
    model {
    
    for (i in 1:Birds){
    
      for (j in 1:Plants){

        # True state model for the only partially observed true state    
        presentp[i,j] <- occ[i,j]
        present[i,j] ~ dbern(presentp[i,j]) 
      
        for (k in 1:Months) {   

          # Observation model for the actual observations
          sightp[i,j,k] <- present[i,j] * detect[i]
          Y[i,j,k] ~ dbern(sightp[i,j,k]) 
          }
      }
    }
    
    for (i in 1:Birds){
      detect[i] ~ dunif(0,1) # Detection for each bird species
      for (j in 1:Plants){
        occ[i,j] ~ dunif(0,1) # Occupancy for each Bird-Plant combination
        }
      }
    }
    
    ",fill=TRUE)

sink()

setwd("C:/Users/Ben/Documents/Occupy/Bayesian")

sink("Simulation1.jags")

cat("
    model {
    
      for (i in 1:Birds){
        
        for (j in 1:Plants){
        # True state model for the only partially observed true state    
          presentp[i,j] <- occ[i]
          present[i,j] ~ dbern(presentp[i,j]) 

          for (k in 1:Months) {    
            # Observation model for the actual observations
            sightp[i,j,k] <- present[i,j] * detect[i]
            Y[i,j,k] ~ dbern(sightp[i,j,k]) 
          }
        }
      }
    
    for (i in 1:Birds){
    occ[i] ~ dunif(0,1)
    detect[i] ~ dunif(0,1)
    }
    
    }
    
    ",fill=TRUE)

sink()

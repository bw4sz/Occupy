setwd("C:/Users/Ben/Documents/Occupy/Bayesian")

sink("Simulation.jags")

cat("
    model {
    
    for (b in 1:Birds){
    
      for (i in 1:Plants){
        # True state model for the only partially observed true state    
        presentp[i,b] <- occ[b]
        present[i,b] ~ dbern(presentp[i,b]) # True occupancy z at site 
        
        for (j in 1:Months) {    
          # Observation model for the actual observations
          sightp[((b-1)*Plants)+i,j] <- present[i,b] * detect[b]
          Y[((b-1)*Plants)+i,j] ~ dbern(sightp[((b-1)*Plants)+i,j]) 
        }}
      }
    
    for (b in 1:Birds){
    occ[b] ~ dunif(0,1)
    detect[b] ~ dunif(0,1)
    }
    
    }
    
    ",fill=TRUE)

sink()
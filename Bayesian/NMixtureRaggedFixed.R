
sink("Bayesian/NMixtureRaggedFixed.jags")

cat("
    model {
    
  for (i in 1:Birds){
    for (j in 1:Plants){
      #Process Model
      log(lambda[i,j])<-alpha[i] + beta[i] * traitmatch[i,j]

      #True state model  
      N[i,j] ~ dpois(lambda[i,j])
      }
    }

    #Observation Model
    for (i in 1:Nobs){
      Y[i] ~ dbin(detect[Bird[i]],N[Bird[i],Plant[i]]) 
      
      #Fit discrepancy statistics
      eval[i]<-detect[Bird[i]]*N[Bird[i],Plant[i]]
      E[i]<-pow((Y[i]-eval[i]),2)/(eval[i]+0.5)

      y.new[i]~dbin(detect[Bird[i]],N[Bird[i],Plant[i]])
      E.new[i]<-pow((y.new[i]-eval[i]),2)/(eval[i]+0.5)
    }
    
    for (k in 1:Birds){
    detect[k] ~ dunif(0,1) # Detection for each bird species
    alpha[k] ~ dnorm(0.001,0.001)
    beta[k] ~ dnorm(0.001,0.001)    
    }
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()

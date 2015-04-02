
sink("Bayesian/NMixture.jags")

cat("
    model {
    for (i in 1:Birds){
    for (j in 1:Plants){    
    for (k in 1:Months) {   
    
    #True state model for the only partially observed true state    
    log(lambda[i,j])<-alpha + beta * traitmatch[i,j]
    Y[i,j,k] ~ dpois(lambda[i,j])
    
    #Fit discrepancy statistics
    eval[i,j,k]<-detect[i]*N[i,j]
    E[i,j,k]<-pow((Y[i,j,k]-eval[i,j,k]),2)/(eval[i,j,k]+0.5)
    
    y.new[i,j,k]~dbin(detect[i],N[i,j])
    E.new[i,j,k]<-pow((y.new[i,j,k]-eval[i,j,k]),2)/(eval[i,j,k]+0.5)
    }
    }
    }
    
    
    alpha ~ dnorm(0.001,0.001)
    beta ~ dnorm(0.001,0.01)    
      
    #derived posterior check
    fit<-sum(E[,,]) #Discrepancy for the observed data
    fitnew<-sum(E.new[,,])
    
    }
    ",fill=TRUE)

sink()

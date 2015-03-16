library(R2jags)

### 21.2. Data generation
n.site <- 200
vege <- sort(runif(n = n.site, min = -1.5, max =1.5))

alpha.lam <- 2  			# Intercept
beta1.lam <- 2				# Linear effect of vegetation
beta2.lam <- -2				# Quadratic effect of vegetation
lam <- exp(alpha.lam + beta1.lam * vege + beta2.lam * (vege^2))

par(mfrow = c(2,1))
plot(vege, lam, main = "", xlab = "", ylab = "Expected abundance", las = 1)

N <- rpois(n = n.site, lambda = lam)
table(N)				# Distribution of abundances across sites
sum(N > 0) / n.site			# Empirical occupancy

plot(vege, N, main = "", xlab = "Vegetation cover", ylab = "Realized abundance")
points(vege, lam, type = "l", lwd = 2)

par(mfrow = c(2,1))
alpha.p <- 1				# Intercept
beta.p <- -4				# Linear effect of vegetation
det.prob <- exp(alpha.p + beta.p * vege) / (1 + exp(alpha.p + beta.p * vege))
plot(vege, det.prob, ylim = c(0,1), main = "", xlab = "", ylab = "Detection probability")

expected.count <- N * det.prob
plot(vege, expected.count, main = "", xlab = "Vegetation cover", ylab = 
       "Apparent abundance", ylim = c(0, max(N)), las = 1)
points(vege, lam, type = "l", col = "black", lwd = 2) # Truth

R <- n.site
T <- 3					# Number of replicate counts at each site
y <- array(dim = c(R, T))

for(j in 1:T){
  y[,j] <- rbinom(n = n.site, size = N, prob = det.prob)
}
y

sum(apply(y, 1, sum) > 0)		# Apparent distribution (proportion occupied sites)
sum(N > 0)				# True occupancy

C <- c(y)

site = 1:R				# ‘Short’ version of site covariate
site.p <- rep(site, T)			# ‘Long’ version of site covariate
vege.p <- rep(vege, T)			# ‘Long’ version of vegetation covariate
cbind(C, site.p, vege.p)		# Check that all went right

max.count <- apply(y, 1, max)
naive.analysis <- glm(max.count ~ vege + I(vege^2), family = poisson)
summary(naive.analysis)
lin.pred <- naive.analysis$coefficients[1] + naive.analysis$coefficients[2] * vege + 
  naive.analysis$coefficients[3] * (vege*vege)

par(mfrow = c(1,1))
plot(vege, max.count, main = "", xlab = "Vegetation cover", ylab = "Abundance or count", 
     ylim = c(0,max(N)), las = 1)
points(vege, lam, type = "l", col = "black", lwd = 2)
points(vege, exp(lin.pred), type = "l", col = "red", lwd = 2)



### 21.3. Analysis using WinBUGS
# Define model
sink("KeryBugs.jags")
cat("
    model {
    
    # Priors
    alpha.lam ~ dnorm(0, 0.1)
    beta1.lam ~ dnorm(0, 0.1)
    beta2.lam ~ dnorm(0, 0.1)
    alpha.p ~ dnorm(0, 0.1)
    beta.p ~ dnorm(0, 0.1)
    
    # Likelihood
    # Biological model for true abundance
    for (i in 1:R) {			# Loop over R sites
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.lam + beta1.lam * vege[i] + beta2.lam * vege2[i]
    }
    
    # Observation model for replicated counts
    for (i in 1:n) {			# Loop over all n observations
    C[i] ~ dbin(p[i], N[site.p[i]])
    logit(p[i]) <- alpha.p + beta.p * vege.p[i]
    }
    
    # Derived quantities
    totalN <- sum(N[])			# Estimate total population size across all sites
    
    }
    ",fill=TRUE)
sink()

###WINBUGS
setwd("C:/Users/Ben/Documents/Occupy/Bayesian")

runs<-40000

#Source model
source("KeryBugs.jags")

#print model
print.noquote(readLines("KeryBugs.txt"))

# Inits function
Nst <- apply(y, 1, max) + 1
inits <- function(){list(N = Nst, alpha.lam=rnorm(1, 0, 1), beta1.lam=rnorm(1, 0, 1), 
                         beta2.lam=rnorm(1, 0, 1), alpha.p=rnorm(1, 0, 1), beta.p=rnorm(1, 0, 1))}

# Parameters to estimate
params <- c("N", "totalN", "alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p")

#Data
Dat<-list(R = R, vege = vege, vege2 = vege2, n = n, C = C, site.p = site.p, vege.p = vege.p)

# MCMC settings
nc <- 3
nb <- 200
ni <- 1200
nt <- 5

#Jags

mkery = jags(inits=inits,
          n.chains=nc,
          model.file="KeryBugs.jags",
          working.directory=getwd(),
          data=Dat,
          parameters.to.save=params,
          n.thin=nt,
          n.iter=ni,
          n.burnin=nb,
          DIC=T)


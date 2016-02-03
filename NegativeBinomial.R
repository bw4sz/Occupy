library(reshape2)
library(foreach)
library(doSNOW)
library(chron)
library(ggplot2)
library(knitr)
library(R2jags)
library(dplyr)
library(stringr)
library(gridExtra)
library(boot)
library(picante)
library(GGally)
library(bipartite)
library(MASS)

set.seed(3)
#source functions

source("Functions.R")
paste("Run Completed at",Sys.time())

#Simulation   

### Parameters

# * 10 hummingbird species
# * Range of hummingbird bill sizes (in mm) ~ Pois(10)/10
# * Twenty plants
# * Range of corolla sizes (in mm) ~ Pois(15)/10
# * Mean frequeny ($\lambda$) for each hummingbird is drawn from U(0,10)  
# * Trait matching (minimizing Bill-Corolla difference) is drawn from a hierarcichal distribution
# $$log(\lambda)<-\alpha_i + \beta_i *traitmatch$$
#   $$\alpha=N(3,0.2)$$
#   $$\beta = N(1,0.2)$$
#   
#   * Imperfect detection 
# * $$ p_i = U(0,0.5) $$ 
#   * 36 camera replicates
# 
# **View simulated strength and form of trait matching **
  
  #Simulation Parameters
  
#Number of hummingbird species
h_species=10
plant_species=20
cameras<-20
days<-3

#set dispersion parameter, get from command line
args <- commandArgs(TRUE)
print(args)
dispersion<-as.numeric(args[1])

#Bill sizes
Bill<-rpois(h_species,10)

#Corolla sizes
Corolla<-rpois(plant_species,15)

#Subtract both and take absolute value, convert cm
traitmatch<-abs(sapply(Corolla,function(x) x - Bill)/10)

#regression slope
gamma<- -1
intercept<- 3
sigma_slope<- 0.2
sigma_intercept<- 0.2

detection= runif(h_species,0,0.5)
beta<-rnorm(h_species,gamma,sigma_slope)
alpha<-rnorm(h_species,intercept,sigma_intercept)

#Compute true interaction matrices

#for each species loop through and create a replicate dataframe
obs<-array(dim=c(h_species,plant_species,cameras,days))
lambda<-array(dim=c(h_species,plant_species))
N<-array(dim=c(h_species,plant_species,cameras))

#draw intensities
for(x in 1:h_species){
  for (y in 1:plant_species){
    #true lambda - interaction intensity
    lambda[x,y]<-exp(alpha[x] + beta[x] * traitmatch[x,y])
  }
}

#draw latent states
for(x in 1:h_species){
  for (y in 1:plant_species){
    for (z in 1:cameras){
      # true latent count
      N[x,y,z]<-rnbinom(n=1,mu=lambda[x,y],size=dispersion)
    }
  }
}

#Observed counts in each day
for(x in 1:h_species){
  for (y in 1:plant_species){
    for (z in 1:cameras){
      for (d in 1:days){
        #true detection rate of that observed count
        obs[x,y,z,d]<-rbinom(1,N[x,y,z],p=detection[x])
      }
    }
  }
}

##View correlation in simulated latent state

mdat<-melt(N)
colnames(mdat)<-c("Bird","Plant","Camera","Interactions")

traitmelt<-melt(traitmatch)
colnames(traitmelt)<-c("Bird","Plant","traitmatch")

mdat<-merge(mdat,traitmelt,c("Bird","Plant"))
ggplot(mdat,aes(x=traitmatch,y=Interactions,col=as.factor(Bird))) + geom_point() + geom_smooth(aes(group=1),method="glm",method.args=list(family="poisson") ) + labs(col="Bird") + xlab("Absolute value of Bill Length - Corolla Length ")

##View Detection Rates

obs.state<-melt(obs)
colnames(obs.state)<-c("Bird","Plant","Camera","Day","Yobs")
obs.state<-merge(mdat,obs.state,by=c("Bird","Plant","Camera"))
ggplot(obs.state,aes(x=Interactions,y=Yobs,col=Camera)) + geom_point() + theme_bw() + geom_abline() + coord_equal()

# Hierarcichal Occupancy Model

# For hummingbird i visiting plant j recorded by camera k on day d:
#   
#   $$ Y_{i,j,k,d} \sim Binom(N_{i,j,k},detect_i)$$
#   $$N_{i,j,k} \sim Pois(\lambda_{i,j}) $$
#   $$log(\lambda_{i,j})<-\alpha_i + \beta_i * abs(Bill_i - Corolla_i) $$
#   $$detect_i \sim U(0,0.5)$$     
#   
#   **Priors**
#   
#   $$\alpha_i \sim N(intercept,\tau_{\alpha})$$
#   $$\beta_i \sim N(\gamma,\tau_{\beta})$$
#   
#   **Hyperpriors**
#   $$gamma \sim N(0,0.0001)$$
#   $$intercept \sim N(0,0.0001)$$
#   
#   $$\tau_{\alpha} \sim Gamma(0.0001,0.0001)$$
#   $$\tau_\beta \sim Gamma(0.0001,0.0001)$$
#   
#   **Derived quantities**
#   
#   $$\sigma_{intercept} = \sqrt[2]{\frac{1}{\tau_\alpha}}$$
#   $$\sigma_{slope} = \sqrt[2]{\frac{1}{\tau_\beta}}$$
  
  #Simulated data without detection
  
#runs<-70000
runs<-1000

#trigger parallel
paralleljags<-T

#Source model
source("Bayesian/NoDetectNmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NoDetectNmixturePoissonRagged.R"))

if(paralleljags){
  
  #for parallel run
  Yobs=obs.state$Yobs
  Bird=obs.state$Bird
  Plant=obs.state$Plant
  Plants=max(obs.state$Plant)
  Camera=obs.state$Camera
  Day=obs.state$Day
  Traitmatch<-traitmatch
  #number of birds to iterate
  Birds=max(obs.state$Bird)
  Nobs=length(obs.state$Yobs)
  
  #A blank Y matrix - all present
  Ninit<-rep(max(obs.state$Yobs)+1,Nobs)
  
  #Inits
  InitStage <- function() {list(beta=rep(0.5,Birds),alpha=rep(0.5,Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 5   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  Dat<-list("Yobs","Bird","Plant","Plants","Traitmatch","Birds","Nobs","Ninit")
  
  sim_niave<-do.call(jags.parallel,list(Dat,InitStage,ParsStage,model.file="Bayesian/NoDetectNmixturePoissonRagged.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc))
  
} else {
  #Input Data
  Dat <- list(
    Yobs=obs.state$Yobs,
    Bird=obs.state$Bird,
    Plant=obs.state$Plant,
    Plants=max(obs.state$Plant),
    Time=obs.state$Time,
    Traitmatch=traitmatch,
    #number of birds to iterate
    Birds=max(obs.state$Bird),
    Nobs=length(obs.state$Yobs))
  
  #A blank Y matrix - all present
  Ninit<-rep(max(Dat$Yobs)+1,Dat$Nobs)
  
  #Inits
  InitStage <- function() {list(beta=rep(0.5,Dat$Birds),alpha=rep(0.5,Dat$Birds),dtrans=rep(0,Dat$Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 5   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  #Jags
  
  sim_niave <- jags(inits=InitStage,
                    n.chains=nc, model.file="Bayesian/NoDetectNmixturePoissonRagged.jags",
                    working.directory=getwd(),
                    data=Dat,
                    parameters.to.save=ParsStage,
                    n.thin=nt,
                    n.iter=ni,
                    n.burnin=nb,
                    DIC=T)
}

pars_niave<-extract_par(sim_niave,data=obs.state)
pars_niave$Model<-c("Poisson GLMM")

##Assess Convergence

ggplot(pars_niave[pars_niave$par %in% c("detect","alpha","beta"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")

ggplot(pars_niave[pars_niave$par %in% c("gamma","sigma_int","sigma_slope"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")

##Posteriors

###Posterior Distributions
p<-ggplot(pars_niave[pars_niave$par %in% c("alpha","beta"),],aes(x=estimate)) + geom_histogram() + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() + ggtitle("Species Posteriors")

#Add true values
tr<-melt(data.frame(species=1:length(detection),alpha=alpha,beta=beta),id.var='species')
colnames(tr)<-c("species","par","value")
psim<-p + geom_vline(data=tr,aes(xintercept=value),col='red',linetype='dashed',size=1)
psim

p<-ggplot(pars_niave[pars_niave$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=estimate)) + geom_histogram() + ggtitle("Hierarchical Posteriors") + facet_grid(~par,scale="free") + theme_bw()

#Add true values
tr<-melt(list(gamma=gamma,intercept=intercept,sigma_int=sigma_intercept,sigma_slope=sigma_slope))

colnames(tr)<-c("value","par")

psim2<-p + geom_vline(data=tr,aes(xintercept=value),linetype='dashed',size=1,col="red")

#**True values are given in the dashed lines.**
  
  ##Predicted Relationship 
  
castdf<-group_by(pars_niave,Chain) %>% dplyr::select(par,estimate) %>% filter(par %in% c("gamma","intercept"))

castdf<-dcast(pars_niave[pars_niave$par %in% c("gamma","intercept"),], Chain + Draw~par,value.var="estimate")

#calculated predicted y
predyniave<-trajF(alpha=castdf$intercept,beta=castdf$gamma,x=as.numeric(traitmatch))

orig<-trajF(alpha=rnorm(2000,intercept,sigma_intercept),beta=rnorm(2000,gamma,sigma_slope),x=as.numeric(traitmatch))

# Simulated data with detection

runs<-900

#trigger parallel
paralleljags<-T

#Source model
source("Bayesian/NmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NmixturePoissonRagged.R"))

if(paralleljags){
  
  #for parallel run
  Yobs=obs.state$Yobs
  Bird=obs.state$Bird
  Plant=obs.state$Plant
  Camera=obs.state$Camera
  Cameras=max(obs.state$Camera)
  Day<-obs.state$Day
  Days<-max(obs.state$Day)
  Traitmatch=traitmatch
  #number of birds to iterate
  Birds=max(obs.state$Bird)
  Plants=max(obs.state$Plant)
  Nobs=length(obs.state$Yobs)
  
  #A blank Y matrix - all present
  Ninit<-array(dim=c(h_species,plant_species,cameras),data=max(obs.state$Yobs)+1)
  
  #Inits
  InitStage <- function() {list(beta=rep(0.5,Birds),alpha=rep(0.5,Birds),dtrans=rep(0,Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 8   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  Dat<-list("Yobs","Bird","Plant","Plants","Traitmatch","Birds","Nobs","Ninit","Day","Days","Camera","Cameras")
  
  system.time(sim_detect<-do.call(jags.parallel,list(Dat,InitStage,ParsStage,model.file="Bayesian/NmixturePoissonRagged.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc)))
} else {
  #Input Data
  Dat <- list(
    Yobs=obs.state$Yobs,
    Bird=obs.state$Bird,
    Plant=obs.state$Plant,
    Camera=obs.state$Camera,
    Cameras=max(obs.state$Camera),
    Day<-obs.state$Day,
    Days<-max(obs.state$Day),
    Traitmatch=traitmatch,
    #number of birds to iterate
    Birds=max(obs.state$Bird),
    Plants=max(obs.state$Plant),
    Nobs=length(obs.state$Yobs)
  )
  #A blank Y matrix - all present
  Ninit<-array(dim=c(h_species,plant_species,cameras),data=max(obs.state$Yobs)+1)
  
  #Inits
  InitStage <- function() {list(beta=rep(0.5,Dat$Birds),alpha=rep(0.5,Dat$Birds),dtrans=rep(0,Dat$Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 5   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  #Jags
  
  system.time(sim_detect <- jags(inits=InitStage,
                                 n.chains=nc, model.file="Bayesian/NmixturePoissonRagged.jags",
                                 working.directory=getwd(),
                                 data=Dat,
                                 parameters.to.save=ParsStage,
                                 n.thin=nt,
                                 n.iter=ni,
                                 n.burnin=nb,
                                 DIC=T))
}

pars<-extract_par(sim_detect,data=obs.state)
pars$Model<-"Occupancy"

pO<-rbind_all(list(pars_niave,pars))
filename<-paste("Dispersion/Allpars",dispersion,".csv",sep="")
write.csv(pO,filename)

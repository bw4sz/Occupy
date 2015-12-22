# Hierarchical occupancy models for species interactions - Simulations
Ben Weinstein - Stony Brook University  




```
## [1] "Run Completed at 2015-12-21 19:21:49"
```

#Simulation   

### Parameters

* 10 hummingbird species
* Range of hummingbird bill sizes (mm) ~ Pois(2)
* Twenty plants
* Range of corolla sizes (mm) ~ Pois(2)
* Mean frequeny ($\lambda$) for each hummingbird is drawn from U(0,10)  
* Trait matching (minimizing Bill-Corolla difference) is drawn from a hierarcichal distribution
$$log(\lambda)<-\alpha_i + \beta_i *traitmatch$$
$$\alpha=N(3,0.01)$$
$$\beta = N(-0.01,0.001)$$

* Imperfect detection 
* $$ p_i = U(0,0.5) $$ 
* 36 camera replicates

**View simulated strength and form of trait matching **


```r
load("Abundance.RData")
```

#Simulation Parameters


```r
#Number of hummingbird species
h_species=10
plant_species=20
cameras<-20
days<-3

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
```

#Compute true interaction matrices


```r
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
      N[x,y,z]<-rpois(1,lambda[x,y])
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
```

##View correlation in simulated latent state


```r
mdat<-melt(N)
colnames(mdat)<-c("Bird","Plant","Camera","Interactions")

traitmelt<-melt(traitmatch)
colnames(traitmelt)<-c("Bird","Plant","traitmatch")

mdat<-merge(mdat,traitmelt,c("Bird","Plant"))
ggplot(mdat,aes(x=traitmatch,y=Interactions,col=as.factor(Bird))) + geom_point() + geom_smooth(aes(group=1),method="glm",family="poisson") + labs(col="Bird") + xlab("Absolute value of Bill Length - Corolla Length ")
```

<img src="figure/unnamed-chunk-6-1.png" title="" alt="" style="display: block; margin: auto;" />

##View Detection Rates


```r
obs.state<-melt(obs)
colnames(obs.state)<-c("Bird","Plant","Camera","Day","Yobs")
obs.state<-merge(mdat,obs.state,by=c("Bird","Plant","Camera"))
ggplot(obs.state,aes(x=Interactions,y=Yobs,col=Camera)) + geom_point() + theme_bw() + geom_abline() + coord_equal()
```

<img src="figure/unnamed-chunk-7-1.png" title="" alt="" style="display: block; margin: auto;" />

# Hierarcichal Occupancy Model

For hummingbird i visiting plant j recorded by camera k on day d:

$$ Y_{i,j,k,d} \sim Binom(N_{i,j,k},detect_i)$$
$$N_{i,j,k} \sim Pois(\lambda_{i,j}) $$
$$log(\lambda_{i,j})<-\alpha_i + \beta_i * abs(Bill_i - Corolla_i) $$
$$detect_i \sim U(0,0.5)$$     

**Priors**

$$\alpha_i \sim N(intercept,\tau_{\alpha})$$
$$\beta_i \sim N(\gamma,\tau_{\beta})$$

**Hyperpriors**
$$gamma \sim N(0,0.0001)$$
$$intercept \sim N(0,0.0001)$$

$$\tau_{\alpha} \sim Gamma(0.0001,0.0001)$$
$$\tau_\beta \sim Gamma(0.0001,0.0001)$$

**Derived quantities**

$$\sigma_{intercept} = \sqrt[2]{\frac{1}{\tau_\alpha}}$$
$$\sigma_{slope} = \sqrt[2]{\frac{1}{\tau_\beta}}$$
$$\sigma_{detect} = \sqrt[2]{\frac{1}{\tau_detect}}$$

#Simulated data without detection


```r
runs<-40000
#runs<-1000

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
```


```r
#recompile if needed
load.module("dic")
runs<-30000
recompile(sim_niave)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 86685
## 
## Initializing model
## 
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 86685
## 
## Initializing model
```

```r
sim_niave<-update(sim_niave,n.iter=runs,n.burnin=runs*.95)
```


```r
pars_niave<-extract_par(sim_niave,data=obs.state)
pars_niave$Model<-c("Poisson GLMM")
```

##Assess Convergence


```r
ggplot(pars_niave[pars_niave$par %in% c("detect","alpha","beta"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figure/unnamed-chunk-11-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
ggplot(pars_niave[pars_niave$par %in% c("gamma","sigma_int","sigma_slope"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figure/unnamed-chunk-12-1.png" title="" alt="" style="display: block; margin: auto;" />

##Posteriors


```r
###Posterior Distributions
p<-ggplot(pars_niave[pars_niave$par %in% c("alpha","beta"),],aes(x=estimate)) + geom_histogram() + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() + ggtitle("Species Posteriors")

#Add true values
tr<-melt(data.frame(species=1:length(detection),alpha=alpha,beta=beta),id.var='species')
colnames(tr)<-c("species","par","value")
psim<-p + geom_vline(data=tr,aes(xintercept=value),col='red',linetype='dashed',size=1)
psim
```

<img src="figure/unnamed-chunk-13-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/SimulationPosteriorsNoDetect.jpg",dpi=300,height=8,width=8)
```


```r
p<-ggplot(pars_niave[pars_niave$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=estimate)) + geom_histogram() + ggtitle("Hierarchical Posteriors") + facet_grid(~par,scale="free") + theme_bw()

#Add true values
tr<-melt(list(gamma=gamma,intercept=intercept,sigma_int=sigma_intercept,sigma_slope=sigma_slope))

colnames(tr)<-c("value","par")

psim2<-p + geom_vline(data=tr,aes(xintercept=value),linetype='dashed',size=1,col="red")
```

**True values are given in the dashed lines.**

##Predicted Relationship 


```r
castdf<-group_by(pars_niave,Chain) %>% select(par,estimate) %>% filter(par %in% c("gamma","intercept"))

castdf<-dcast(pars_niave[pars_niave$par %in% c("gamma","intercept"),], Chain + Draw~par,value.var="estimate")

#calculated predicted y
predyniave<-trajF(alpha=castdf$intercept,beta=castdf$gamma,x=as.numeric(traitmatch))

orig<-trajF(alpha=rnorm(2000,intercept,sigma_intercept),beta=rnorm(2000,gamma,sigma_slope),x=as.numeric(traitmatch))
```

# Simulated data with detection


```r
runs<-60000

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
```


```r
#recompile if needed
load.module("dic")
runs<-10000
recompile(sim_detect)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 135200
## 
## Initializing model
## 
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 135200
## 
## Initializing model
```

```r
sim_detect<-update(sim_detect,n.iter=runs,n.burnin=runs*.95,n.thin=5)
```


```r
pars<-extract_par(sim_detect,data=obs.state)
pars$Model<-"Occupancy"
```

##Assess Convergence


```r
ggplot(pars[pars$par %in% c("detect","alpha","beta"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figure/unnamed-chunk-19-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
ggplot(pars[pars$par %in% c("gamma","sigma_int","sigma_slope"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figure/unnamed-chunk-20-1.png" title="" alt="" style="display: block; margin: auto;" />

##Posteriors


```r
###Posterior Distributions
p<-ggplot(pars[pars$par %in% c("detect","alpha","beta"),],aes(x=estimate)) + geom_histogram() + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() + ggtitle("Species Posteriors")

#Add true values
tr<-melt(data.frame(species=1:length(detection),detect=detection,alpha=alpha,beta=beta),id.var='species')
colnames(tr)<-c("species","par","value")
psim<-p + geom_vline(data=tr,aes(xintercept=value),col='red',linetype='dashed',size=1)
#ggsave("Figures/SimulationPosteriors.jpg",dpi=300,height=8,width=8)
```


```r
p<-ggplot(pars[pars$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=estimate)) + geom_histogram() + ggtitle("Hierarchical Posteriors") + facet_wrap(~par,scale="free",nrow=2) + theme_bw() 

#Add true values
tr<-melt(list(gamma=gamma,intercept=intercept,sigma_int=sigma_intercept,sigma_slope=sigma_slope))

colnames(tr)<-c("value","par")

psim2<-p + geom_vline(data=tr,aes(xintercept=value),linetype='dashed',size=1,col="black")
#ggsave("Figures/SimulationH.jpg",dpi=300,height=4,width=10)
```

**True values are given in the dashed lines.**

##Compare simulation posteriors with and without detection


```r
#Bind to other dataset
parsall<-rbind.data.frame(pars[!pars$par %in% "ynew",],pars_niave[!pars_niave$par %in% "ynew",])
parsall$Model<-as.factor(parsall$Model)

###Posterior Distributions
p<-ggplot(parsall[parsall$par %in% c("detect","alpha","beta"),],aes(x=estimate,fill=Model)) + geom_histogram(position="identity") + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() 

#Add true values
tr<-melt(data.frame(species=1:length(detection),detect=detection,alpha=alpha,beta=beta),id.var='species')
colnames(tr)<-c("species","par","value")
psim<-p + geom_vline(data=tr,aes(xintercept=value),col='black',linetype='dashed',size=1)
psim
```

<img src="figure/unnamed-chunk-23-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
#ggsave("Figures/SimulationPosteriorsBoth.jpg",dpi=300,height=8,width=8)
```


```r
p<-ggplot(parsall[parsall$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=estimate,fill=Model)) + geom_histogram(position="identity") + ggtitle("Hierarchical Posteriors") + facet_wrap(~par,scale="free",nrow=2) + theme_bw() 

#Add true values
tr<-melt(list(gamma=gamma,intercept=intercept,sigma_int=sigma_intercept,sigma_slope=sigma_slope))

colnames(tr)<-c("value","par")

psim2<-p + geom_vline(data=tr,aes(xintercept=value),linetype='dashed',size=1,col="black")
psim2
```

<img src="figure/unnamed-chunk-24-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
#ggsave("Figures/SimulationHBoth.jpg",dpi=300,height=4,width=10)
```

##Correlation in posteriors for Occupancy Model


```r
castdf<- pars %>% filter(Model =="Occupancy") %>% group_by(Chain) %>% select(par,estimate,Draw) %>% filter(par %in% c("gamma","intercept")) %>% dcast(Chain+Draw~par,value.var="estimate")
head(castdf)
```

```
##   Chain Draw      gamma intercept
## 1     1    1 -1.0533785  2.981538
## 2     1    2 -1.0902280  2.800631
## 3     1    3 -0.9295978  2.900946
## 4     1    4 -0.9989492  2.986151
## 5     1    5 -0.8629959  2.946017
## 6     1    6 -0.9382878  2.915441
```

```r
ggpairs(castdf[,3:4],title="Correlation in Group-Level Posteriors")
```

<img src="figure/unnamed-chunk-25-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
castdf<- pars %>% filter(Model =="Occupancy") %>% group_by(Chain) %>% select(par,estimate,Draw,species) %>% filter(par %in% c("alpha","beta","detect")) %>% dcast(species+Chain+Draw~par,value.var="estimate")
head(castdf)
```

```
##   species Chain Draw    alpha       beta    detect
## 1       1     1    1 2.978369 -0.6186920 0.2722375
## 2       1     1    2 2.966052 -0.6728209 0.2797062
## 3       1     1    3 2.977678 -0.6534394 0.2730205
## 4       1     1    4 2.951012 -0.5925107 0.2717407
## 5       1     1    5 2.921249 -0.6008911 0.2777064
## 6       1     1    6 3.014426 -0.7611605 0.2942833
```

```r
ggpairs(castdf[,4:6],title="Correlation in Species-Level Posteriors")
```

<img src="figure/unnamed-chunk-25-2.png" title="" alt="" style="display: block; margin: auto;" />

##Predicted Relationship 

Does accounting for non-independence and detection change our estimate of trait matching?


```r
castdf<-group_by(pars,Chain) %>% select(par,estimate) %>% filter(par %in% c("gamma","intercept"))

castdf<-dcast(pars[pars$par %in% c("gamma","intercept"),], Chain + Draw~par,value.var="estimate")

trajF<-function(alpha,beta,x){
  indat<-cbind(alpha,beta)
  
  #fit regression for each input estimate
  sampletraj<-apply(indat,1,function(s){
    data.frame(x=x,y=exp(s['alpha'] + s['beta'] *x))})
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
  }

#calculated predicted y
predy<-trajF(alpha=castdf$intercept,beta=castdf$gamma,x=as.numeric(traitmatch))

orig<-trajF(alpha=rnorm(2000,intercept,sigma_intercept),beta=rnorm(2000,gamma,sigma_slope),x=as.numeric(traitmatch))


#plot and compare to original data
psim3<-ggplot(data=predy,aes(x=x)) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.1,fill="red")  + geom_line(aes(y=mean),size=.8,col="red",linetype="dashed") + theme_bw() + ylab("Interactions") + geom_line(data=orig,aes(x=x,y=mean),col='black',size=1) + geom_ribbon(data=predyniave,aes(ymin=lower,ymax=upper),alpha=0.1,fill="blue") +  geom_line(data=predyniave,aes(y=mean),size=.8,col="blue",linetype="dashed") + xlab("Difference between Bill and Corolla Length") 

psim3
```

<img src="figure/unnamed-chunk-26-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/SimulationResults.jpg",height=5,width=6,dpi=300)
```

Black line is the true relationship. The red line is the posterior mean with confidible intervals in shaded grey for the proposed bayesian model. The blue line is the same model, but assuming perfect detection rates.

**Conclusion:** Accounting for detection and non-independence greatly increases the accuracy of the predicted state. The perfect detection model underestimates the strength of trait matching among hummingbirds and their foodplants.

##Posterior Predictive Check

Since I have simualted the data, it should fit as well as any random dataset drawn from the estimated parameters. An ideal fit would be posterior values sitting along the 1:1 line.


```r
fitstat<-droplevels(parsall[parsall$par %in% c("fit","fitnew"),])
fitstat<-dcast(fitstat,Draw+Chain+Model~par,value.var="estimate")

#add 1:1 line
ymin<-round(min(c(fitstat$fit,fitstat$fitnew)))
ymax<-round(max(c(fitstat$fit,fitstat$fitnew)))
ab<-data.frame(x=ymin:ymax,y=ymin:ymax)
p<-ggplot(fitstat,aes(x=fit,y=fitnew)) + geom_point(aes(col=Model)) + theme_bw() + coord_equal()
psim4<-p  + labs(x="Discrepancy of observed data",y="Discrepancy of replicated data",col="Model") + geom_line(data=ab,aes(x=x,y=y)) + ggtitle("Simulated Data")
psim4
```

<img src="figure/unnamed-chunk-27-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/SimulationDisc.jpeg",height=5,width=5)
```


```
## png 
##   2
```

```
## png 
##   2
```

#Compare Occupancy Models to a Multinomial fit

Several studies follow Vasquez 2009 in fitting a multinomial relationship to frequency and interaction probabilities. See ??mgen in package bipartite. Since traitmatching is a distance, take 1-probabilities to correct weight the multinomial


```r
#true number of observed interactions
true_state<-obs.state %>% group_by(Bird,Plant) %>% summarize(n=sum(Yobs)) %>% acast(.$Bird~.$Plant)

m<-max(traitmatch)-traitmatch
m<-m/sum(m)

#True latent interactions, since the simulated data is known.
print(paste("Correlation coefficient is:", round(cor(c(true_state),c(m),method="spearman"),2)))
```

```
## [1] "Correlation coefficient is: 0.53"
```

What is the discrepancy of a multinomial approach?

Chisquared statistic is (observed-expected)^2/(expected + 0.5), we add the 0.5 to avoid dividing by 0. For each combination of birds and plants, predicted versus mean observed.


```r
#define discrep function
chisq<-function(o,e){(o-e)^2/(e+0.5)}

nullm<-function(){
r<-mgen(m,sum(true_state),keep.species = F)
}

#same number of draws as bayesian
cl<-makeCluster(20,"SOCK")
registerDoSNOW(cl)
mult_matrix<-foreach(x=1:length(pars[pars$par %in% "fit","estimate"]),.packages=c("bipartite","reshape2")) %dopar% nullm()
stopCluster(cl)

mats<-lapply(mult_matrix,function(r){
rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))

#for each position what is the chisq
for (x in 1:nrow(r)){
  for (y in 1:ncol(r)){
   rmerge[x,y]<-chisq(e=r[x,y],o=true_state[x,y])
  }
}

#sum of chisq driscrepancy
return(rmerge)
})

multi_disc<-sapply(mats,function(x) sum(x))

qplot(multi_disc)+ xlab("Chi-squared Discrepancy for Multimonial Liklihood") + geom_vline(xintercept=mean(multi_disc),col='red',linetype='dashed')
```

<img src="figure/unnamed-chunk-30-1.png" title="" alt="" style="display: block; margin: auto;" />

#Compare Bayesian Occupancy and Multinomial using true known interactions

## No Detection Occupancy Model


```r
N_niave<-pars_niave[ pars_niave$par %in% "ynew",]

bydraw<-split(N_niave,list(N_niave$Chain,N_niave$Draw))

occ_nodetect_matrix<-lapply(bydraw,function(x){
  r<-acast(x,species ~ plant,value.var = "estimate",fun.aggregate = sum)
})

#calculate discrep on those deviates
occ_nodetect<-lapply(occ_nodetect_matrix,function(r){
    #for each position what is the chisq
  rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))
  for (x in 1:nrow(r)){
    for (y in 1:ncol(r)){
     rmerge[x,y]<-chisq(e=r[x,y],o=true_state[x,y])
      }
    }
  return(rmerge)
})

names(occ_nodetect)<-1:length(occ_nodetect)
names(occ_nodetect_matrix)<-1:length(occ_nodetect_matrix)
```

##With Detection


```r
Ndetect<-pars[pars$par %in% "ynew",]
bydraw<-split(Ndetect,list(Ndetect$Chain,Ndetect$Draw))
occ_matrix<-lapply(bydraw,function(x){
  r<-acast(x,species ~ plant,value.var = "estimate",fun.aggregate = sum)
  })

#calculate discrep for those aggregated matrices
occ<-lapply(occ_matrix,function(r){
    #for each position what is the chisq
  rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))
  for (x in 1:nrow(r)){
    for (y in 1:ncol(r)){
     rmerge[x,y]<-chisq(e=r[x,y],o=true_state[x,y])
      }
    }
  return(rmerge)
})

names(occ)<-1:length(occ)
names(occ_matrix)<-1:length(occ_matrix)
```

##Compare to observed data


```r
mmat<-melt(true_state)
colnames(mmat)<-c("Bird","Plant","True_State")

#append to predicted matrices

#occupancy with detection
mocc<-melt(occ_matrix)
colnames(mocc)<-c("Bird","Plant","Occupancy","Iteration")
simdat<-merge(mocc,mmat,by=c("Bird","Plant"),all.x=T)

#occupancy with nodetection
moccd<-melt(occ_nodetect_matrix)
colnames(moccd)<-c("Bird","Plant","Poisson GLMM","Iteration")

simdat<-merge(simdat,moccd,by=c("Bird","Plant","Iteration"))

#multinomial
multmats<-melt(mult_matrix)
colnames(multmats)<-c("Bird","Plant","Multinomial","Iteration")
simdat<-merge(simdat,multmats,by=c("Bird","Plant","Iteration"))


simdat<-melt(simdat,measure.vars = c("Occupancy","Poisson GLMM","Multinomial"))

ggplot(simdat,aes(x=True_State,y=value,col=variable),alpha=.4) + geom_point() + geom_abline() + labs(col="Model") + ylab("Predicted State") + xlab("True State") + theme_bw() + facet_wrap(~variable)
```

<img src="figure/unnamed-chunk-33-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/PredictedState.jpeg",height=4,width=7)

#difference in the middle
simd<-dcast(simdat,...~variable)
simd$Diff<-simd$Occupancy-simd$`Poisson GLMM`
ggplot(simd,aes(x=True_State,y=Diff)) + geom_point() + ylab("Difference in Occupancy and Poisson GLMM")
```

<img src="figure/unnamed-chunk-33-2.png" title="" alt="" style="display: block; margin: auto;" />

## View predicted trait-matching relationship with multinomial included.

```r
#merge multinomial with trait relationship
multmats<-merge(multmats,traitmelt)

multline<-multmats %>% group_by(traitmatch) %>% summarize(y=mean(Multinomial),Upper=quantile(Multinomial,0.95),Lower=quantile(Multinomial,0.05)) 

psim3 + geom_ribbon(data=multline,aes(x=traitmatch,ymin=Lower,ymax=Upper),alpha=0.5,fill='gray') + geom_line(data=multline,aes(x=traitmatch,y=y),linetype='dashed')
```

<img src="figure/unnamed-chunk-34-1.png" title="" alt="" style="display: block; margin: auto;" />

View a couple example data points from across the type of interactions.


```r
h<-simdat[which.max(simdat$True_State),c("Bird","Plant")]
d<-simdat[simdat$Bird %in% h$Bird & simdat$Plant %in% h$Plant,]

ggplot(data=d,aes(x=value,fill=variable))+ geom_histogram(position="identity") + labs(fill="Model") + geom_vline(aes(xintercept=True_State)) + ggtitle("High Visitation Example")
```

<img src="figure/unnamed-chunk-35-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
h<-simdat[which.min(simdat$True_State),c("Bird","Plant")]
d<-simdat[simdat$Bird %in% h$Bird & simdat$Plant %in% h$Plant,]

ggplot(data=d,aes(x=value,fill=variable))+ geom_histogram(position="identity") + labs(fill="Model") + geom_vline(aes(xintercept=True_State)) + ggtitle("Low Visitation Example")
```

<img src="figure/unnamed-chunk-35-2.png" title="" alt="" style="display: block; margin: auto;" />

##Summary of discrepancy of predicted matrices


```r
#Multinomial
multi_disc<-sapply(mats,function(x) median(x))

#occupancy without detection
occno_disc<-sapply(occ_nodetect,function(x) median(x))

#Occupancy with detection
occ_disc<-sapply(occ,function(x) median(x))

#compared to bayesian
ggplot(data.frame(multi_disc)) + geom_histogram(aes(x=multi_disc),fill="blue",alpha=.6)+ xlab("Chi-squared Discrepancy") + geom_histogram(data=data.frame(occ_disc),aes(x=occ_disc),fill="red",alpha=.6) + theme_bw() +geom_vline(aes(xintercept=mean(occ_disc)),linetype="dashed",col="red")+ geom_vline(xintercept=mean(multi_disc),linetype="dashed",col="blue") + geom_histogram(data=data.frame(occno_disc),aes(x=occno_disc),fill="orange",alpha=.6) + geom_vline(aes(xintercept=mean(occno_disc)),linetype="dashed",col="orange")
```

<img src="figure/unnamed-chunk-36-1.png" title="" alt="" style="display: block; margin: auto;" />

##Comparison of summary statistics for all three approaches


```r
d<-list(Occupancy=occ,Multinomial=mats,Poisson_GLM=occ_nodetect)
d<-melt(d)
colnames(d)<-c("Bird","Plant","value","Iteration","Model")

d %>% group_by(Model,Iteration) %>% summarize(mean=mean(value),sd=sd(value),sum=sum(value)) %>% group_by(Model) %>% summarize(mean_mean=round(mean(mean),2),mean_sd=round(sd(mean),2),mean_sum=round(mean(sum),2))
```

```
## Source: local data frame [3 x 4]
## 
##         Model mean_mean mean_sd mean_sum
##         (chr)     (dbl)   (dbl)    (dbl)
## 1 Multinomial    184.47    1.56 36893.33
## 2   Occupancy      1.67    0.21   334.05
## 3 Poisson_GLM      2.65    0.26   530.23
```

#DIC


```r
sim_niave$BUGSoutput$DIC
```

```
## [1] 39459.54
```

```r
sim_detect$BUGSoutput$DIC
```

```
## [1] 41537.42
```



```r
save.image("Abundance.Rdata")
```

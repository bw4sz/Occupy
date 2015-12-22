# Hierarchical occupancy Models for species interactions: Empirical Data
Ben Weinstein - Stony Brook University  




```
## [1] "Run Completed at 2015-12-22 11:18:44"
```


```r
#reload if needed
load("Observed.Rdata")
```

#Observed dataset


```r
#read in flower morphology data, comes from Nectar.R
droppath<-"C:/Users/Ben/Dropbox/"
fl.morph<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/FlowerMorphology.csv",sep=""))

#First row is empty
fl.morph<-fl.morph[-1,]

#Bring in Hummingbird Morphology Dataset, comes from
hum.morph<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/HummingbirdMorphology.csv",sep=""))

#Bring in Interaction Matrix
int<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/Network/HummingbirdInteractions.csv",sep=""),row.names=1)

#just use camera data
int<-int[is.na(int$TransectID),]

#Correct known taxonomic disagreements, atleast compared to traits
int[int$Iplant_Double=="Alloplectus purpureus","Iplant_Double"]<-"Glossoloma purpureum"
int[int$Iplant_Double=="Capanea affinis","Iplant_Double"]<-"Kohleria affinis"
int[int$Iplant_Double=="Columnea cinerea","Iplant_Double"]<-"Columnea mastersonii"
int[int$Iplant_Double=="Alloplectus teuscheri","Iplant_Double"]<-"Drymonia teuscheri"
int[int$Iplant_Double=="Drymonia collegarum","Iplant_Double"]<-"Alloplectus tetragonoides"

#Some reasonable level of presences, 25 points
keep<-names(which(table(int$Hummingbird) > 20))

int<-int[int$Hummingbird %in% keep,]

m.dat<-droplevels(int[colnames(int) %in% c("ID","Video","Time","Hummingbird","Sex","TransectID","Transect_R","Iplant_Double","Pierce","DateP","Month","ele","Type")])

#Does the data come from camera or transect?
m.dat$Type<-(is.na(m.dat$TransectID))*1

m.dat$Year<-years(as.Date(m.dat$DateP))
#one missing date
m.dat$Year[m.dat$Year %in% 2012]<-2013

#Number of bird species
h_species<-nlevels(m.dat$Hummingbird)

#Number of plant species
plant_species<-nlevels(m.dat$Iplant_Double)

#Get english name
dath<-merge(m.dat,hum.morph, by.x="Hummingbird",by.y="English",keep=all)

#Merge to flowers
int.FLlevels<-levels(factor(dath$Iplant_Double))

#Which flowers are we missing info for?
missingTraits<-int.FLlevels[!int.FLlevels %in% fl.morph$X]

#print(paste("Missing Trait Information:",missingTraits))
dath<-merge(dath,fl.morph, by.x="Iplant_Double",by.y="X")

#Drop piercing events, since they don't represent correlation
dath<-dath[!dath$Pierce %in% c("y","Y"),]
```

##Match Species to Morphology


```r
#remove species with less than  10 observations
keep<-names(which(table(dath$Hummingbird) > 20))

dath<-droplevels(dath[dath$Hummingbird %in% keep,])

#observed traitmatching
traitmatchF<-abs(t(sapply(hum.morph$Bill,function(x){x-fl.morph$TotalCorolla})))

rownames(traitmatchF)<-hum.morph$English
colnames(traitmatchF)<-fl.morph$Group.1
```


```r
#match names #Round to 2 decimals #Convert to cm for winbugs, avoids numerical underflow
traitmatchT<-round(traitmatchF[rownames(traitmatchF) %in% dath$Hummingbird,colnames(traitmatchF) %in% dath$Iplant_Double],2)/10

traitmatchT<-traitmatchT[sort(rownames(traitmatchT)),sort(colnames(traitmatchT))]
```

##Elevation ranges

Create a binary variable whether each observation was in a low elevation or high elevation transect. We have some species that just occur at the top of the gradient, and are not present in the sampling window of flowers at the low elevation.

Accounting for non-availability.
We have to figure out which plants were sampled in which periods, and if it was sampled, the non-detection are 0 if it wasn't the non-detection are NA. then remove all the Na's.


```r
elevH<-read.csv("InputData/HummingbirdElevation.csv",row.names=1)
head(elevH)
```

```
##                 Hummingbird  Low        m High Index
## 1            Andean Emerald 1378 1378.632 1380     1
## 2    White-whiskered Hermit 1331 1430.091 1621     1
## 3    Stripe-throated Hermit 1340 1462.730 1558     1
## 4         Crowned Woodnymph 1360 1523.477 2049     1
## 5 Rufous-tailed Hummingbird 1370 1532.000 1862     3
## 6  Wedge-billed Hummingbird 1331 1606.158 1966     3
```

```r
colnames(elevH)[5]<-"Elevation"
elevH$Bird<-1:nrow(elevH)

#high elevation or low elevation
elevP<-read.csv("InputData/PlantElevation.csv",row.names=1)
colnames(elevP)[5]<-"Elevation"
elevP$Plant<-1:nrow(elevP)
elevP$Iplant_Double<-as.character(elevP$Iplant_Double)

#Correct known taxonomic errors
elevP[elevP$Iplant_Double %in% "Alloplectus purpureus","Iplant_Double"]<-"Glossoloma purpureum"
elevP[elevP$Iplant_Double %in% "Capanea affinis","Iplant_Double"]<-"Kohleria affinis"
elevP[elevP$Iplant_Double %in% "Alloplectus teuscheri","Iplant_Double"]<-"Drymonia teuscheri"
elevP[elevP$Iplant_Double %in% "Columnea cinerea","Iplant_Double"]<-"Columnea mastersonii"
elevP[elevP$Iplant_Double %in% "Alloplectus tenuis","Iplant_Double"]<-"Drymonia tenuis"

#Merge to observed Data
#plants
dathp<-merge(dath,elevP,by="Iplant_Double")

#birds
datph<-merge(dathp,elevH,by="Hummingbird")
```

### Summarize Observations


```r
#ID for NA is holger transects, make the id's 1:n for each day of transect at each elevation, assuming no elevation was split across days.
datph$ID<-as.character(datph$ID)

noid<-datph[is.na(datph$ID),]

id_topaste<-paste(noid$Transect_R,noid$DateP,"Transect",sep="_")
datph[which(is.na(datph$ID)),"ID"]<-id_topaste
  
indatraw<- datph %>% group_by(Bird,Plant,ID,DateP) %>% summarize(Yobs=n(),Elev=mean(ele,na.rm=T),Transect_R=unique(Transect_R)) 

indatraw[order(indatraw$Yobs,decreasing=T),]
```

```
## Source: local data frame [562 x 7]
## Groups: Bird, Plant, ID [467]
## 
##     Bird Plant     ID      DateP  Yobs  Elev Transect_R
##    (int) (int)  (chr)     (fctr) (int) (dbl)      (lgl)
## 1     11   101  FH709 2014-02-27    21  1990         NA
## 2     20   109  FH616 2014-05-07    17  2450         NA
## 3     14   101  FH414 2014-04-22    13  1990         NA
## 4     18   123  FL083 2013-07-29    13  2350         NA
## 5      2   120  FL061 2013-07-10    12  1325         NA
## 6     11    28  FH303 2013-11-19    12  1940         NA
## 7     11   101 FH1213 2015-02-11    11  1990         NA
## 8     18   123  FL083 2013-07-28    11  2350         NA
## 9      3    13  FL057 2013-07-05    10  1390         NA
## 10     8   108  NF084 2014-05-24    10  1513         NA
## ..   ...   ...    ...        ...   ...   ...        ...
```

```r
#add unique Camera ID
indatraw$Camera<-as.numeric(factor(indatraw$ID))

#add day ID
sdat<-split(indatraw,list(indatraw$Camera,indatraw$Plant),drop = T)

sdat<-lapply(sdat,function(x){
  x<-droplevels(x)
  x$Day<-as.numeric(as.factor(x$DateP))
  return(x)
})

indatraw<-rbind_all(sdat)

#add months and years
indatraw$Month<-as.numeric(factor(months(strptime(indatraw$DateP,format="%Y-%m-%d")),levels=month.name))
indatraw$Year<-years(indatraw$DateP)

#Species names
for (x in 1:nrow(indatraw)){
  indatraw$Hummingbird[x]<-as.character(elevH[elevH$Bird %in% indatraw$Bird[x],"Hummingbird"])
  indatraw$Iplant_Double[x]<-as.character(elevP[elevP$Plant %in% indatraw$Plant[x],"Iplant_Double"])
}
```

What elevation transect is each observation in?
The camera data need to be inferred from the GPS point.


```r
#cut working best on data.frame
indatraw<-as.data.frame(indatraw)

#which elevation bin is each observation within
labs<-paste(seq(1300,2500,200),seq(1500,2700,200),sep="_")
indatraw$Transect_R[is.na(indatraw$Transect_R)]<-as.character(cut(indatraw[is.na(indatraw$Transect_R),]$Elev,seq(1300,2700,200),labels=labs))
```


```r
#match the traits
traitmelt<-melt(traitmatchT)
colnames(traitmelt)<-c("Hummingbird","Iplant_Double","Traitmatch")
```

##Absences - accounting for non-detection

We have more information than just the presences, given species elevation ranges, we have absences as well. Absences are birds that occur at the elevation of the plant sample, but were not recorded feeding on the flower.


```r
indatlong<-acast(indatraw,Bird~Plant~Camera~Day,value.var="Yobs")
indatlong[is.na(indatlong)]<-0
```


```r
#Only non-detections are real 0's, the rest are NA's and are removed.
#Plants not surveyed in that time period
#Hummingbirds not present at that site

for(x in 1:dim(indatlong)[3]){
  #Remove non sampled plants 
  a<-indatlong[,,x,]
  
  #No observations at that plant
  toNA<-as.numeric(names(which(apply(a,2,sum)==0)))
  pres<-as.numeric(names(which(!apply(a,2,sum)==0)))
  indatlong[,colnames(a) %in% toNA,x,]<-NA
  
  #match elevation if there were presences, slightly terse code.
  if(length(pres)>0){
    for (y in 1:length(pres)){
      #Plant range
      pr<-elevP[elevP$Plant %in% pres[y],"Elevation"]
      #if there is not enough elevation area, make NA
      if(is.na(pr)|!pr==3){
        helim<-elevH[!elevH$Elevation %in% c(pr,3),"Bird"]
        #sanity check, make sure rows don't have values
        sums<-indatlong[rownames(a) %in% helim,colnames(a) %in% pres[y],x,]==0
        helim<-helim[sums]
        #set to NA if they are outside of elevation and are blank
        indatlong[rownames(a) %in% helim ,colnames(a) %in% pres[y],x,]<-NA
        }
    }
  }
}

### There can't be absences in days that weren't sampled.
for (x in 1:dim(indatlong)[3]){
  cam<-indatlong[,,x,]
  for (y in 1:dim(cam)[3]){
    sc<-sum(cam[,,y],na.rm=T)
    if (sc ==0){
      cam[,,y]<-NA
    }
  }
}

#melt and remove Na's
indat<-melt(indatlong)
indat<-indat[!is.na(indat$value),]

colnames(indat)<-c("Bird","Plant","Camera","Day","Yobs")
```


```r
#remerge the time period data
Timelookup<-indatraw %>% dplyr::select(Camera,DateP,Transect_R,Month,Year) %>% group_by(Camera,DateP,Transect_R,Month,Year) %>% distinct() %>% arrange(Camera)

#Get time information
indat<-merge(indat,Timelookup,by=c("Camera"))

#Species names
for (x in 1:nrow(indat)){
  indat$Hummingbird[x]<-as.character(elevH[elevH$Bird %in% indat$Bird[x],"Hummingbird"])
  indat$Iplant_Double[x]<-as.character(elevP[elevP$Plant %in% indat$Plant[x],"Iplant_Double"])
}

#Get trait information
#match the traits
indat<-merge(indat,traitmelt,by=c("Hummingbird","Iplant_Double"))
```

Reformat index for jags.
Jags needs a vector of input species 1:n with no breaks.


```r
indat$Hummingbird<-as.factor(indat$Hummingbird)
indat$Iplant_Double<-as.factor(indat$Iplant_Double)
indat$jBird<-as.numeric(indat$Hummingbird)
indat$jPlant<-as.numeric(indat$Iplant_Double)

jagsIndexBird<-data.frame(Hummingbird=levels(indat$Hummingbird),jBird=1:length(levels(indat$Hummingbird)))
 
jagsIndexPlants<-data.frame(Iplant_Double=levels(indat$Iplant_Double),jPlant=1:length(levels(indat$Iplant_Double)))

#Similiarly, the trait matrix needs to reflect this indexing.
jTraitmatch<-traitmatchT[rownames(traitmatchT) %in% unique(indat$Hummingbird),colnames(traitmatchT) %in% unique(indat$Iplant_Double)]

write.csv(indat,"InputData/ObservedData.csv")
```

# Hierarcichal Occupancy Model

$$ Y_{i,j,k} \sim Binom(N_{i,j,k},detect_i)$$
$$N_{i,j,k} \sim Pois(\lambda_{i,j}) $$
$$log(\lambda_{i,j})<-\alpha_i + \beta_i * abs(Bill_i - Corolla_i) $$
$$detect_i \sim U(0,0.5)$$     

**Priors**

$$\alpha_i \sim N(intercept,\tau_{\alpha})$$
$$\beta_i \sim N(\gamma,\tau_{detect})$$

**Hyperpriors**
$$gamma \sim N(0,0.0001)$$
$$intercept \sim N(0,0.0001)$$

$$\tau_{\alpha} \sim Gamma(0.0001,0.0001)$$
$$\tau_\beta \sim Gamma(0.0001,0.0001)$$
$$\tau_detect \sim Gamma(0.0001,0.0001)$$

**Derived quantities**

$$\sigma_{intercept} = \sqrt[2]{\frac{1}{\tau_\alpha}}$$
$$\sigma_{slope} = \sqrt[2]{\frac{1}{\tau_\beta}}$$
$$\sigma_{detect} = \sqrt[2]{\frac{1}{\tau_detect}}$$

# Poisson GLMM


```r
runs<-50000

#trigger parallel
paralleljags<-T

#Source model
source("Bayesian/NoDetectNmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NoDetectNmixturePoissonRagged.R"))

if(paralleljags){

  #for parallel run
  Yobs=indat$Yobs
  Bird=indat$jBird
  Plant=indat$jPlant
  Plants=max(indat$jPlant)
  Time=indat$Camera
  Traitmatch=jTraitmatch
  #number of birds to iterate
  Birds=max(indat$jBird)
  Nobs=length(indat$Yobs)

  #A blank Y matrix - all present
  Ninit<-rep(max(indat$Yobs)+1,Nobs)
  
  #Inits
  InitStage <- function() {list(beta=rep(0.5,Birds),alpha=rep(0.5,Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 1   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains

  Dat<-list("Yobs","Bird","Plant","Plants","Time","Traitmatch","Birds","Nobs","Ninit")

    m2_niave<-do.call(jags.parallel,list(Dat,InitStage,ParsStage,model.file="Bayesian/NoDetectNmixturePoissonRagged.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc))
    
  } else {
 #Input Data
  Dat <- list(
    Yobs=indat$Yobs,
    Bird=indat$jBird,
    Plant=indat$jPlant,
    Plants=max(indat$jPlant),
    Time=indat$Camera,
    Traitmatch=jTraitmatch,
    #number of birds to iterate
    Birds=max(indat$jBird),
    Nobs=length(indat$Yobs))
  
    #A blank Y matrix - all present
    Ninit<-rep(max(Dat$Yobs)+1,Dat$Nobs)
  
    #Inits
  InitStage <- function() {list(beta=rep(0.5,Dat$Birds),alpha=rep(0.5,Dat$Birds),dtrans=rep(0,Dat$Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 2   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  #Jags
  
  m2_niave <- jags(inits=InitStage,
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
runs<-20000
recompile(m2_niave)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 87109
## 
## Initializing model
## 
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 87109
## 
## Initializing model
```

```r
m2_niave<-update(m2_niave,n.iter=runs,n.burnin=runs*.95)
```


```r
pars_dniave<-extract_par(m2_niave,data=indat,Bird="jBird",Plant="jPlant")
pars_dniave$Model<-"Poisson GLMM"
```

##Assess Convergence


```r
###Chains
ggplot(pars_dniave[pars_dniave$par %in% c("detect","alpha","beta"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-18-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
ggplot(pars_dniave[pars_dniave$par %in% c("gamma","sigma_int","sigma_slope","intercept"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figureObserved/unnamed-chunk-19-1.png" title="" alt="" style="display: block; margin: auto;" />

# Observed Data With Detection


```r
runs<-70000

#trigger parallel
paralleljags<-T

#Source model
source("Bayesian/NmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NmixturePoissonRagged.R"))

if(paralleljags){

  #for parallel run
  Yobs=indat$Yobs
  Bird=indat$jBird
  Plant=indat$jPlant
  Camera=indat$Camera
  Cameras=max(indat$Camera)
  Day<-indat$Day
  Days<-max(indat$Day)
  Traitmatch=jTraitmatch
  #number of birds to iterate
  Birds=max(indat$jBird)
  Plants=max(indat$jPlant)
  Nobs=length(indat$Yobs)

  #A blank Y matrix - all present
  Ninit<-array(dim=c(Birds,Plants,Cameras),data=max(indat$Yobs)+1)

  #Inits
  InitStage <- function() {list(beta=rep(0.5,Birds),alpha=rep(0.5,Birds),dtrans=rep(0,Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew","lambda")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 8   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains

  Dat<-list("Yobs","Bird","Plant","Plants","Traitmatch","Birds","Nobs","Ninit","Day","Days","Camera","Cameras")

    system.time(m2<-do.call(jags.parallel,list(Dat,InitStage,ParsStage,model.file="Bayesian/NmixturePoissonRagged.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc)))
  } else {
  #Input Data
  Dat <- list(
     Yobs=indat$Yobs,
    Bird=indat$jBird,
    Plant=indat$jPlant,
    Camera=indat$Camera,
    Cameras=max(indat$Camera),
    Day<-indat$Day,
    Days<-max(indat$Day),
    Traitmatch=traitmatch,
    #number of birds to iterate
    Birds=max(indat$jBird),
    Plants=max(indat$jPlant),
    Nobs=length(indat$Yobs)
  )
    #A blank Y matrix - all present
  Ninit<-array(dim=c(Dat$Birds,Dat$Plants,Dat$Cameras),data=max(indat$Yobs)+1)
  
    #Inits
  InitStage <- function() {list(beta=rep(0.5,Dat$Birds),alpha=rep(0.5,Dat$Birds),dtrans=rep(0,Dat$Birds),intercept=0,tau_alpha=0.1,tau_beta=0.1,N=Ninit,gamma=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta","intercept","sigma_int","sigma_slope","ynew","gamma","fit","fitnew","lambda")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 5   #thinning rate
  nb <- runs*.95 # number to discard for burn-in
  nc <- 2  # number of chains
  
  #Jags
  
  system.time(m2 <- jags(inits=InitStage,
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
runs<-30000
recompile(m2)
m2<-update(m2,n.iter=runs,n.burnin=runs*.8,n.thin=2)
```


```r
#extract par to data.frame
pars_detect<-extract_par(m2,data=indat,Bird="jBird",Plant="jPlant")

#name
pars_detect$Model<-"Occupancy"
```

##Assess Convergence


```r
###Chains
ggplot(pars_detect[pars_detect$par %in% c("detect","alpha","beta"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-23-1.png" title="" alt="" style="display: block; margin: auto;" />

###Hierarcichal Posteriors


```r
ggplot(pars_detect[pars_detect$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figureObserved/unnamed-chunk-24-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/BothObs.svg",height=5,width=7)
ggsave("Figures/BothObs.jpg",height=5,width=7,dpi=300)
```


```r
#Bind together the two models
parsObs<-rbind(pars_detect,pars_dniave)
```

##Posteriors


```r
###Posterior Distributions
ggplot(parsObs[parsObs$par %in% c("detect","alpha","beta"),],aes(x=estimate,fill=Model)) + geom_histogram(position='identity') + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-26-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
#Detection figure
ggplot(parsObs[parsObs$par %in% c("detect"),],aes(x=as.factor(species),y=estimate,fill=Model)) + geom_violin() + ggtitle("Estimate of parameters") + theme_bw() + ggtitle("Detection Probability") +facet_wrap(~Model,scales="free") 
```

<img src="figureObserved/unnamed-chunk-27-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
pars_detect<-merge(pars_detect,jagsIndexBird,by.x="species",by.y="jBird",all.x=T)

ggplot(pars_detect[pars_detect$par %in% c("detect"),],aes(x=estimate)) + geom_histogram() + ggtitle("Posterior Distribution") + theme_bw() + facet_wrap(~Hummingbird,ncol=5) + xlab("Probability of Detection")
```

<img src="figureObserved/unnamed-chunk-27-2.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/DetectionProb.jpg",dpi=300,height=7,width=11)
```


```r
ggplot(parsObs[parsObs$par %in% c("gamma","intercept","sigma_int","sigma_slope"),],aes(x=estimate,fill=Model)) + geom_histogram() + ggtitle("Trait matching regression parameters") + facet_wrap(~par,scale="free",nrow=2) + theme_bw() 
```

<img src="figureObserved/unnamed-chunk-28-1.png" title="" alt="" style="display: block; margin: auto;" />

##Species Predictions


```r
castdf<-dcast(parsObs[parsObs$par %in% c("beta","alpha"),], species +Chain +Model+ Draw~par,value.var="estimate")

#Turn to species level
castdf$species<-factor(castdf$species,levels=1:max(as.numeric(castdf$species)))

species.split<-split(castdf,list(castdf$species,castdf$Model))

species.traj<-list()

for(d in 1:length(species.split)){
  x<-species.split[[d]]
  #species name
  index<-jagsIndexBird[unique(x$species),"Hummingbird"]
  
  #range of trait distances
  tsp<-indat %>% filter(Hummingbird==index) %>% distinct(Traitmatch) %>% .$Traitmatch
  species.traj[[d]]<-trajF(alpha=x$alpha,beta=x$beta,x=tsp)
  }

names(species.traj)<-names(species.split)

species.traj<-melt(species.traj,id.var=colnames(species.traj[[1]]))

#split out names and model
species.traj[,c("Index","Model")]<-colsplit(species.traj$L1,"\\.",c("Index","Model"))

spe<-merge(species.traj,jagsIndexBird,by.x="Index",by.y="jBird")

#match colnames

#plot and compare to original data
ggplot(data=spe[,],aes(x=x)) + geom_point(data=indat,aes(x=Traitmatch,y=Yobs)) + geom_ribbon(aes(ymin=lower,ymax=upper,fill=Model),alpha=0.2)  + geom_line(aes(y=mean,col=Model),size=1) + theme_bw() + ylab("Interactions") + xlab("Difference between Bill and Corolla Length") + facet_wrap(~Hummingbird,scales="free",ncol=3)+ labs(fill="Model")  + ylab("Interactions per day")
```

<img src="figureObserved/unnamed-chunk-29-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/SpeciesPredictionsBoth.jpg",dpi=300,height=8,width=10)
```

###Overall predicted relationship 

Does accounting for non-independence and detection change our estimate of trait matching?


```r
castdf<-dcast(parsObs[parsObs$par %in% c("gamma","intercept"),], Model+Chain + Draw~par,value.var="estimate")

castdf<-split(castdf,castdf$Model)
#calculated predicted y
predy<-rbind_all(lapply(castdf,function(i){
  #calculate trajectory and append model
  pr<-trajF(alpha=i$intercept,beta=i$gamma,x=indat$Traitmatch)  
  pr$Model<-unique(i$Model)
  return(pr)
  }))

#plot and compare to original data
ggplot(data=predy,aes(x=x,col=Model,fill=Model)) + geom_point(data=indat,aes(x=Traitmatch,y=Yobs),col="black",fill="black",size=2) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.6)  + geom_line(aes(y=mean),size=1) + theme_bw() + ylab("Interactions") + xlab("Difference between Bill and Corolla Length") + ylab("Interactions per day") + ylim(0,13)
```

<img src="figureObserved/unnamed-chunk-30-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/ObsResults.jpg",dpi=300,height=5,width=7)
```

##Discrepancy 

The goodness of fit is a measured as chi-squared. The expected value for each day is the detection rate * the estimate intensity of interactions. The expected value is compared to the observed value of the actual data. In addition, a replicate dataset is generated from the posterior predicted intensity. Better fitting models will have lower discrepancy values and be 
Better fitting models are smaller values and closer to the 1:1 line. A perfect model would be 0 discrepancy. This is unrealsitic given the stochasticity in the sampling processes. Rather, its better to focus on relative discrepancy. In addition, a model with 0 discrepancy would likely be seriously overfit and have little to no predictive power.


```r
fitstat<-parsObs[parsObs$par %in% c("fit","fitnew"),]
fitstat<-dcast(fitstat,Model+Draw+Chain~par,value.var="estimate")

ymin<-round(min(fitstat$fit))
ymax<-round(max(fitstat$fit))
ab<-data.frame(x=0:ymax,y=0:ymax)
disc_obs<-ggplot(fitstat,aes(x=fit,y=fitnew)) + geom_point(aes(col=Model)) + theme_bw() + labs(x="Discrepancy of observed data",y="Discrepancy of replicated data",col="Model")  + ggtitle("Empirical Data") + geom_line(data=ab,aes(x=x,y=y)) + coord_fixed() + ylim(ymin=0,ymax=max(max(c(fitstat$fit,fitstat$fitnew)))) + xlim(xmin=0,xmax=max(max(c(fitstat$fit,fitstat$fitnew))))
disc_obs
```

<img src="figureObserved/unnamed-chunk-31-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/ObservedDiscrepancy.jpeg",width = 5,height=10)
```

##Detection table


```r
dp<-group_by(pars_detect[pars_detect$par %in% c("detect"),],species) %>% summarise(mean=round(mean(estimate,na.rm=T),3)*100,lower=round(quantile(estimate,0.025,na.rm=T),3)*100,upper=round(quantile(estimate,0.975,na.rm=T),3)*100)

tab<-merge(dp,jagsIndexBird,by.x="species",by.y="jBird")[,-1]
tab[,c(4,1,2,3)]
```

```
##                Hummingbird mean lower upper
## 1       Booted Racket-tail 25.5  19.3  32.0
## 2               Brown Inca 25.4  20.7  30.4
## 3      Buff-tailed Coronet 42.3  32.4  49.4
## 4            Collared Inca 19.3  11.8  27.5
## 5        Gorgeted Sunangel 40.2  32.3  47.9
## 6  Green-fronted Lancebill 17.1   7.3  28.1
## 7     Speckled Hummingbird  7.6   2.6  15.0
## 8   Stripe-throated Hermit 34.1  26.9  41.5
## 9     Tawny-bellied Hermit 17.0  13.0  21.3
## 10     Violet-tailed Sylph 18.4  14.5  22.6
## 11  White-whiskered Hermit 28.9  22.6  35.3
```

```r
write.csv(tab[,c(4,1,2,3)],"Figures/Table1.csv")
```

##Sampling intensity and detection for each hummingbird species

The probability of missing a species at each daily visit is 1 - detection probability.

The probability of missing a species on sequential visits is (1- detection probability) * (1 - detection probability).

We are interested in the number of sampling events that minimize this value to some reasonable threshold. I have chosen 0.05 by convention. 

The following figure represesent the estimated number of daily surveys to capture a hummingbird event given that we know it occurs. These data can be thought of as successful draws from a negative binomial distribution.

It is easiest to interpret as the number of days until you are likely to see an interaction, so i prefer to calculate:

$$ 1-(1-p)^n $$
$$ p = probability of detection$$
$$ n = days sample $$


```r
dp<-function(n,p){
  1-((1-p)^n)
}

ts<-split(tab,tab$Hummingbird,drop=T)
detectd<-lapply(ts,function(x){
  meanD<-dp(n=1:10,p=x$mean/100)
  lowerD<-dp(n=1:10,p=x$lower/100)
  upperD<- dp(n=1:10,p=x$upper/100)
  data.frame(Days=1:10,mean=meanD,lower=lowerD,upper=upperD)
})

md<-melt(detectd,id.var="Days")
md<-dcast(md,...~variable)
ggplot(md,aes(x=Days,fill=L1,y=mean,ymin=lower,ymax=upper)) + geom_ribbon(alpha=.5) + geom_line() + facet_wrap(~L1,ncol=3) + geom_hline (yintercept=.5,linetype="dashed") + ylab("Probability of Detection") + scale_fill_discrete(guide="none") + theme_bw() + scale_x_continuous(breaks=seq(0,10,1))
```

<img src="figureObserved/unnamed-chunk-33-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/DetectionDays.jpeg",height=9,width=5,dpi=300) 
```

The number of days it would take to have 50% confidence you have sampled enough to capture known interactions is the x axis value where the dotted line hits the curve.

Number of samples per species

##DIC Table


```r
m2_niave$BUGSoutput$DIC
```

```
## [1] 11835.49
```

```r
m2$BUGSoutput$DIC
```

```
## [1] 10750.06
```

#Predicted versus Observed Data


```r
m<-max(traitmatchT)-traitmatchT

mat<-indat %>% group_by(jBird,jPlant) %>% summarize(n=sum(Yobs))
true_state<-acast(mat,jBird~jPlant,fill=0)

#get species names?

#trait liklihood
dmultinom(true_state,prob=m,log=T)
```

```
## [1] -3298.957
```

```r
paste("Correlation coefficient is:", round(cor(c(true_state),c(m),method="spearman"),2))
```

```
## [1] "Correlation coefficient is: 0.25"
```

###Test Statistic
Chisquared statistic is (observed-expected)^2/(expected + 0.5), we add the 0.5 to avoid dividing by 0. For each combination of birds and plants, predicted versus  observed for each data point.

##Multinomial prediction

Based on trait differences.


```r
#define discrep function
chisq<-function(o,e){(o-e)^2/(e+0.5)}

nullm<-function(){
r<-mgen(m,sum(true_state),keep.species = F)
}

#same number of draws as bayesian
cl<-makeCluster(20,"SOCK")
registerDoSNOW(cl)
mult_matrix<-foreach(x=1:length(pars_detect[pars_detect$par %in% "fit","estimate"]),.packages=c("bipartite","reshape2")) %dopar% nullm()
stopCluster(cl)

mats<-lapply(mult_matrix,function(r){
rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))
colnames(rmerge)<-colnames(true_state)
  rownames(rmerge)<-rownames(true_state)
#for each position what is the chisq
for (x in 1:nrow(r)){
  for (y in 1:ncol(r)){
   rmerge[x,y]<-chisq(o=r[x,y],e=true_state[x,y])
  }
}

#sum of chisq driscrepancy
return(rmerge)
})

multi_disc<-sapply(mats,function(x) sum(x))

qplot(multi_disc)+ xlab("Chisquared Discrepancy for Multimonial Liklihood") + geom_vline(xintercept=mean(multi_disc),col='red',linetype='dashed')
```

<img src="figureObserved/unnamed-chunk-36-1.png" title="" alt="" style="display: block; margin: auto;" />

#Compare Bayesian Occupancy and Multinomial using true known interactions

## Poisson GLMM


```r
N_niave<-pars_dniave[ pars_dniave$par %in% "ynew",]

bydraw<-split(N_niave,list(N_niave$Chain,N_niave$Draw))

occ_nodetect_matrix<-lapply(bydraw,function(x){
  r<-acast(x,species ~ plant,value.var = "estimate",fun.aggregate = sum)
})

#calculate discrep on those deviates
occ_nodetect<-lapply(occ_nodetect_matrix,function(r){
    #for each position what is the chisq
  rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))
  colnames(rmerge)<-colnames(true_state)
  rownames(rmerge)<-rownames(true_state)
  for (x in 1:nrow(r)){
    for (y in 1:ncol(r)){
     rmerge[x,y]<-chisq(o=r[x,y],e=true_state[x,y])
      }
    }
  return(rmerge)
})

names(occ_nodetect)<-1:length(occ_nodetect)
names(occ_nodetect_matrix)<-1:length(occ_nodetect_matrix)
```

##With Detection


```r
N<-pars_detect[pars_detect$par %in% "ynew",]
bydraw<-split(N,list(N$Chain,N$Draw))
occ_matrix<-lapply(bydraw,function(x){
  r<-acast(x,species ~ plant,value.var = "estimate",fun.aggregate = sum)
  })

#calculate discrep for those aggregated matrices
occ<-lapply(occ_matrix,function(r){
    #for each position what is the chisq
  rmerge<-matrix(nrow = nrow(true_state),ncol=ncol(true_state))
  colnames(rmerge)<-colnames(true_state)
  rownames(rmerge)<-rownames(true_state)
  for (x in 1:nrow(r)){
    for (y in 1:ncol(r)){
     rmerge[x,y]<-chisq(o=r[x,y],e=true_state[x,y])
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
simdat<-merge(mocc,mmat,by=c("Bird","Plant"))

#occupancy with nodetection
moccd<-melt(occ_nodetect_matrix)
colnames(moccd)<-c("Bird","Plant","Poisson GLMM","Iteration")

simdat<-merge(simdat,moccd,by=c("Bird","Plant","Iteration"))

#multinomial
multmats<-melt(mult_matrix)
colnames(multmats)<-c("Bird","Plant","Multinomial","Iteration")
simdat<-merge(simdat,multmats,by=c("Bird","Plant","Iteration"))


simdat<-melt(simdat,measure.vars = c("Occupancy","Poisson GLMM","Multinomial"))

ggplot(simdat,aes(x=True_State,y=value,col=variable)) + geom_point() + geom_abline() + labs(col="Model") + ylab("Predicted State") + xlab("True State") + theme_bw() + facet_wrap(~variable)
```

<img src="figureObserved/unnamed-chunk-39-1.png" title="" alt="" style="display: block; margin: auto;" />

##Summary of discrepancy of predicted matrices


```r
#Multinomial
multi_disc<-sapply(mats,function(x) mean(x))

#occupancy without detection
occno_disc<-sapply(occ_nodetect,function(x) mean(x))

#Occupancy with detection
occ_disc<-sapply(occ,function(x) mean(x))

#compared to bayesian
ggplot(data.frame(multi_disc)) + geom_histogram(aes(x=multi_disc),fill="blue",alpha=.6)+ xlab("Chi-squared Discrepancy") + geom_histogram(data=data.frame(occ_disc),aes(x=occ_disc),fill="red",alpha=.6) + theme_bw() +geom_vline(aes(xintercept=mean(occ_disc)),linetype="dashed",col="red")+ geom_vline(xintercept=mean(multi_disc),linetype="dashed",col="blue") + geom_histogram(data=data.frame(occno_disc),aes(x=occno_disc),fill="orange",alpha=.6) + geom_vline(aes(xintercept=mean(occno_disc)),linetype="dashed",col="orange")
```

<img src="figureObserved/unnamed-chunk-40-1.png" title="" alt="" style="display: block; margin: auto;" />

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
## 1 Multinomial     32.33    1.09 14580.59
## 2   Occupancy      4.46    0.79  2010.57
## 3 Poisson_GLM     27.22    2.77 12274.50
```

Merge with morphological data.


```r
jT<-indat %>% group_by(Bird=jBird,Plant=jPlant) %>% summarize(Traitmatch=unique(Traitmatch))
simdat<-merge(simdat,jT,by=c("Bird","Plant"))

mmat<-merge(mmat,jT,by=c("Bird","Plant"))
```

#Predicted total number of visits based on morphology


```r
simT<-simdat %>% group_by(variable,Traitmatch) %>% summarize(Lower=quantile(value,0.05),Upper=quantile(value,0.95),y=mean(value))

ggplot(simT,aes(x=Traitmatch)) + geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=variable),alpha=0.6) + geom_line(aes(y=y,col=variable),linetype='dashed') + theme_bw() + facet_wrap(~variable,nrow=3) + geom_point(data=mmat,aes(y=True_State))
```

<img src="figureObserved/unnamed-chunk-43-1.png" title="" alt="" style="display: block; margin: auto;" />



```r
gc()
```

```
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells   50257384  2684.1   90464440  4831.4   90464440  4831.4
## Vcells 1479361904 11286.7 2782421106 21228.2 2782419736 21228.2
```

```r
save.image("Observed.Rdata")
```



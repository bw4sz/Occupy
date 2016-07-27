# Hierarchical Nmixture Models for species interactions: Empirical Data
Ben Weinstein - Stony Brook University  




```
## [1] "Run Completed at 2016-07-27 11:00:57"
```


```r
#reload if needed
#load("Observed.Rdata")
```

#Load in data


```r
#read in flower morphology data, comes from Nectar.R
droppath<-"C:/Users/Ben/Dropbox/"
fl.morph<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/FlowerMorphology.csv",sep=""))


#use effective corolla where possible.
fl.morph$Corolla<-fl.morph$EffectiveCorolla

fl.morph[is.na(fl.morph$Corolla),"Corolla"]<-fl.morph[is.na(fl.morph$Corolla),"TotalCorolla"]

#fuchsia macrostigma has an undue influence on this analysis, being 3x longer than other flowers, its not clear that birds really have to reach down the full corolla lenghth, use effective corolla length.
fl.morph[fl.morph$Group.1 %in% "Fuchsia macrostigma","Corolla"]<-50

#First row is empty
fl.morph<-fl.morph[-1,]

#Bring in Hummingbird Morphology Dataset, comes from
hum.morph<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/HummingbirdMorphology.csv",sep=""))

#taxonomy change, we are calling them Crowned Woodnymph's now.
hum.morph$English<-as.character(hum.morph$English)

hum.morph$English[hum.morph$English %in% "Green-crowned Woodnymph"]<-"Crowned Woodnymph"

#Bring in Interaction Matrix
int<-read.csv(paste(droppath,"Thesis/Maquipucuna_SantaLucia/Results/Network/HummingbirdInteractions.csv",sep=""),row.names=1)

#one date error
int[int$DateP %in% '2013-07-25',"Month"]<-7

#one duplicate camera error, perhaps two GPS records.
int<-int[!(int$ID %in% "FH1108" & int$Date_F %in% '2014-12-01'),]

#Correct known taxonomic disagreements, atleast compared to traits
int[int$Iplant_Double=="Alloplectus purpureus","Iplant_Double"]<-"Glossoloma purpureum"
int[int$Iplant_Double=="Capanea affinis","Iplant_Double"]<-"Kohleria affinis"
int[int$Iplant_Double=="Columnea cinerea","Iplant_Double"]<-"Columnea mastersonii"
int[int$Iplant_Double=="Alloplectus teuscheri","Iplant_Double"]<-"Drymonia teuscheri"
int[int$Iplant_Double=="Drymonia collegarum","Iplant_Double"]<-"Alloplectus tetragonoides"

#Some reasonable level of presences, 25 points
keep<-names(which(table(int$Hummingbird) > 10))

int<-int[int$Hummingbird %in% keep & !int$Hummingbird %in% c("Sparkling Violetear"),]

m.dat<-droplevels(int[colnames(int) %in% c("ID","Video","Time","Hummingbird","Sex","TransectID","Transect_R","Iplant_Double","Pierce","DateP","Month","ele","Type")])

#Does the data come from camera or transect?
m.dat$Type<-(is.na(m.dat$TransectID))*1

m.dat$Year<-years(as.Date(m.dat$DateP))
#one missing date
m.dat$Year[m.dat$Year %in% 2012]<-2013
m.dat$Year[m.dat$Year %in% 2106]<-2016

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
#observed traitmatching
traitmatchF<-abs(t(sapply(hum.morph$Bill,function(x){x-fl.morph$Corolla})))
rownames(traitmatchF)<-hum.morph$English
colnames(traitmatchF)<-fl.morph$Group.1
```


```r
#match names #Round to 2 decimals #Convert to cm for winbugs, avoids numerical underflow
traitmatchT<-round(traitmatchF[rownames(traitmatchF) %in% dath$Hummingbird,colnames(traitmatchF) %in% dath$Iplant_Double],2)
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
## 2         Crowned Woodnymph 1367 1383.190 1380     1
## 3    White-whiskered Hermit 1360 1399.436 1400     1
## 4    Stripe-throated Hermit 1370 1408.319 1400     1
## 5 Rufous-tailed Hummingbird 1370 1435.700 1380     1
## 6  Wedge-billed Hummingbird 1331 1658.059 1950     3
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

What elevation transect is each observation in?
The camera data need to be inferred from the GPS point.


```r
#cut working best on data.frame
datph<-as.data.frame(datph)

#which elevation bin is each observation within
labs<-paste(seq(1300,2500,200),seq(1500,2700,200),sep="_")

#for the couple points that have 1290 elevation, round up to 300 for convienance
datph$ele[datph$ele < 1300]<-1301
#make sure transect is a character
datph$Transect_R<-as.character(datph$Transect_R)
datph$Transect_R[is.na(datph$Transect_R)]<-as.character(cut(datph[is.na(datph$Transect_R),]$ele,seq(1300,2700,200),labels=labs))

#Elev for the transects is the midpoint
tran_elev<-datph[datph$Survey_Type=='Transect',"Transect_R"]
datph[datph$Survey_Type=='Transect',"ele"]<-sapply(tran_elev,function(x){
  mean(as.numeric(str_split(x,"_")[[1]]))
})
```

### Define Time Events


```r
#ID for NA is holger transects, make the id's 1:n for each day of transect at each elevation, assuming no elevation was split across days.
datph$ID<-as.character(datph$ID)
noid<-datph[is.na(datph$ID),]

id_topaste<-paste(noid$Month,noid$Year,"Transect",sep="_")
datph[which(is.na(datph$ID)),"ID"]<-id_topaste

#Create year month combination
datph$Time<-paste(datph$Month,datph$Year,sep="_")

#Label survey type
datph$Survey_Type<-NA

mt<-!is.na(datph$TransectID)*1
datph$Survey_Type[mt==1]<-"Transect"
datph$Survey_Type[!datph$Survey_Type %in% "Transect"]<-"Camera"

datph<-datph[datph$Survey_Type=="Camera",]

#Day level
#add day ID
sdat<-split(datph,list(datph$ID),drop = T)

sdat<-lapply(sdat,function(x){
  x<-droplevels(x)
  x$Day<-as.numeric(as.factor(x$DateP))
  return(x)
})

indatraw<-rbind_all(sdat)

#Species names
for (x in 1:nrow(indatraw)){
  indatraw$Hummingbird[x]<-as.character(elevH[elevH$Bird %in% indatraw$Bird[x],"Hummingbird"])
  indatraw$Iplant_Double[x]<-as.character(elevP[elevP$Plant %in% indatraw$Plant[x],"Iplant_Double"])
}

#match the traits
traitmelt<-melt(traitmatchT)
colnames(traitmelt)<-c("Hummingbird","Iplant_Double","Traitmatch")

#dummy presence variable
indatraw$Yobs<-1

#prune columsn to make more readable
indatraw<-indatraw[,c("Hummingbird","Iplant_Double","ID","Time","Month","Year","Transect_R","ele","DateP","Yobs","Day","Survey_Type","Pierce")]
```

##Summarize daily interactions
To estimate the daily detectability, there can only be a max of one interaction per day.
We use mean elevation to average across observations within a transect

```r
indatraw<-indatraw %>% group_by(Hummingbird,Iplant_Double,ID,Day) %>% summarize(Yobs=sum(Yobs),Time=unique(Time),Transect_R=unique(Transect_R),Month=unique(Month),Year=unique(Year),ele=mean(ele),DateP=unique(DateP),Survey_Type=unique(Survey_Type)) %>% ungroup()
```

##Absences - accounting for non-detection

We have more information than just the presences, given species elevation ranges, we have absences as well. Absences are birds that occur at the elevation of the plant sample, but were not recorded feeding on the flower.


```r
#Only non-detections are real 0's, the rest are NA's and are removed.
#Plants not surveyed in that time period
#Hummingbirds not present at that elevation

#For each ID
Time<-unique(indatraw$Time)

#absences data frame
absences<-list()

for(t in Time){
  IDlist<-unlist(unique(indatraw[indatraw$Time ==t,"ID"]))

  for (j in IDlist){
  #Which plants were sampled
  a<-indatraw %>% filter(Time==t,ID==j)
  
  #For each sampled transect
  trans<-unique(a$Transect_R)
  
  if(!length(trans)==0){
    for(transect in trans){

    #for each date 
    datec<-a %>% filter(Transect_R %in% transect)
    datecam<-unique(datec$DateP)
    }} else{
      datecam<-a %>% distinct(DateP) %>% .$DateP
    }
    for(Date in datecam){
      
    #for each plant along that transect at that date
    pres<-a %>% filter(DateP %in% Date) %>% distinct(Iplant_Double) %>% .$Iplant_Double
    
    #Which day in sampling
    dday<-a %>% filter(Transect_R %in% transect,DateP %in% Date) %>% distinct(Day) %>% .$Day

      for (plant in pres){
        #Get mean elevation of that plant record
        camelev<- a %>% filter(Transect_R %in% transect,DateP %in% Date,Iplant_Double %in% plant) %>% .$ele %>% mean()
        
        #Which birds are present at that observation
        predh<-elevH[((elevH$Low < camelev) & (camelev < elevH$High)),"Hummingbird"]
        
        #remove the ones seen on that plant
        hum_present<-a %>% filter(Transect_R %in% transect,DateP %in% Date,Iplant_Double %in% plant) %>% .$Hummingbird
        abbh<-predh[!predh %in% hum_present]
        if(length(abbh)==0){next}
        
        #Make absences from those )(cat not the best)
        add_absences<-data.frame(Hummingbird=abbh,Iplant_Double=plant,Time=t,ID=j,DateP=Date,Month=min(a$Month),Year=unique(a$Year),Transect_R=transect,ele=camelev,Day=unique(dday),Survey_Type=unique(a$Survey_Type),Yobs=0)
        absences<-append(absences,list(add_absences))
      }
    }
  }
}
    
indatab<-rbind_all(absences)

#merge with original data
indat<-rbind_all(list(indatraw,indatab))
```


```r
#Get trait information
#match the traits
indat<-merge(indat,traitmelt,by=c("Hummingbird","Iplant_Double"))
```

#Resources at each month

In our model the covariate is indexed at the scale at which the latent count is considered fixed. This means we need the resource availability per month across the entire elevation gradient for each point.


```r
#Get flower transect data
full.fl<-read.csv("C:/Users/Ben/Dropbox/Thesis/Maquipucuna_SantaLucia/Results/FlowerTransects/FlowerTransectClean.csv")[,-1]

 #month should be capital 
colnames(full.fl)[colnames(full.fl) %in% "month"]<-"Month"

#group by month and replicate, remove date errors by making a max of 10 flowers, couple times where the gps places it in wrong transect by 1 to 2 meters. 
flower.month<-group_by(full.fl,Month,Year,Transect_R,Date_F) %>% dplyr::summarise(Flowers=sum(Total_Flowers,na.rm=TRUE))  %>% filter(Flowers>20)
  
#Make month abbreviation column, with the right order
flower.month$Month.a<-factor(month.abb[flower.month$Month],month.abb[c(1:12)])

#Make year factor column
flower.month$Year<-as.factor(flower.month$Year)

#get quantile for each transect
#thresh<-melt(group_by(flower.month) %>% summarize(Threshold=quantile(Flowers,0.5)))
flower.month$R<-cut(flower.month$Flowers,breaks=c(0,quantile(flower.month$Flowers,0.33),quantile(flower.month$Flowers,0.66),max(flower.month$Flowers)),label=c("Low","Medium","High"))

#fix the levels
flower.month$PTransect_R<-flower.month$Transect_R
levels(flower.month$PTransect_R)<-c("1300m - 1500m", "1500m - 1700m","1700m - 1900m","1900m - 2100m","2100m - 2300m","2300m - 2500m")
#plot

ggplot(flower.month,aes(x=Month.a,log(Flowers),col=R,shape=as.factor(Year))) + geom_point(size=3) + theme_bw()  + geom_smooth(aes(group=1)) + ylab("Flowers") + xlab("Month") + facet_wrap(~PTransect_R) + labs(shape="Year", y= "Log Available Flowers") + scale_x_discrete(breaks=month.abb[seq(1,12,2)]) + scale_color_manual(labels=c("Low","Medium","High"),values=c("black","blue","red")) + labs(col="Resource Availability")
```

<img src="figureObserved/unnamed-chunk-13-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/FlowerMonth.jpeg",dpi=600,height=5,width=9)

#turn min and max elvation into seperate columns for the range
flower.month$minElev<-as.numeric(str_extract(flower.month$Transect_R,"(\\d+)"))
flower.month$maxElev<-as.numeric(str_match(flower.month$Transect_R,"(\\d+)_(\\d+)")[,3])
```


```r
indat$All_Flowers<-NA
indat$Used_Flowers<-NA
indat$FlowerA<-NA

#Resource list for each species.
slist<-int %>% group_by(Hummingbird,Iplant_Double) %>% distinct() %>% dplyr::select(Hummingbird,Iplant_Double) %>% arrange(Hummingbird)

#Create time ID for flower transects
full.fl$Time<-paste(full.fl$Month,full.fl$Year,sep="_")

#all flowers for each ID period
allF<-full.fl %>% group_by(Month,Year,Transect_R,Date_F) %>% summarize(n=sum(Total_Flowers,na.rm=T)) %>% summarize(mn=mean(n)) %>% summarize(F=sum(mn)) %>% as.data.frame()

#Individual flowers for each ID period
indF<-full.fl %>% group_by(Iplant_Double,Month,Year,Transect_R,Date_F) %>% summarize(n=sum(Total_Flowers,na.rm=T)) %>% summarize(mn=mean(n)) %>% summarize(F=sum(mn)) %>% as.data.frame()

for (x in 1:nrow(indat)){

#All flowers
 indat$All_Flowers[x]<-allF[allF$Month %in% indat$Month[x] & allF$Year %in% indat$Year[x],"F"]
 
 #filter by species used by hummingbird
 sp_list<-slist[slist$Hummingbird %in% indat$Hummingbird[x],"Iplant_Double"]

 indat$Used_Flowers[x]<-sum(indF[indF$Iplant_Double %in% sp_list$Iplant_Double & indF$Month %in% indat$Month[x] & indF$Year %in% indat$Year[x],"F"])
  
  #just the abundance of that species
  indat$FlowerA[x]<-sum(indF[indF$Iplant_Double %in% indat$Iplant_Double[x] & indF$Month %in% indat$Month[x] & indF$Year %in% indat$Year[x],"F"])

}
```

###Relationship between resource measures


```r
ggplot(indat,aes(x=All_Flowers,y=Used_Flowers)) + geom_point() + facet_wrap(~Hummingbird,scales="free")
```

<img src="figureObserved/unnamed-chunk-15-1.png" title="" alt="" style="display: block; margin: auto;" />

##Binary Measures of Resources


```r
#All Resources
#indat$BAll_Flowers<-(indat$Month  %in% c("6","7","8","9","10"))*1

indat$BAll_Flowers<-(indat$All_Flowers > quantile(indat$All_Flowers,0.5))*1

qthresh<-indat %>% group_by(Hummingbird) %>% summarize(UThresh=quantile(Used_Flowers,0.75))

indat<-merge(indat,qthresh)
indat$BUsed_Flowers<-(indat$Used_Flowers > indat$UThresh)*1

fthresh<-indat %>% group_by(Hummingbird) %>% summarize(FThresh=mean(FlowerA))
indat<-merge(indat,fthresh)
indat$BFlowerA<-(indat$FlowerA > indat$FThresh)*1

#merge with flower month, split by elevation, mean per month
sflowers<-flower.month %>% group_by(Transect_R,Month,Year) %>% summarize(Flowers=mean(Flowers))
sflowers$R<-cut(sflowers$Flowers,breaks=c(0,quantile(sflowers$Flowers,0.33),quantile(sflowers$Flowers,0.66),max(sflowers$Flowers)),label=c("Low","Medium","High"))
 
indat<-merge(indat,sflowers,c("Transect_R","Month","Year"))
```


```r
#Combine resources with observed data
f<-(indat$Survey_Type=="Camera")*1
f[f==0]<-NA
indat$Camera<-indat$Yobs * f

f<-(indat$Survey_Type=="Transect")*1
f[f==0]<-NA
indat$Transect<-indat$Yobs * f
```


Reformat index for jags.
Jags needs a vector of input species 1:n with no breaks.


```r
#Easiest to work with jags as numeric ordinal values
indat$Hummingbird<-as.factor(indat$Hummingbird)
indat$Iplant_Double<-as.factor(indat$Iplant_Double)
indat$jBird<-as.numeric(indat$Hummingbird)
indat$jPlant<-as.numeric(indat$Iplant_Double)

jagsIndexBird<-data.frame(Hummingbird=levels(indat$Hummingbird),jBird=1:length(levels(indat$Hummingbird)))
 
jagsIndexPlants<-data.frame(Iplant_Double=levels(indat$Iplant_Double),jPlant=1:length(levels(indat$Iplant_Double)))

#Similiarly, the trait matrix needs to reflect this indexing.
jTraitmatch<-traitmatchT[rownames(traitmatchT) %in% unique(indat$Hummingbird),colnames(traitmatchT) %in% unique(indat$Iplant_Double)]
```


```r
indat<-droplevels(indat)

#Turn Time and ID into numeric indexes
indat$jTime<-as.numeric(as.factor(indat$Time))
indat$jID<-as.numeric(as.factor(indat$ID))

#index resources
indat$scaledR<-(indat$FlowerA>0)*1
resourcemat<-indat %>% group_by(jBird,jPlant,jID) %>% summarize(v=max(scaledR))  %>% acast(jBird ~ jPlant ~ jID,value.var='v',fill=0)
```

# Hierarcichal Nmixture Model

For hummingbird i visiting plant j recorded by camera k on day d:

$$ Y_{i,j,k,d} \sim Binom(N_{i,j,k},\omega_i)$$
$$N_{i,j,k} \sim Pois(\lambda_{i,j,k} * \Resource_{i,j,k} $$
$$log(\lambda_{i,j})<-\alpha_i + \beta_{1,i} * |Bill_i - Corolla_j|$$ 


**Priors**

Please recall that jags parameterizes models using precision, not sd (precision = 1/sd^2)

$$\omega_i \sim (\mu_{\omega},\tau_{\omega})$$  $$\mu_{\omega} \sim Normal(0,0.5)   
$$\tau_{\omega} \sim Uniform(0,10)

$$\alpha_i \sim Normal(\mu_{\alpha},\tau_{\alpha})$$
$$\beta_{1,i} \sim Normal(\mu_{\beta_1},\tau_{\beta_1})$$
$$\beta_{2,i} \sim Normal(\mu_{\beta_2},\tau_{\beta_2})$$

**Hyperpriors**
$$\mu_{\alpha} \sim Normal(0,0.0001)$$
$$\mu_{\beta_1} \sim Normal(0,0.0001)$$

$$\tau_{\alpha} \sim Half-T(0.0001,0.0001)$$
$$\tau_{\beta_1} = \sqrt[2]{\frac{1}{\sigma_{\beta_1}}}$$
$$\tau_{\beta_2} = \sqrt[2]{\frac{1}{\sigma_{\beta_2}}}$$

$$\sigma_{\alpha} = \sqrt[2]{\frac{1}{\tau_\alpha}}$$
$$\sigma_{\beta_1} \sim Half-T(0,1)$$
$$\sigma_{\beta_2} \sim Half-T(0,1)$$

# Poisson GLMM


```r
runs<-5000

#Source model
source("Bayesian/NoDetectNmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NoDetectNmixturePoissonRagged.R"))
```

```
##  [1]                                                                                          
##  [2] sink("Bayesian/NoDetectNmixturePoissonRagged.jags")                                      
##  [3]                                                                                          
##  [4] cat("                                                                                    
##  [5]     model {                                                                              
##  [6]     #Compute intensity for each pair of birds and plants                                 
##  [7]     for (i in 1:Birds){                                                                  
##  [8]     for (j in 1:Plants){                                                                 
##  [9]     for (k in 1:Cameras){                                                                
## [10]                                                                                          
## [11]     #Process Model                                                                       
## [12]     log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]                            
## [13]     }                                                                                    
## [14]     }                                                                                    
## [15]     }                                                                                    
## [16]                                                                                          
## [17]     for (x in 1:Nobs){                                                                   
## [18]                                                                                          
## [19]       #expected intensity                                                                
## [20]       N[x]<-lambda[Bird[x],Plant[x],Camera[x]] * resources[Bird[x],Plant[x],Camera[x]]   
## [21]                                                                                          
## [22]       # Observed State                                                                   
## [23]       Yobs[x] ~ dpois(N[x]+0.00000001)                                                   
## [24]                                                                                          
## [25]       #Assess Model Fit                                                                  
## [26]                                                                                          
## [27]       #Fit discrepancy statistics                                                        
## [28]       eval[x]<-lambda[Bird[x],Plant[x],Camera[x]] * resources[Bird[x],Plant[x],Camera[x]]
## [29]       E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)                                       
## [30]                                                                                          
## [31]       ynew[x]~dpois(lambda[Bird[x],Plant[x],Camera[x]])                                  
## [32]       E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)                                   
## [33]                                                                                          
## [34]       }                                                                                  
## [35]                                                                                          
## [36]       for (i in 1:Birds){                                                                
## [37]       alpha[i] ~ dnorm(alpha_mu,alpha_tau)                                               
## [38]       beta1[i] ~ dnorm(beta1_mu,beta1_tau)                                               
## [39]       }                                                                                  
## [40]                                                                                          
## [41]     #Hyperpriors                                                                         
## [42]                                                                                          
## [43]     #Intercept grouping                                                                  
## [44]     alpha_mu~dnorm(0,0.0001)                                                             
## [45]                                                                                          
## [46]     # Group intercept variance                                                           
## [47]     alpha_sigma ~ dt(0,1,1)I(0,)                                                         
## [48]     alpha_tau <- pow(alpha_sigma,-2)                                                     
## [49]                                                                                          
## [50]     #Trait Slope                                                                         
## [51]     #Mean                                                                                
## [52]     beta1_mu~dnorm(0,0.0001)                                                             
## [53]                                                                                          
## [54]     #Variance                                                                            
## [55]     beta1_sigma ~ dt(0,1,1)I(0,)                                                         
## [56]     beta1_tau <- pow(beta1_sigma,-2)                                                     
## [57]                                                                                          
## [58]     #derived posterior check                                                             
## [59]     fit<-sum(E[]) #Discrepancy for the observed data                                     
## [60]     fitnew<-sum(E.new[])                                                                 
## [61]                                                                                          
## [62]     }                                                                                    
## [63]     ",fill=TRUE)                                                                         
## [64]                                                                                          
## [65] sink()
```

```r
  #Data objects for parallel run
  Yobs=indat$Yobs
  Bird=indat$jBird
  Birds=max(indat$jBird)
  Plant=indat$jPlant
  Plants=max(indat$jPlant)
  Camera=indat$jID
  Cameras=max(indat$jID)
  Traitmatch=jTraitmatch
  Nobs=length(indat$Yobs)
  resources=resourcemat

  #A blank Y matrix - all present
  Ninit<-array(dim=c(Nobs),1)
  
  #Inits
  InitStage <- function() {list(beta1=rep(0.5,Birds),alpha=rep(0.5,Birds),alpha_mu=0,beta1_mu=0,beta2_mu=0)}
  
  #Parameters to track
  ParsStage <- c("alpha","beta1","alpha_mu","beta1_mu","ynew","fit","fitnew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 4   #thinning rate
  nb <- runs-2000 # number to discard for burn-in
  nc <- 2  # number of chains

  Dat<-list("Yobs","Bird","Plant","Plants","Camera","Cameras","Traitmatch","Birds","Ninit","Nobs","resources","nb","nt","nc","ni")
    
  system.time(m2_niave<-jags.parallel(Dat,InitStage,ParsStage,model.file="Bayesian/NoDetectNmixturePoissonRagged.jags", n.iter=ni,n.burnin=nb,n.chains=nc,n.thin=nt))
```

```
##    user  system elapsed 
##    7.25    0.80  152.78
```


```r
#recompile if needed
load.module("dic")
runs<-100000
recompile(m2_niave)
m2_niave<-update(m2_niave,n.iter=runs,n.burnin=runs*.9,n.thin = 5)
```


```r
pars_dniave<-extract_par(m2_niave,data=indat,Bird="jBird",Plant="jPlant")
pars_dniave$Model<-"Poisson GLMM"
```

##Assess Convergence


```r
###Chains
ggplot(pars_dniave[pars_dniave$par %in% c("alpha","beta1","beta2"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-23-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
ggplot(pars_dniave[pars_dniave$par %in% c("beta1_mu","beta2_mu","sigma_alpha","beta1_sigma","sigma_beta2","alpha_mu"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figureObserved/unnamed-chunk-24-1.png" title="" alt="" style="display: block; margin: auto;" />

# Observed Data With Detection

## Traits


```r
runs<-5000

#Source model
source("Bayesian/NmixturePoissonRagged_traits.R")

#print model
print.noquote(readLines("Bayesian//NmixturePoissonRagged_traits.R"))
```

```
##  [1]                                                                      
##  [2] sink("Bayesian/NmixturePoissonRagged_traits.jags")                   
##  [3]                                                                      
##  [4] cat("                                                                
##  [5]     model {                                                          
##  [6]     #Compute intensity for each pair of birds and plants             
##  [7]     for (i in 1:Birds){                                              
##  [8]     for (j in 1:Plants){                                             
##  [9]     for (k in 1:Cameras){                                            
## [10]                                                                      
## [11]     #Process Model                                                   
## [12]     log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]        
## [13]                                                                      
## [14]     #For each camera - there is a latent count                       
## [15]     N[i,j,k] ~ dpois(lambda[i,j,k])                                  
## [16]     }                                                                
## [17]     }                                                                
## [18]     }                                                                
## [19]                                                                      
## [20]     #Observed counts for each day of sampling at that camera         
## [21]     for (x in 1:Nobs){                                               
## [22]                                                                      
## [23]     #Observation Process                                             
## [24]     Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])    
## [25]                                                                      
## [26]     #Assess Model Fit                                                
## [27]                                                                      
## [28]     #Fit discrepancy statistics                                      
## [29]     eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Camera[x]]           
## [30]     E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)                     
## [31]                                                                      
## [32]     ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])      
## [33]     E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)                 
## [34]                                                                      
## [35]     }                                                                
## [36]                                                                      
## [37]     for (i in 1:Birds){                                              
## [38]     logit(detect[i]) <- dtrans[i]                                    
## [39]     dtrans[i] ~ dnorm(dprior,tau_detect)                             
## [40]     alpha[i] ~ dnorm(alpha_mu,alpha_tau)                             
## [41]     beta1[i] ~ dnorm(beta1_mu,beta1_tau)                             
## [42]     }                                                                
## [43]                                                                      
## [44]     #Hyperpriors                                                     
## [45]                                                                      
## [46]     #Intercept grouping                                              
## [47]     alpha_mu~dnorm(0,0.0001)                                         
## [48]                                                                      
## [49]     # Group intercept variance                                       
## [50]     alpha_tau ~ dgamma(0.0001,0.0001)                                
## [51]     alpha_sigma<-pow(1/alpha_tau,0.5)                                
## [52]                                                                      
## [53]     #Detect grouping                                                 
## [54]     dprior ~ dnorm(0,0.5)                                            
## [55]                                                                      
## [56]     # Detect variance                                                
## [57]     tau_detect ~ dunif(0,5)                                          
## [58]     sigma_detect<-pow(1/tau_detect,0.5)                              
## [59]                                                                      
## [60]     #Trait Slope                                                     
## [61]     beta1_mu~dnorm(0,0.0001)                                         
## [62]                                                                      
## [63]     #Slope variance, turning precision to sd                         
## [64]     beta1_tau ~ dgamma(0.0001,0.0001)                                
## [65]     beta1_sigma<-pow(1/beta1_tau,0.5)                                
## [66]                                                                      
## [67]     #derived posterior check                                         
## [68]     fit<-sum(E[]) #Discrepancy for the observed data                 
## [69]     fitnew<-sum(E.new[])                                             
## [70]                                                                      
## [71]     }                                                                
## [72]     ",fill=TRUE)                                                     
## [73]                                                                      
## [74] sink()
```

```r
  #for parallel run
  Yobs=indat$Yobs
  Bird=indat$jBird
  Plant=indat$jPlant
  Camera=indat$jID
  Cameras=max(indat$jID)
  Traitmatch=jTraitmatch
  Birds=max(indat$jBird)
  Plants=max(indat$jPlant)
  Nobs=length(indat$Yobs)

  #A blank Y matrix - all present
  Ninit<-array(dim=c(Birds,Plants,Cameras),data=max(indat$Yobs)+1)

  #Inits
  InitStage <- function() {list(beta1=rep(0.5,Birds),alpha=rep(0.5,Birds),N=Ninit,beta1_mu=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta1","alpha_mu","sigma_alpha","beta1_sigma","beta1_mu","dprior")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 4   #thinning rate
  nb <- runs-2000 # number to discard for burn-in
  nc <- 2  # number of chains

  Dat<-list("Yobs","Bird","Plant","Plants","Traitmatch","Birds","Nobs","Ninit","Camera","Cameras","nb","nc","ni","nt")

    system.time(traits<-jags.parallel(Dat,InitStage,ParsStage,model.file="Bayesian/NmixturePoissonRagged_traits.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc))
```

```
##    user  system elapsed 
##    0.26    0.68  751.77
```


```r
#recompile if needed
load.module("dic")
runs<-100000
recompile(traits)
traits<-update(traits,n.iter=runs,n.burnin=runs*.9,n.thin=5)
```


```r
#extract par to data.frame
pars_detect_traits<-extract_par(traits,data=indat,Bird="jBird",Plant="jPlant",ynew=F)

#name
pars_detect_traits$Model<-"Nmixture"
```

###Assess Convergence


```r
###Chains
ggplot(pars_detect_traits[pars_detect_traits$par %in% c("detect","alpha","beta1"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-28-1.png" title="" alt="" style="display: block; margin: auto;" />

###Hierarcichal Posteriors


```r
ggplot(pars_detect_traits[pars_detect_traits$par %in% c("beta1_mu","alpha_mu","sigma_alpha","beta1_sigma","dprior"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figureObserved/unnamed-chunk-29-1.png" title="" alt="" style="display: block; margin: auto;" />

## Traits + Abundance


```r
runs<-5000

#Source model
source("Bayesian/NmixturePoissonRagged.R")

#print model
print.noquote(readLines("Bayesian//NmixturePoissonRagged.R"))
```

```
##  [1]                                                                                                   
##  [2] sink("Bayesian/NmixturePoissonRagged.jags")                                                       
##  [3]                                                                                                   
##  [4] cat("                                                                                             
##  [5]     model {                                                                                       
##  [6]     #Compute intensity for each pair of birds and plants                                          
##  [7]     for (i in 1:Birds){                                                                           
##  [8]     for (j in 1:Plants){                                                                          
##  [9]     for (k in 1:Cameras){                                                                         
## [10]                                                                                                   
## [11]     #Process Model                                                                                
## [12]     log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j]                                     
## [13]                                                                                                   
## [14]                                                                                                   
## [15]     #For each camera - there is a latent count                                                    
## [16]     N[i,j,k] ~ dpois(lambda[i,j,k] * resources[i,j,k] + 0.0000001)                                
## [17]     }                                                                                             
## [18]     }                                                                                             
## [19]     }                                                                                             
## [20]                                                                                                   
## [21]                                                                                                   
## [22]     #Observed counts for each day of sampling at that camera                                      
## [23]     for (x in 1:Nobs){                                                                            
## [24]                                                                                                   
## [25]     #Observation Process                                                                          
## [26]     Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])                                 
## [27]                                                                                                   
## [28]     #Assess Model Fit                                                                             
## [29]                                                                                                   
## [30]     #Fit discrepancy statistics                                                                   
## [31]     eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Camera[x]] * resources[Bird[x],Plant[x],Camera[x]]
## [32]     E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)                                                  
## [33]                                                                                                   
## [34]     ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])                                   
## [35]     E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)                                              
## [36]                                                                                                   
## [37]     }                                                                                             
## [38]                                                                                                   
## [39]     for (i in 1:Birds){                                                                           
## [40]     logit(detect[i]) <- dtrans[i]                                                                 
## [41]     dtrans[i] ~ dnorm(dprior,tau_detect)                                                          
## [42]     alpha[i] ~ dnorm(alpha_mu,alpha_tau)                                                          
## [43]     beta1[i] ~ dnorm(beta1_mu,beta1_tau)                                                          
## [44]     }                                                                                             
## [45]                                                                                                   
## [46]     #Hyperpriors                                                                                  
## [47]                                                                                                   
## [48]     #Intercept grouping                                                                           
## [49]     alpha_mu~dnorm(0,0.0001)                                                                      
## [50]                                                                                                   
## [51]     #Group intercept variance                                                                     
## [52]     alpha_sigma ~ dt(0,1,1)I(0,)                                                                  
## [53]     alpha_tau <- pow(alpha_sigma,-2)                                                              
## [54]                                                                                                   
## [55]     #Detect grouping                                                                              
## [56]     dprior ~ dnorm(0,0.5)                                                                         
## [57]                                                                                                   
## [58]     #Detect variance                                                                              
## [59]     tau_detect ~ dunif(0,5)                                                                       
## [60]     sigma_detect<-pow(1/tau_detect,0.5)                                                           
## [61]                                                                                                   
## [62]     #Trait Slope                                                                                  
## [63]                                                                                                   
## [64]     #Mean                                                                                         
## [65]     beta1_mu~dnorm(0,0.0001)                                                                      
## [66]                                                                                                   
## [67]     #Variance                                                                                     
## [68]     beta1_sigma ~ dt(0,1,1)I(0,)                                                                  
## [69]     beta1_tau <- pow(beta1_sigma,-2)                                                              
## [70]                                                                                                   
## [71]     #derived posterior check                                                                      
## [72]                                                                                                   
## [73]     fit<-sum(E[]) #Discrepancy for the observed data                                              
## [74]     fitnew<-sum(E.new[])                                                                          
## [75]                                                                                                   
## [76]     }                                                                                             
## [77]     ",fill=TRUE)                                                                                  
## [78]                                                                                                   
## [79] sink()
```

```r
  #for parallel run
  Yobs=indat$Yobs
  Bird=indat$jBird
  Plant=indat$jPlant
  Camera=indat$jID
  Cameras=max(indat$jID)
  Traitmatch=jTraitmatch
  Birds=max(indat$jBird)
  Plants=max(indat$jPlant)
  Nobs=length(indat$Yobs)
  resources=resourcemat

  #A blank Y matrix - all present
  Ninit<-array(dim=c(Birds,Plants,Cameras),data=max(indat$Yobs)+1)

  #Inits
  InitStage <- function() {list(beta1=rep(0.5,Birds),alpha=rep(0.5,Birds),alpha_mu=0,sigma_alpha=0.1,sigma_beta2=0.1,beta1_sigma=0.1,N=Ninit,beta1_mu=0)}
  
  #Parameters to track
  ParsStage <- c("detect","alpha","beta1","beta2","alpha_mu","sigma_alpha","beta1_sigma","sigma_beta2","beta1_mu","beta2_mu","fit","fitnew","dprior","sigma_detect","ynew")
  
  #MCMC options
  ni <- runs  # number of draws from the posterior
  nt <- 4   #thinning rate
  nb <- runs-2000 # number to discard for burn-in
  nc <- 2  # number of chains

  Dat<-list("Yobs","Bird","Plant","Plants","Traitmatch","Birds","Nobs","Ninit","Camera","Cameras","resources","nc","nb","ni","nt")

    system.time(m2<-jags.parallel(Dat,InitStage,parameters.to.save=ParsStage,model.file="Bayesian/NmixturePoissonRagged.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc))
```

```
##    user  system elapsed 
##    4.20    2.37  856.86
```


```r
#recompile if needed
load.module("dic")
runs<-100000
recompile(m2)
m2<-update(m2,n.iter=runs,n.burnin=runs*.8,n.thin=5,parameters.to.save=ParsStage)
```


```r
#extract par to data.frame
pars_detect<-extract_par(m2,data=indat,Bird="jBird",Plant="jPlant")

#name
pars_detect$Model<-"Nmixture"
```

###Assess Convergence


```r
###Chains
ggplot(pars_detect[pars_detect$par %in% c("detect","alpha","beta1","beta2"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

<img src="figureObserved/unnamed-chunk-33-1.png" title="" alt="" style="display: block; margin: auto;" />

###Hierarcichal Posteriors


```r
ggplot(pars_detect[pars_detect$par %in% c("beta1_mu","beta2_mu","alpha_mu","sigma_alpha","beta1_sigma","sigma_beta2","dprior","sigma_detect"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Trait-matching regression") + facet_wrap(~par,scales="free")
```

<img src="figureObserved/unnamed-chunk-34-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
#Bind together the two models
parsObs<-rbind(pars_detect,pars_dniave)
```

##Posteriors


```r
###Posterior Distributions
ggplot(parsObs[parsObs$par %in% c("detect","alpha","beta1","beta2"),],aes(x=estimate,fill=Model)) + geom_histogram(position='identity') + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() 
```

<img src="figureObserved/unnamed-chunk-36-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
#Detection figure
ggplot(parsObs[parsObs$par %in% c("detect"),],aes(x=as.factor(species),y=estimate,fill=Model)) + geom_violin() + ggtitle("Estimate of parameters") + theme_bw() + ggtitle("Detection Probability") +facet_wrap(~Model,scales="free") 
```

<img src="figureObserved/unnamed-chunk-37-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
pars_detect<-merge(pars_detect,jagsIndexBird,by.x="species",by.y="jBird",all.x=T)

ggplot(pars_detect[pars_detect$par %in% c("detect"),],aes(x=estimate)) + geom_histogram() + ggtitle("Posterior Distribution") + theme_bw() + facet_wrap(~Hummingbird,ncol=5) + xlab("Probability of Detection")
```

<img src="figureObserved/unnamed-chunk-37-2.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/DetectionProb.jpg",dpi=300,height=7,width=11)
```


```r
ggplot(parsObs[parsObs$par %in% c("beta1_mu","beta2_mu","alpha_mu","sigma_alpha","beta1_sigma","sigma_beta2","dprior","sigma_detect"),],aes(x=estimate,fill=Model)) + geom_histogram() + ggtitle("Trait matching regression parameters") + facet_wrap(~par,scale="free",nrow=2) + theme_bw() 
```

<img src="figureObserved/unnamed-chunk-38-1.png" title="" alt="" style="display: block; margin: auto;" />

###Overall predicted relationship 

Calculated predicted visitation rates

Does accounting for non-independence and detection change our estimate of trait matching?


```r
castdf<-dcast(parsObs[parsObs$par %in% c("beta1_mu","beta2_mu","alpha_mu"),], Model+Chain + Draw~par,value.var="estimate")

castdf<-split(castdf,castdf$Model)
```

## Trait+Abundance


```r
predy<-rbind_all(lapply(castdf,function(i){
  #calculate trajectory and append model
  pr<-trajF(alpha=i$alpha_mu,beta1=i$beta1_mu,beta2=0,trait=indat$Traitmatch,resources=indat$scaledR)
  pr$Model<-unique(i$Model)
  return(pr)
  }))

fplot<-ggplot(data=predy,aes(x=trait)) + geom_ribbon(aes(ymin=lower,ymax=upper,fill=Model),alpha=0.5)  + geom_line(aes(y=mean,col=Model),size=.4,linetype="dashed") + theme_bw() + ylab("Daily Interactions") + xlab("Difference between Bill and Corolla Length") + geom_point(data=indat,aes(x=Traitmatch,y=Yobs),size=.5,alpha=.5) + labs(fill="Model",col="Model") + ggtitle("Traits+Abundance") + scale_fill_manual(values=c("grey70","black"))+ scale_color_manual(values=c("grey70","black"))
fplot
```

<img src="figureObserved/unnamed-chunk-40-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/BothObs.svg",height=5,width=7)
ggsave("Figures/BothObs.jpg",heigh=5,width=7,dpi=300)
```

# Does including abundance change the biological inference?

We are only interested in the inference for the observed detection model data.

## Trait only


```r
castdf<-dcast(pars_detect_traits[pars_detect_traits$par %in% c("beta1_mu","alpha_mu"),], Chain + Draw~par,value.var="estimate")

predy_traits<-trajF(alpha=castdf$alpha_mu,beta1=castdf$beta1_mu,beta2=0,trait=indat$Traitmatch,resources=indat$Flowers)

tplot<-ggplot(data=predy_traits,aes(x=trait)) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)  + geom_line(aes(y=mean),size=.4,linetype="dashed") + theme_bw() + ylab("Daily Interactions") + xlab("Difference between Bill and Corolla Length") + geom_point(data=indat,aes(x=Traitmatch,y=Yobs),size=.5,alpha=.5) + labs(fill="Model",col="Model") + ggtitle("Traits")
tplot
```

<img src="figureObserved/unnamed-chunk-41-1.png" title="" alt="" style="display: block; margin: auto;" />



```r
allpred<-list('Traits+Abundance'=predy[predy$Model=="Nmixture",!colnames(predy) %in% "Model"],Traits=predy_traits)

allpred<-melt(allpred,id.vars=colnames(predy_traits))

allplot<-ggplot(data=allpred,aes(x=trait)) + geom_ribbon(aes(ymin=lower,ymax=upper,fill=L1),alpha=0.7)  + geom_line(aes(y=mean,col=L1),size=.4,linetype="dashed") + theme_bw() + ylab("Daily Interactions") + xlab("Difference between Bill and Corolla Length") + geom_point(data=indat,aes(x=Traitmatch,y=Yobs),size=.5,alpha=.7) + labs(fill="Model",col="Model") + ggtitle("N-mixture Model Comparison") + scale_fill_manual(values=c("black","grey70")) + scale_color_manual(values=c("black","grey80"))
allplot
```

<img src="figureObserved/unnamed-chunk-42-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/AbundanceBothPlot.jpeg",height=4,width=7,dpi=300)
```

##Species Predictions


```r
castdf<-dcast(parsObs[parsObs$par %in% c("beta1","beta2","alpha"),], species +Chain +Model+ Draw~par,value.var="estimate")

#Turn to species level
castdf$species<-factor(castdf$species,levels=1:max(as.numeric(castdf$species)))

species.split<-split(castdf,list(castdf$species,castdf$Model))

species.traj<-list()

for(d in 1:length(species.split)){
  x<-species.split[[d]]
  #species name
  index<-jagsIndexBird[unique(x$species),"Hummingbird"]
  
  #range of trait distances
  tsp<-indat %>% filter(Hummingbird==index) %>% .$Traitmatch
  
  #Range of abundances
    fsp<-indat %>% filter(Hummingbird==index) %>% .$Flowers
    
  species.traj[[d]]<-trajF(alpha=x$alpha,beta1=x$beta1,beta2=0,trait=tsp,resources=fsp)
}

names(species.traj)<-names(species.split)

species.traj<-melt(species.traj,id.var=colnames(species.traj[[1]]))

#split out names and model
species.traj[,c("Index","Model")]<-colsplit(species.traj$L1,"\\.",c("Index","Model"))

spe<-merge(species.traj,jagsIndexBird,by.x="Index",by.y="jBird")

#match colnames

#plot and compare to original data
ggplot(data=spe[,],aes(x=trait)) + geom_point(data=indat,aes(x=Traitmatch,y=Yobs)) + geom_ribbon(aes(ymin=lower,ymax=upper,fill=Model),alpha=0.2)  + geom_line(aes(y=mean,col=Model),size=1) + theme_bw() + ylab("Interactions") + xlab("Difference between Bill and Corolla Length") + facet_wrap(~Hummingbird,scales="free",ncol=3)+ labs(fill="Model")  + ylab("Interactions per day")
```

<img src="figureObserved/unnamed-chunk-43-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/SpeciesPredictionsBoth.jpg",dpi=300,height=8,width=10)
```

##Species Predictions - Abundance


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

<img src="figureObserved/unnamed-chunk-44-1.png" title="" alt="" style="display: block; margin: auto;" />

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
##                  Hummingbird mean lower upper
## 1             Andean Emerald 34.5   4.5  78.8
## 2         Booted Racket-tail 15.7   4.1  33.8
## 3                 Brown Inca 81.9  72.4  89.6
## 4        Buff-tailed Coronet  0.0   0.0   0.0
## 5              Collared Inca  0.1   0.1   0.3
## 6          Crowned Woodnymph 44.2   0.3  85.7
## 7    Fawn-breasted Brilliant  0.0   0.0   0.0
## 8          Gorgeted Sunangel  0.0   0.0   0.0
## 9    Green-crowned Brilliant  0.0   0.0   0.0
## 10   Green-fronted Lancebill 63.5   1.1  86.0
## 11             Hoary Puffleg  0.0   0.0   0.0
## 12    Purple-bibbed Whitetip  6.5   0.0  63.7
## 13 Rufous-tailed Hummingbird  0.0   0.0   0.1
## 14      Speckled Hummingbird  0.0   0.0   0.0
## 15    Stripe-throated Hermit 67.1  50.9  80.6
## 16      Tawny-bellied Hermit  0.2   0.1   0.4
## 17       Violet-tailed Sylph  0.0   0.0   0.0
## 18  Wedge-billed Hummingbird  0.1   0.0   0.4
## 19    White-whiskered Hermit 63.9  54.7  72.7
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

#get the 0.5 line
dpn<-function(t,p){
  n<-(1 - (1-t))/(p/100)
  return(n)
}

#for each bird get the upper and middle estimate for 50% chance.
daydf<-list()
for (x in 1:nrow(tab)){
  mean_day=dpn(t=0.5,tab$mean[x])
  lower_day=dpn(t=0.5,tab$lower[x])
  upper_day=dpn(t=0.5,tab$upper[x])
  daydf[[x]]<-data.frame(L1=tab$Hummingbird[x],mean=mean_day,lower=lower_day,upper=upper_day)
}
daydf<-rbind_all(daydf)

ggplot(md) + geom_ribbon(aes(x=Days,y=mean,ymin=lower,ymax=upper)) + geom_line(aes(x=Days,fill=L1,y=mean,ymin=lower,ymax=upper)) + facet_wrap(~L1,nrow=4,scale="free_x")  + ylab("Probability of detecting a interaction") + scale_fill_discrete(guide="none") + theme_bw() + scale_x_continuous(breaks=seq(0,10,2),limits=c(0,10))+ geom_rect(fill='grey',data=daydf,alpha=0.4,aes(xmax=upper,xmin=lower,ymin=0,ymax=Inf)) + ylim(0,1)
```

<img src="figureObserved/unnamed-chunk-46-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/DetectionDays.jpeg",height=7,width=9,dpi=300) 
```

The number of days it would take to have 50% confidence you have sampled enough to capture known interactions is the x axis value where the dotted line hits the curve.


```r
sampling<-indatraw %>% group_by(Hummingbird) %>% summarize(Obs=length(Hummingbird))

tabD<-merge(tab,sampling,by="Hummingbird")
ggplot(tabD,aes(x=Obs,ymin=lower,ymax=upper,y=mean)) + geom_pointrange() + labs(y="Detectability",x="Detections") + geom_text(aes(label=Hummingbird),vjust=2) + theme_bw() + xlim(0,175)
```

<img src="figureObserved/unnamed-chunk-47-1.png" title="" alt="" style="display: block; margin: auto;" />

##DIC Table


```r
m2_niave$BUGSoutput$DIC
```

```
## [1] 12971.63
```

```r
m2$BUGSoutput$DIC
```

```
## [1] 93415.18
```

#Predicted versus Observed Data


```r
m<-max(jTraitmatch)-jTraitmatch

mat<-indat %>% group_by(jBird,jPlant) %>% summarize(n=sum(Yobs))
true_state<-acast(mat,jBird~jPlant,fill=0)
```

###Test Statistic

Chisquared statistic is (observed-expected)^2/(expected + 0.5), we add the 0.5 to avoid dividing by 0. For each combination of birds and plants, predicted versus  observed for each data point.

#Compare using true known interactions

## Poisson GLMM


```r
#Discrepancy function
chisq<-function(o,e){(o-e)^2/(e+0.5)}

N_niave<-pars_dniave[ pars_dniave$par %in% "ynew",]

bydraw<-split(N_niave,list(N_niave$Chain,N_niave$Draw))

#remove large N matrix
rm(pars_dniave,N_niave)
gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  5891341 314.7    9306876  497.1   9306876  497.1
## Vcells 87440519 667.2  175292416 1337.4 174607305 1332.2
```

```r
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

#remove large N matrix
rm(N)
gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  5908420 315.6    9306876  497.1   9306876  497.1
## Vcells 90130508 687.7  175292416 1337.4 175291958 1337.4
```

```r
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

#Nmixture with detection
mocc<-melt(occ_matrix)
colnames(mocc)<-c("Bird","Plant","Nmixture","Iteration")
simdat<-merge(mocc,mmat,by=c("Bird","Plant"))

#Nmixture with nodetection
moccd<-melt(occ_nodetect_matrix)
colnames(moccd)<-c("Bird","Plant","Poisson GLMM","Iteration")

simdat<-merge(simdat,moccd,by=c("Bird","Plant","Iteration"))

simdat<-melt(simdat,measure.vars = c("Nmixture","Poisson GLMM"))

ggplot(simdat,aes(x=True_State,y=value,col=variable)) + geom_point() + geom_abline() + labs(col="Model") + ylab("Predicted State") + xlab("True State") + theme_bw() + facet_wrap(~variable)
```

<img src="figureObserved/unnamed-chunk-52-1.png" title="" alt="" style="display: block; margin: auto;" />

##Summary of discrepancy of predicted matrices


```r
#Nmixture without detection
occno_disc<-sapply(occ_nodetect,function(x) mean(x))

#Nmixture with detection
occ_disc<-sapply(occ,function(x) mean(x))

#compared to bayesian
ggplot() + xlab("Chi-squared Discrepancy") + geom_histogram(data=data.frame(occ_disc),aes(x=occ_disc),fill="red",alpha=.6) + theme_bw() +geom_vline(aes(xintercept=mean(occ_disc)),linetype="dashed",col="red") + geom_histogram(data=data.frame(occno_disc),aes(x=occno_disc),fill="orange",alpha=.6) + geom_vline(aes(xintercept=mean(occno_disc)),linetype="dashed",col="orange")
```

<img src="figureObserved/unnamed-chunk-53-1.png" title="" alt="" style="display: block; margin: auto;" />

##Comparison of summary statistics for all three approaches


```r
d<-list(Nmixture=occ,Poisson_GLM=occ_nodetect)
d<-melt(d)
colnames(d)<-c("Bird","Plant","value","Iteration","Model")

d %>% group_by(Model,Iteration) %>% summarize(mean=mean(value),sd=sd(value),sum=sum(value)) %>% group_by(Model) %>% summarize(mean_mean=round(mean(mean),2),mean_sd=round(sd(mean),2),mean_sum=round(mean(sum),2))
```

```
## Source: local data frame [2 x 4]
## 
##         Model mean_mean mean_sd mean_sum
##         (chr)     (dbl)   (dbl)    (dbl)
## 1    Nmixture     66.66    7.46 51929.29
## 2 Poisson_GLM      3.82    0.49  2979.44
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

<img src="figureObserved/unnamed-chunk-56-1.png" title="" alt="" style="display: block; margin: auto;" />

##Generate Networks


```r
castdf<-dcast(parsObs[parsObs$par %in% c("beta1","alpha"),], species +Model+Chain + Draw~par,value.var="estimate")

#Turn to 
castdf$species<-factor(castdf$species,levels=1:max(as.numeric(castdf$species)))

species.split<-split(castdf,list(castdf$species,castdf$Model),drop = T)

species.traj<-lapply(species.split,function(dat){
  index<-unique(dat$species)
  
  #get data for those species
  billd<-indat[indat$jBird %in% index,]
  
  d<-data.frame(alpha=dat$alpha,beta1=dat$beta1,beta2=0)
  
  #fit regression for each input estimate
  sampletraj<-list()
  
  for (y in 1:nrow(d)){
    v=exp(d$alpha[y] + d$beta1[y] * billd$Traitmatch) * billd$scaledR
    
    sampletraj[[y]]<-data.frame(x=as.numeric(billd$Traitmatch),y=as.numeric(v),jBird=billd$jBird,jPlant=billd$jPlant,Model=unique(dat$Model))
  }
  
  sample_all<-rbind_all(sampletraj)
})

species.traj<-rbind_all(species.traj)



species.mean<-species.traj %>% group_by(jBird,jPlant,Model) %>% summarize(Traitmatch=unique(x),phi=mean(y))

species.mean<-merge(species.mean,indat[,colnames(indat) %in% c("jBird","jPlant","jTime","Hummingbird","Iplant_Double")])

#get corolla sizes
species.mean<-merge(species.mean,fl.morph,by.x="Iplant_Double", by.y="Group.1")

#bill order
ord<-hum.morph %>% arrange(Total_Culmen) %>% .$English
species.mean$Hummingbird<-factor(species.mean$Hummingbird,levels=ord)

#add level to hum.morph to match naming convention
species.mean<-merge(species.mean,hum.morph[,c("English","Total_Culmen")],by.x="Hummingbird",by.y="English")


#ggplot(species.mean) + geom_density2d(aes(x=TotalCorolla,y=lambda,col=as.factor(BAll_Flowers))) + theme_bw() + facet_wrap(~Hummingbird,scales="free",ncol=3)+ scale_color_manual("Resources Availability",labels=c("Low","High"),values=c("Blue","Red")) + ggtitle("2D Density Plots") + geom_vline(aes(xintercept=Total_Culmen),linetype='dashed')


#Niche Breadth 


species.mean<-species.traj %>% group_by(jBird,jPlant,Model) %>% summarize(Traitmatch=unique(x),phi=mean(y),phi_low=quantile(y,0.05),phi_high=quantile(y,0.95))

#merge names
species.mean<-merge(species.mean,jagsIndexBird)
species.mean<-merge(species.mean,jagsIndexPlants)

#get corolla sizes
species.mean<-merge(species.mean,fl.morph,by.x="Iplant_Double", by.y="Group.1")

#bill order
ord<-hum.morph %>% arrange(Total_Culmen) %>% .$English
species.mean$Hummingbird<-factor(species.mean$Hummingbird,levels=ord)

#add level to hum.morph to match naming convention
species.mean<-merge(species.mean,hum.morph[,c("English","Total_Culmen")],by.x="Hummingbird",by.y="English")

ggplot(species.mean) + geom_ribbon(alpha=0.4,aes(x=TotalCorolla,ymin=phi_low,ymax=phi_high,fill=as.factor(Model))) + theme_bw() + facet_wrap(~Hummingbird,scales="free",ncol=4)+ ggtitle("Niche Breadth") + geom_vline(aes(xintercept=Total_Culmen),linetype='dashed') + geom_line(aes(x=TotalCorolla,y=phi,fill=as.factor(Model))) + ylab("Daily Interaction Rate") + xlab("Corolla Length (mm)") + scale_fill_discrete("Resource Availability")
```

<img src="figureObserved/unnamed-chunk-57-1.png" title="" alt="" style="display: block; margin: auto;" />

```r
ggsave("Figures/NicheBreadth.jpeg",height=6,width=9)
```


#Network Statistics

Given the uncertainty in species interactions, what do emergant statistics look like?


```r
#Split by resource
nsplit<-split(species.mean,species.mean$Model)

makeN<-function(x){
  
  #input matrix
  aggm<-matrix(nrow=nrow(jagsIndexBird),ncol=nrow(jagsIndexPlants),data=0)
  for (j in 1:nrow(x)){
    aggm[x[j,"jBird"],x[j,"jPlant"]]<-rpois(1,lambda=x[j,"phi"])
  }
  #calculate network statistic
  nstat<-networklevel(aggm,index=c("connectance","nestedness"))
}

nstat<-lapply(nsplit,function(x){
  netstat<-melt(t(sapply(1:100,function(k) makeN(x)))) 
  colnames(netstat)<-c("Iteration","Metric","value")
  return(netstat)
})

nstat<-melt(nstat,colnames(nstat[[1]]))

ggplot(nstat,aes(x=value,fill=L1)) + geom_density(alpha=0.6) + facet_wrap(~Metric,scales='free',nrow=2) + scale_fill_discrete("Model")
```

<img src="figureObserved/unnamed-chunk-58-1.png" title="" alt="" style="display: block; margin: auto;" />


```r
gc()
```

```
##             used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   5930771 316.8    9306876  497.1   9306876  497.1
## Vcells 127259497 971.0  210430899 1605.5 210430615 1605.5
```

```r
save.image("Observed.Rdata")
```


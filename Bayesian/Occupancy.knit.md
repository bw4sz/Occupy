---
title: "Combining automated monitering and bayesian occupancy models for detecting rare plant-animal interactions"
author: "Ben Weinstein"
date: "Tuesday, February months, 2015"
output: html_document
---



#Aim

Network ecology is a rapidly developing field that uses graph theory to represents interactions among plants and pollinators, hosts and parasites, and herbivores and hosts. One aim of network ecology is to quantify the rate of interaction among partners, their relative specialization, and the robustness of these interactions to perturbation. Considerable work has focused on the effect sampling on network structure and how best to measure network properties (Bluthgen, Bascompte, Jordano). Considerably less work has focused on a more fundamental problem of network ecology: Our ability to detect interactions is imperfect in time and space. This limitation is a fundamental assumption in wildlife ecology, where occupancy models are common(). In occupancy modeling, repeated surveys of the same site are used to estimate the probability of detection given a species probability of occurrence. The probability of occurrence is a latent variable which cannot be directly infered from the data. If we did not see a species at a site, was it absent or undetected? By analogy, if we failed to find a pollinator interacting with a plant, is it due to detection or is it a 'forbidden link' due mismatch in species ecology? The goal of this paper is to use a large dataset on plant-hummingbird interactins from a tropical montane forest to evaluate the effect of imperfect detection on estimating species interactions and network properties. We will then seperate estimate environmental and morphological covariates to both the probability of interaction ($\psi$) as well as the probability of detection (p).

#Background

  Collecting network data is time consuming and few studies have the resources to adopt a repeated measures sampling design. The hummingbird dataset in this paper was collected using time-lapse cameras which turn on at dawn and dusk automatically. In addition, more traditional hummingbird transects were used to survey flower visitation. Combining these sampling types, we have a very large network with over 3,000 interactions (Weinstein XXXX).

#Similiar work

Unfortunately someone already did this:
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069200

How can i take it a step further to actually make ecological inference?

#Approach

Potential comparisons

* Transects versus cameras
* Corrected versus uncorrected network statistics
* Simulations - networks with imperfect detection - what are the biases?


**In this study we want to estimate three quantitites:**

$$\hat{Y}_{i,j} = \text{Detection  of  Hummingbird i  at  Flower j}$$

$$ \psi_{i,j} = \text{the probability a flower is visited by a hummigbird}  $$

$$ p_{i,j} = \text{the probability of detecting a bird-flower interaction given that it occurs} $$

#Simulations

Consider a matrix of interactions, where we assume interactions among plants and hummingbirds are independent.

### Parameters

  * 2 hummingbird
  * Twenty plants
  * Known occupancy $\psi = U(.2,1)$
    * At very low occupancies, the interaction is too rare to estimate.
  * Imperfect detection $p = U(0,1)$ 
  * 24 Months replicates

Each species gets their own occupancy and imperfect detection which is drawn from uniform distribution. The goal is **recover these parameters**.


```r
#Number of hummingbird species
h_species=5
plant_species=10
months=24
detection=runif(h_species,0,1)
occupancy=runif(h_species,0,1)
dat<-sapply(occupancy,function(x){rbinom(plant_species,1,occupancy)})
```

### True Interaction Matrix


```r
#Reshape into a nicer format
mdat<-melt(dat)
colnames(mdat)<-c("Plant","Hummingbird","True_State")
ggplot(mdat,aes(x=as.factor(Plant),y=as.factor(Hummingbird),fill=as.factor(True_State))) + geom_tile() + scale_fill_manual(labels=c(0,1),values=c("White","Black")) + labs(x="Plant",y="Hummingbird",fill="Presence") + ggtitle("True State")
```

<img src="figure/unnamed-chunk-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="576" style="display: block; margin: auto;" />

##Simulate detection 


```r
#for each species loop through and create a replicate dataframe
obs<-list()

for (x in 1:h_species){
  
  hdat<-mdat[mdat$Hummingbird==x,]
  
  #How many detections?
  size<-nrow(hdat[hdat$True_State==1,])
  
  #Simualte detections for species with presences
  pres<-replicate(months,rbinom(size,1,detection[x]))
  #fill in 0's for species with absences
  ab<-replicate(months,rep(0,plant_species-size))
  obs[[x]]<-rbind(pres,ab)
}

mobs<-melt(obs)
colnames(mobs)<-c("Plants","Month","Detection","Hummingbird")

mobs<-acast(mobs,Plants~Hummingbird~Month,value.var="Detection")
```

Data is stored in a multidimensional array
Hummingbird 1 on Plant 1 - True State is Present


```r
dat[1,1]
```

```
## [1] 1
```

Observed state from the simulated transect data
**Format is [Plants,Hummingbird,Month]**


```r
mobs[1,,1]
```

```
## 1 2 
## 1 1
```


##Hierarchical model for the interaction of each hummingbird on each plant species

$$ Y_{i,j} \sim Bern(sightp_{i,j})$$
$$sightp_{i,j} = present_{i,j} * detect_i$$
$$present_{i,j} \sim Bern(presentp) $$
$$presentp_{i,j}<-occ_i $$
















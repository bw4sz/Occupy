#extract and create a dataframe of posteriors

extract_par<-function(x){
#extract desired info from the models
parsO<-melt(x$BUGSoutput$sims.array)
colnames(parsO)<-c("Draw","Chain","parameter","estimate")
parsO<-parsO[!parsO$Draw %in% 1:(max(parsO$Draw)-3000/x$BUGSoutput$n.chains),]

#label species and plants
l<-levels(parsO$parameter)

#parameters to save
totrack<-x$parameters.to.save

sp_pl<-data.frame(parameter=l,str_match(l,pattern="(\\d+),(\\d+)")[,-1],par=str_extract(l,"\\w+"))

colnames(sp_pl)<-c("parameter","species","plant","par")

sp_pl[is.na(sp_pl$species),c("species")]<-str_extract(sp_pl[is.na(sp_pl$species),"parameter"],"\\d+")

#merge levels
pars<-merge(parsO,sp_pl)

#take out deviance
pars<-pars[!pars$par %in% "deviance",]
return(pars)
}

#fits a curve for given poisson function

trajF<-function(alpha,beta,x){
  fdat<-data.frame(alpha=alpha,beta=beta)
  
  #fit regression for each input estimate
  sampletraj<-list()
  for (s in 1:nrow(fdat)){
    a<-fdat$alpha[s]
    b<-fdat$beta[s]
    b2<-fdat$beta2[s]
    yp=exp(a + (b*x))
    
    #compute pred value
    sampletraj[[s]]<-data.frame(x=x,y=yp)
  }
  
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
  return(predy)
}

#fits a chisquared residual for a given poisson function

trajState<-function(alpha,beta,x,observed){
  
  #Bind together
  fdat<-data.frame(alpha=alpha,beta=beta)
  
  #fit regression for each input estimate
  sampletraj<-list()
  for (s in 1:nrow(fdat)){
    a<-fdat$alpha[s]
    b<-fdat$beta[s]
    yp=exp(a + (b*x$value))
    
    #compute pred value
    state<-data.frame(x,State=rpois(length(yp),yp))
    
    #merge with observed state
    mstate<-merge(state,observed,by=c("Bird","Plant"))
    
    #Compute chisquared
    csq<-sum((mstate$Y-mstate$State)^2/(mstate$State+0.5))
    
    sampletraj[[s]]<-csq
  }
  
  #return as a vector
  return(unlist(sampletraj))
}

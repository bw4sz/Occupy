#extract and create a dataframe of posteriors

extract_par<-function(x,data=obs,Bird="Bird",Plant="Plant"){
  #extract desired info from the models
  n<-dim(x$BUGSoutput$sims.array)[1]
  parsO<-melt(x$BUGSoutput$sims.array[max(0,(n-500)):n,,])
  colnames(parsO)<-c("Draw","Chain","parameter","estimate")
  
  #label species and plants
  l<-levels(parsO$parameter)
  
  #parameters to save
  totrack<-x$parameters.to.save
  
  #assign species index to ragged frame.
  sp_pl<-data.frame(parameter=l,species=as.numeric(str_match(l,pattern="\\[(\\d+)]")[,2]),par=str_extract(l,"\\w+"))
  
  #correct N samples
  i<-sp_pl$par %in% "ynew"
  
  #Species
  sp_pl[i,][,"species"]<-data[as.numeric(str_match(sp_pl[i,][,"parameter"],pattern="\\[(\\d+)]")[,2]),Bird]
  
  
  #Plant
  #add a NA plant columns
  sp_pl$plant<-NA
  sp_pl[i,][,"plant"]<-
    data[as.numeric(str_match(sp_pl[i,][,"parameter"],pattern="\\[(\\d+)]")[,2]),Plant]
  
  #merge levels, can be very large, do in pieces. 
  parsO<-inner_join(parsO,sp_pl) %>% filter(!par == "deviance")
  
  return(parsO)
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

trajHDI<-function(alpha,beta,x){
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
  predy<-group_by(sample_all,x) %>% summarise(lower=hdi(y)[[1]],upper=hdi(y)[[2]],mean=mean(y,na.rm=T))
  return(predy)
}


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

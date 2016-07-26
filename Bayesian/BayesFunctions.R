#extract and create a dataframe of posteriors

extract_par<-function(x,data=obs,Bird="Bird",Plant="Plant"){
  #extract desired info from the models
  parsO<-melt(x$BUGSoutput$sims.array)
  colnames(parsO)<-c("Draw","Chain","parameter","estimate")
  
  #label species and plants
  l<-levels(parsO$parameter)
  
  #parameters to save
  totrack<-x$parameters.to.save
  
  #assign species index to ragged frame.
  sp_pl<-data.frame(parameter=l,species=as.numeric(str_match(l,pattern="\\[(\\d+)]")[,2]),par=str_extract(l,"\\w+"))
  
  #correct N samples
  i<-sp_pl$par %in% "E"

  #Species
  sp_pl[i,][,"species"]<-data[as.numeric(str_match(sp_pl[i,][,"parameter"],pattern="\\[(\\d+)]")[,2]),Bird]
  
  #Plant
  #add a NA plant columns
  sp_pl$plant<-NA
  sp_pl[i,][,"plant"]<-data[as.numeric(str_match(sp_pl[i,][,"parameter"],pattern="\\[(\\d+)]")[,2]),Plant]
  
  #merge levels
  pars<-merge(parsO,sp_pl)
  
  #take out deviance
  pars<-pars[!pars$par %in% "deviance",]
  return(pars)
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
    yp=exp(a + b*x$value)
    
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

#sample trajectory for a given posterior
trajF<-function(alpha,beta1,beta2,beta3,trait,resources){
  g<-data.frame(alpha,beta1,beta2,beta3)
  
  #label rows
  g$id<-1:nrow(g)
  
  sampletraj<-g %>% group_by(id) %>% do(traj(.$alpha,.$beta1,trait=trait,resources,.$beta2,.$beta3)) %>% group_by(trait) %>% summarize(mean=mean(y),lower=quantile(y,0.05),upper=quantile(y,0.95))
  return(sampletraj)
}

#sample trajectory for a given posterior using quantile or hdi interval
traj<-function(alpha,beta1,beta2,beta3,trait,resources){

    #fit regression for each input estimate
    
    v=exp(alpha + beta1 * trait + beta2 * resources + beta3 * trait*resources)
    
    sampletraj<-data.frame(trait=trait,y=as.numeric(v))
  
  #Compute CI intervals
  return(sampletraj)
}

#sample trajectory for a given posterior using quantile or hdi interval
trajLog<-function(alpha,beta1,beta2,beta3,x,resources,type='quantile'){
  indat<-data.frame(alpha,beta1,beta2,beta3)
  
  #fit regression for each input estimate
  sampletraj<-list()
  
  for (y in 1:nrow(indat)){
    v=exp(indat$alpha[y] + indat$beta1[y] * x + indat$beta2[y] * resources + indat$beta3[y] * x*resources)
    
    sampletraj[[y]]<-data.frame(x=as.numeric(x),y=as.numeric(v))
  }
  
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  if(type=='quantile'){
    predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
  }
  if(type=='hdi'){
    predy<-group_by(sample_all,x) %>% summarise(lower=hdi(y)[[1]],upper=hdi(y)[[2]],mean=mean(y,na.rm=T))
  }
  return(predy)
}

#predicted y for logistic
trajLogistic<-function(alpha,beta1,beta2,beta3,x,resources){
  indat<-data.frame(alpha,beta1,beta2,beta3)
  
  #fit regression for each input estimate
  sampletraj<-list()
  
  for (y in 1:nrow(indat)){
    v=exp(indat$alpha[y] + indat$beta1[y] * x + indat$beta2[y] * resources + indat$beta3[y] * x*resources)
    
    sampletraj[[y]]<-data.frame(x=as.numeric(x),y=as.numeric(v))
  }
  
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
}

#calculate poisson interactions

intF<-function(alpha,beta1,beta2,beta3,x,resources){
  indat<-data.frame(alpha,beta1,beta2,beta3)
  
  #fit regression for each input estimate
  sampletraj<-list()
  
  for (y in 1:nrow(indat)){
    v=indat$beta2[y] + indat$beta3[y]  * x
    sampletraj[[y]]<-data.frame(x=as.numeric(x),y=as.numeric(v))
  }
  
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
  return(predy)
}

#plots
#converge of chains
chainplot<-function(pars,param,title){
  ggplot(pars[pars$par %in% param,],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_wrap(~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle(title)  
}

#posteriors
tracegplot<-function(pars,param,title){
  ggplot(pars[pars$par %in% param,],aes(x=estimate)) + geom_histogram() + ggtitle("Estimate of Intercept") + theme_bw() + ggtitle(title)
}



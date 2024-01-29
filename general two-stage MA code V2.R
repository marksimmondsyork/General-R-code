# general two-stage analysis by period

two.stage.ma = function(data,param){
    
  pds = unique(data$period)
  npds = length(pds)
  trialnames = as.character(unique(data$trial))
  ntrials = length(trialnames)
  
  num.pats = tapply(param,factor(data$trtgrp):factor(data$period):factor(data$trial),function(x) length(x[is.na(x)==F]))
  num.pats = ifelse(is.na(num.pats)==T,0,num.pats)
  
  num.evs = tapply(param,factor(data$trtgrp):factor(data$period):factor(data$trial),function(x) sum(x[is.na(x)==F]))
  num.evs = ifelse(is.na(num.evs)==T,0,num.evs)
  
  n.t = matrix(num.pats[grep("Inv",names(num.pats))],ntrials,npds)
  n.c = matrix(num.pats[grep("Con",names(num.pats))],ntrials,npds)
  n.et = matrix(num.evs[grep("Inv",names(num.evs))],ntrials,npds)
  n.ec = matrix(num.evs[grep("Con",names(num.evs))],ntrials,npds)
  
  TE = seTE = matrix(NA,ntrials,npds)
  
  for (pd in 1:npds){
    
    m1 = metabin(n.et[,pd],n.t[,pd],n.ec[,pd],n.c[,pd],sm="RR")
    TE[,pd] = m1$TE
    seTE[,pd] = m1$seTE
  }

  out = data.frame(trial=rep(trialnames,npds),period=rep(pds,each=ntrials),RR=as.numeric(TE),seRR=as.numeric(seTE),n.t=as.numeric(n.t),n.c=as.numeric(n.c),n.et=as.numeric(n.et),n.ec=as.numeric(n.ec))
  return(out)
  
}

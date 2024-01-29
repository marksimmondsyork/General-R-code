# general two-stage MD analysis

two.stage.md.ma = function(data,param,name,forest=TRUE, fname="two-stage MD "){
  
  pds = unique(data$period)
  npds = length(pds)
  trialnames = as.character(unique(data$trial))
  ntrials = length(trialnames)
  
  num.pats = tapply(param,factor(data$trtgrp):factor(data$period):factor(data$trial),length)
  num.pats = ifelse(is.na(num.pats)==T,0,num.pats) 
  ch = tapply(param,factor(data$trtgrp):factor(data$period):factor(data$trial),function(x) mean(x[is.na(x)==F]))  
  sd = tapply(param,factor(data$trtgrp):factor(data$period):factor(data$trial),function(x) sd(x[is.na(x)==F]))
  
  num.pats.bmp = matrix(num.pats[1:(ntrials*npds)],ntrials,npds)
  num.pats.ic = matrix(num.pats[(ntrials*npds+1):(2*ntrials*npds)],ntrials,npds)
  ch.bmp = matrix(ch[1:(ntrials*npds)],ntrials,npds)
  ch.ic = matrix(ch[(ntrials*npds+1):(2*ntrials*npds)],ntrials,npds)
  sd.bmp = matrix(sd[1:(ntrials*npds)],ntrials,npds)
  sd.ic = matrix(sd[(ntrials*npds+1):(2*ntrials*npds)],ntrials,npds)
  
  TE = seTE = tau = I2 = I2low = I2high = rep(NA,npds)
  
  for (pd in 1:npds){
    
    m1 = metacont(num.pats.bmp[,pd],ch.bmp[,pd],sd.bmp[,pd],num.pats.ic[,pd],ch.ic[,pd],sd.ic[,pd])
    TE[pd] = m1$TE.random
    seTE[pd] = m1$seTE.random
    tau[pd] = m1$tau
    I2[pd] = summary(m1)$I2$TE
    I2low[pd] = summary(m1)$I2$lower
    I2[pd] = summary(m1)$I2$upper
    if (forest==TRUE){
      filename = paste("figures/",  fname, name," ",pds[pd]," forest.png",sep="")
      png(filename,width=1000,height=500)
      forest(m1,studlab=trialnames,xlab="Mean difference")
      dev.off()
    }
  }
  
  out = data.frame(pds,TE,seTE,tau,I2,I2low,I2high)
  return(out)
  
}
  

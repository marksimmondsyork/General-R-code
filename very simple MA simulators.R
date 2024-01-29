library(tidyverse)
library(magrittr)
library(meta)

##############################################
# MA of single estimate (cont.)

simple_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,thetain=NULL){
  
  # parameters
  if (is.null(thetain))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- theta.in

  # IPD
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  epsilon.i <- rnorm(Nt, epsilon, 0.1*epsilon)
  Y <- pmap(list(trialsize,theta.i,epsilon.i), rnorm)
  
  # summary data
  Ybar <- map_dbl(Y,mean)
  Ysd <- map_dbl(Y,sd)
  Ysem <- Ysd / sqrt(trialsize)

  out <- data.frame(N=trialsize,Ybar,Ysd,Ysem,mu.true=mu,tau2.true=tau2,epsilon)
  return(out)
}



#################################################
# MA of proportions (binary)

simple_binary_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0.1,tau2=0.05,LOin=NULL){

  if (is.null(LOin))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- LOin

  p.i <- exp(theta.i) / (1+exp(theta.i))
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  Y <- map2_dbl(trialsize, p.i, ~rbinom(1,.x,.y))
  
  # summaries (Log odds)
  p.est <- Y / trialsize
  Ybar <- log(p.est / (1-p.est))
  Ysem <- sqrt( 1 / (trialsize * p.est * (1-p.est)) )

  out <- data.frame(N=trialsize,Ybar,Ysem,mu.true=mu,tau2.true=tau2)
  return(out)
}

#########################################################
# version with "small sample" effects

bias_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,trend=0.001){
  
  require(purrr)

  theta.i <- rnorm(Nt, mu, sqrt(tau2))

  # IPD
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  epsilon.i <- rnorm(Nt, epsilon, 0.1*epsilon)
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.i - trend * trialsize.c
 
  Y <- pmap(list(trialsize,theta.i.bias,epsilon.i), rnorm)
  
  # summary data
  Ybar <- map_dbl(Y,mean)
  Ysd <- map_dbl(Y,sd)
  Ysem <- Ysd / sqrt(trialsize)
  
  out <- data.frame(N=trialsize,Ybar,Ysd,Ysem,mu.true=mu,tau2.true=tau2,epsilon)
  return(out)
}


################################################
# run standard MAs on data

run_MA <- function(data){
  
  Nt <- length(data$Ybar)
  
  # standard MAs
  ma.std <- metagen(Ybar, Ysem, data=data, prediction=T)
  ma.hakn <- metagen(Ybar, Ysem, data=data, hakn=T)
  wt.fe <- ma.std$w.fixed / sum(ma.std$w.fixed)
  wt.re <- ma.std$w.random / sum(ma.std$w.random)
  cov.fe <- cov(wt.fe,ma.std$TE)
  cov.re <- cov(wt.re,ma.std$TE)
  
  # equal weighting
  tau2.est <- ma.std$tau^2
  mu.eqw <- mean(data$Ybar)
  se.eqw <- sqrt(sum(data$Ysem^2) / Nt^2)
  se.eqwre <- sqrt(sum(data$Ysem^2 + tau2.est) / Nt^2)
  
  # sample size weighting
  mu.ssw <- sum(data$N * data$Ybar) / sum(data$N)
  se.ssw <- sqrt( sum(data$N^2 * data$Ysem^2) / (sum(data$N))^2 )
  se.sswre <- sqrt( sum(data$N^2 * (data$Ysem^2 + tau2.est)) / (sum(data$N))^2 )
  ssw.w <- data$N / sum(data$N)
  cov.ssw <- cov(ssw.w,ma.std$TE)
  
  out.std <- ma.std %$% data.frame(mu.fe=TE.fixed,se.fe=seTE.fixed,
                                   mu.re=TE.random,se.re=seTE.random,
                                   tau2.ma=tau^2, cov.fe,cov.re,
                                   lower.fixed,upper.fixed,
                                   lower.random,upper.random,
                                   lower.predict,upper.predict)
  
  out.hakn <- ma.hakn %$% data.frame(mu.hakn=TE.random,se.hakn=seTE.random,
                                     tau2.hakn=tau^2,
                                     lower.hakn=lower.random,
                                     upper.hakn=upper.random)
  
  out.other <- data.frame(mu.eqw,se.eqw,se.eqwre,
                          lower.eqw=(mu.eqw-1.96*se.eqw),upper.eqw=(mu.eqw+1.96*se.eqw),
                          lower.eqwre=(mu.eqw-1.96*se.eqwre),
                          upper.eqwre=(mu.eqw+1.96*se.eqwre),
                          mu.ssw,se.ssw,se.sswre,cov.ssw,
                          lower.ssw=(mu.ssw-1.96*se.ssw),upper.ssw=(mu.ssw+1.96*se.ssw),
                          lower.sswre=(mu.ssw-1.96*se.sswre),
                          upper.sswre=(mu.ssw+1.96*se.sswre)
                          )
  
  params <- data.frame(mu=data$mu.true[1], tau2=data$tau2.true[1])
  
  out <- data.frame(params,out.std,out.hakn,out.other)
  return(out)
  
}

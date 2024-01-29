# meta-analysis with random weights

library(tidyverse)
library(magrittr)
library(meta)

source("very simple MA simulators.R")

# MA function

ma_general <- function(theta,se,w){
  
  theta.est <- sum(w*theta) / sum(w)
  var.est.fe <- sum(w^2 * se^2) / (sum(w))^2
  cov.est <- cov(theta,w)
  corr.est <- cor(theta,w)
  
  Q <- sum(w*(theta-theta.est)^2)
  tau.num2 <- sum(w*se^2)-(sum(w^2*se^2)/sum(w))
  tau.denom <- sum(w) - (sum(w^2)/sum(w))
  tau2.est <- (Q - tau.num2) / tau.denom
  tau2.cut <- ifelse(tau2.est<0, 0, tau2.est)
  
  var.est.re <- sum(w^2 * (se^2 + tau2.cut)) / (sum(w))^2
  
  sd.est.fe <- sqrt(var.est.fe)
  sd.est.re <- sqrt(var.est.re)
  
  out <- data.frame(theta.est,
                    sd.est.fe,sd.est.re,
                    var.est.fe,var.est.re,
                    tau2.est,Q,cov.est,corr.est)
  return(out)
}

# Simulate Normal-Normal data (optional small sample effect)

ma_normnorm <- function(Nt,Nmin,Nmax,theta,tau2,SD=1,trend=0){
  
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  
  theta.true <- rnorm(Nt, theta, sqrt(tau2))
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.true - trend * trialsize.c
  
  se.theta <- SD / sqrt(trialsize)
  theta.est <- rnorm(Nt, theta.i.bias, se.theta)
  out <- data.frame(Trial=1:Nt, N=trialsize, theta=theta.est, se=se.theta)
  return(out)
}


# Simulate Binomial-Normal data (optional small sample effect)

ma_binnorm <- function(Nt,Nmin,Nmax,theta,tau2,p.c=0.5,trend=0){
  
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  armsize  <- floor(trialsize/2)
  
  theta.true <- rnorm(Nt, theta, sqrt(tau2))
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.true - trend * (0.01*trialsize.c)
  
  ev.c <- rbinom(Nt,armsize,p.c)
  odds.c <- p.c /(1-p.c)
  odds.e <- exp(theta.i.bias) * odds.c
  p.e <- odds.e / (1+odds.e)
  ev.e <- rbinom(Nt,armsize,p.e)
  res <- metabin(ev.e,armsize,ev.c,armsize)
  out <- data.frame(Trial=1:Nt, N=2*armsize, ev.e, N.e=armsize, ev.c, 
                    N.c=armsize,theta=res$TE,se=res$seTE)
  return(out)
}

# generic norm-norm

data.samp <- ma_normnorm(10,100,200,0,0,1)
data.samp <- ma_normnorm(10,100,200,0,1,1)

data.samp <- ma_normnorm(5,10,20,0,0,1)
data.samp <- ma_normnorm(5,10,20,0,1,1)

data.samp <- ma_normnorm(10,100,200,0,0,1,0.1)

w.rand <- map(1:100, ~runif(10, 0, 1))
w.rand.ss <- map(1:100, ~runif(10, 0.5*data.samp$N, 1.5*data.samp$N))
ma.rand.res <- map_dfr(w.rand, ~data.samp %$% ma_general(theta,se,.x))
ma.rand.ss.res <- map_dfr(w.rand.ss, ~data.samp %$% ma_general(theta,se,.x))
ma.classic.res <- data.samp %$% metagen(theta,se,sm="MD",method.tau="DL")

ggplot(ma.rand.res,aes(x=corr.est,y=theta.est)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=sd.est.fe)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=sd.est.re)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=tau2.est)) + geom_point()

ma.classic.res$seTE.fixed
ma.classic.res$seTE.random

ggplot(ma.rand.ss.res,aes(x=corr.est,y=theta.est)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=sd.est.fe)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=sd.est.re)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=tau2.est)) + geom_point()


# binom-norm
data.samp <- ma_binnorm(10,100,200,0,0,0.5,0)
data.samp <- ma_binnorm(10,100,200,0,1,0.5,0)

data.samp <- ma_binnorm(10,100,1000,0,0,0.5,0.1)
data.samp <- ma_binnorm(10,100,1000,0,1,0.5,0.1)

w.rand <- map(1:100, ~runif(10, 0, 1))
w.rand.ss <- map(1:100, ~runif(10, 0.5*data.samp$N, 1.5*data.samp$N))
ma.rand.res <- map_dfr(w.rand, ~data.samp %$% ma_general(theta,se,.x))
ma.rand.ss.res <- map_dfr(w.rand.ss, ~data.samp %$% ma_general(theta,se,.x))
ma.classic.res <- data.samp %$% metagen(theta,se,sm="OR",method.tau="DL")
ma.pipd.res <- data.samp %$% metabin(ev.e,N.e,ev.c,N.c,method="GLMM",sm="OR",model.glmm="UM.RS")
mega.trial <- data.samp %$% metabin(sum(ev.e),sum(N.e),sum(ev.c),sum(N.c),sm="OR")

ggplot(ma.rand.res,aes(x=corr.est,y=theta.est)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=sd.est.fe)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=sd.est.re)) + geom_point()
ggplot(ma.rand.res,aes(x=corr.est,y=tau2.est)) + geom_point()

ma.classic.res$TE.random
ma.classic.res$seTE.fixed
ma.classic.res$seTE.random
ma.pipd.res$TE.random
ma.pipd.res$seTE.fixed
ma.pipd.res$seTE.random

ggplot(ma.rand.ss.res,aes(x=corr.est,y=theta.est)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=sd.est.fe)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=sd.est.re)) + geom_point()
ggplot(ma.rand.ss.res,aes(x=corr.est,y=tau2.est)) + geom_point()

# regresison
Nt.c <- data.samp$N - mean(data.samp$N)
metareg(ma.pipd.res, Nt.c)
metareg(ma.classic.res, Nt.c)


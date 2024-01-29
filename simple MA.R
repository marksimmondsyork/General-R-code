library(tidyverse)
library(magrittr)
library(furrr)
library(meta)

##############################################
# MA of single estimate (cont.)

simple_MA <- function(Nt=10,N=100,mu=1,tau2=1,epsilon=1,thetain=NULL){

  if (is.null(thetain))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- theta.in
  Y <- map_dfc(theta.i, ~rnorm(N, ., epsilon))

  theta.hat <- as.numeric( summarise_all(Y,mean) )
  theta.hat.se <- as.numeric( summarise_all(Y, function(x) sd(x)/sqrt(length(x))) )

  ma.std <- metagen(theta.hat, theta.hat.se, prediction=T)
  ma.hakn <- metagen(theta.hat, theta.hat.se, hakn=T)

  out.std <- ma.std %$% data.frame(mu.fe=TE.fixed,se.fe=seTE.fixed,
                    mu.re=TE.random,se.re=seTE.random,
                    tau2=tau^2,
                    lower.fixed,upper.fixed,
                    lower.random,upper.random,
                    lower.predict,upper.predict,
                    cover.fixed=ifelse(lower.fixed<=mu & upper.fixed>=mu, 1, 0),
                    cover.random=ifelse(lower.random<=mu & upper.random>=mu, 1, 0))

  out.hakn <- ma.hakn %$% data.frame(mu.hakn=TE.random,se.hakn=seTE.random,
                                     tau2.hakn=tau^2,
                                     lower.hakn=lower.random,
                                     upper.hakn=upper.random,
                                     cover.hakn=ifelse(lower.random<=mu & upper.random>=mu, 1, 0))

  out <- data.frame(out.std,out.hakn)

  return(out)
}


#################################################
# MA of proportions (binary)

simple_binary_MA <- function(Nt=10,N=100,mu=0.1,tau2=0.05,LOin=NULL){

  if (is.null(LOin))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- LOin

  p.i <- exp(theta.i) / (1+exp(theta.i))
  Y <- rbinom(Nt, N, p.i)
  tsize <- rep(N, Nt)

  ma.std <- metaprop(Y, tsize, sm="PLOGIT", prediction=T)
  ma.hakn <- metaprop(Y, tsize, sm="PLOGIT", hakn=T)

  out.std <- ma.std %$%
    data.frame(mu.fe=TE.fixed,se.fe=seTE.fixed,
               mu.re=TE.random,se.re=seTE.random,
               tau2=tau^2,
               lower.fixed,upper.fixed,
               lower.random,upper.random,
               lower.predict,upper.predict,
               cover.fixed=ifelse(lower.fixed<=mu & upper.fixed>=mu, 1, 0),
               cover.random=ifelse(lower.random<=mu & upper.random>=mu, 1, 0))


  out.hakn <- ma.hakn %$%
    data.frame(mu.hakn=TE.random,se.hakn=seTE.random,
               tau2.hakn=tau^2,
               lower.hakn=lower.random,
               upper.hakn=upper.random,
               cover.hakn=ifelse(lower.random<=mu & upper.random>=mu, 1, 0))


  out <- data.frame(out.std,out.hakn)

  return(out)
}



#####################################################
# run evaluations

# cont.
ma.all <- replicate(10000, simple_MA(), simplify=F)
ma.all <- do.call(rbind, ma.all)

plan(multiprocess)
ma.all <- future_map_dfr(1:1000, simple_MA, .progress=T)
plan(sequential)

theta.in <- rnorm(10,1,1)
ma.fe <- replicate(10000, simple_MA(thetain=theta.in), simplify=F)
ma.fe <- do.call(rbind, ma.fe)

# binary
ma.all <- replicate(10, simple_binary_MA(), simplify=F)
ma.all <- do.call(rbind, ma.all)

LO.in <- rnorm(10,1,1)
ma.fe <- replicate(10000, simple_binary_MA(LOin=LO.in), simplify=F)
ma.fe <- do.call(rbind, ma.fe)


plan(multiprocess)
ma.all <- future_map_dfr(1:1000, simple_binary_MA, .progress=T)
plan(sequential)

#################################################
# results

ma.all %$% quantile(mu.fe,probs=c(0.025,0.5,0.975))
ma.all %$% quantile(mu.re,probs=c(0.025,0.5,0.975))

ma.all %$% quantile(lower.random,probs=c(0.025,0.5,0.975))
ma.all %$% quantile(upper.random,probs=c(0.025,0.5,0.975))

ma.all %$% quantile(lower.hakn,probs=c(0.025,0.5,0.975))
ma.all %$% quantile(upper.hakn,probs=c(0.025,0.5,0.975))


ma.fe %$% quantile(mu.fe,probs=c(0.025,0.5,0.975))
ma.fe %$% quantile(mu.re,probs=c(0.025,0.5,0.975))

ma.fe %$% quantile(lower.random,probs=c(0.025,0.5,0.975))
ma.fe %$% quantile(upper.random,probs=c(0.025,0.5,0.975))

ma.fe %$% quantile(lower.fixed,probs=c(0.025,0.5,0.975))
ma.fe %$% quantile(upper.fixed,probs=c(0.025,0.5,0.975))

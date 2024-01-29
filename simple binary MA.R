library(tidyverse)
library(magrittr)
library(furrr)
library(meta)


simple_binary_MA <- function(Nt=10,N=100,mu=1,tau2=1,LOin=NULL){

  if (is.null(LOin))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- LOin

  p.i <- exp(theta.i) / (1+exp(theta.i))
  Y <- rbinom(Nt, N, p.i)
  tsize <- rep(N, Nt)

  ma.std <- metaprop(Y, tsize, sm="PLOGIT", prediction=T)
  ma.hakn <- metaprop(Y, tsize, sm="PLOGIT", hakn=T)

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

ma.all <- replicate(10000, simple_MA(), simplify=F)
ma.all <- do.call(rbind, ma.all)

theta.in <- rnorm(10,1,1)
ma.fe <- replicate(10000, simple_MA(thetain=theta.in), simplify=F)
ma.fe <- do.call(rbind, ma.fe)


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

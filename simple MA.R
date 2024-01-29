library(tidyverse)
library(magrittr)
library(furrr)
library(meta)

setwd("~/R/basic principles")

source("very simple MA simulators.R")

#####################################################
# run evaluations

cont_func <- function(Nt,Nmin,Nmax,mu,tau2,epsilon,thetain){
  data.in <- simple_MA_data(Nt,Nmin,Nmax,mu,tau2,epsilon,thetain)
  ma.res <- run_MA(data.in)
  return(ma.res)
}

bin_func <- function(Nt,Nmin,Nmax,mu,tau2,LOin){
  data.in <- simple_binary_MA_data(Nt,Nmin,Nmax,mu,tau2,LOin)
  ma.res <- run_MA(data.in)
  return(ma.res)
}

# continuous
theta.in <- runif(20,-1,1)
theta.in <- (theta.in - mean(theta.in)) / sd(theta.in)

plan(multiprocess(workers=15))

ma.re <- future_map_dfr(1:10000, ~cont_func(Nt=20,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,thetain=NULL))
ma.fe <- future_map_dfr(1:10000, ~cont_func(Nt=20,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,thetain=theta.in))

# total homg.
theta.in <- rep(0,20)
ma.re1 <- future_map_dfr(1:10000, ~cont_func(Nt=20,Nmin=100,Nmax=1000,mu=0,tau2=0,epsilon=1,thetain=NULL))
ma.fe1 <- future_map_dfr(1:10000, ~cont_func(Nt=20,Nmin=100,Nmax=1000,mu=0,tau2=0,epsilon=1,thetain=theta.in))

# binary
LO.in <- rnorm(10,0.1,sqrt(0.05))
ma.bin.re <- future_map_dfr(1:10000, ~bin_func(Nt=10,Nmin=100,Nmax=1000,mu=0.1,tau2=0.05,LOin=NULL))
ma.bin.fe <- future_map_dfr(1:10000, ~bin_func(Nt=10,Nmin=100,Nmax=1000,mu=0.1,tau2=0.05,LOin=LO.in))

plan(sequential)


write_csv(ma.re,"cont re results 1606.csv")
write_csv(ma.fe,"cont fe results 1606.csv")

write_csv(ma.re1,"cont re homog results 1606.csv")
write_csv(ma.fe1,"cont fe homog results 1606.csv")

write_csv(ma.bin.re,"bin re results 1606.csv")
write_csv(ma.bin.fe,"bin fe results 1606.csv")


### Stan Wrapper for the piece-wise exponential model
### With individual level data
### December 8, 2021
### For Revisions of Annals of Internal Medicine manuscript


###########################################
### Libraries
###########################################

library(dplyr)
library(bayesplot)
library(gridExtra)
library(knitr)
library(reshape)
library(xtable)
library(grid)
library(scales)
library(ggplot2)
library(loo)
library("rstan")

#-- Stan settings

Sys.setenv(USE_CXX14 = 1)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


###########################################
### MODEL SETTINGS
###########################################

#-- Set the change points for onset hazard
s_vec<-c(40, 55, 70)


################################################
### LOAD DATA FILES
################################################

load(file="Simulated_cohort_07071219.RData")

dat<-dat.sim

time.tests <- dat[,-c(1:5)] # remove the first 5 columns to keep only the screening times

## Add the onset time
time.tests<-cbind(rep(s_vec[1], dim(time.tests)[1]), time.tests)

## NON-CASES with multiple follow-ups
i0 <- which(dat$case==0)
Xmtimescreen<- as.matrix(time.tests[i0,])
tXmclinical = dat$eventime[i0]
for(i in 1:length(i0)){
  ind<-which(Xmtimescreen[i,]==0);
  Xmtimescreen[i,ind[1]]=tXmclinical[i]
}


## INTERVAL CASES
i1 <- which(dat$case==1 & dat$screen==0)
Ymtimescreen<- as.matrix(time.tests[i1,])
tYmclinical = dat$eventime[i1]
for(i in 1:length(i1)){
  ind<-which(Ymtimescreen[i,]==0);
  Ymtimescreen[i, ind[1]]=tYmclinical[i]
}

## SCREEN CASES 
i2 <- which(dat$case==1 & dat$screen==1)
Zmtimescreen=as.matrix(time.tests[i2,])


data_BCSC = list("ncol"=ncol(time.tests),  "M0"=length(i0), "M1"=length(i1), "M2"=length(i2), "m0"=dat$ni[i0], "m1"=dat$ni[i1], "m2"=dat$ni[i2], 
                 "Xmtimescreen"=Xmtimescreen, "Xm"=rep(1,length(i0)),
                 "Ymtimescreen"=Ymtimescreen, "Ym"=rep(1,length(i1)),
                 "Zmtimescreen"=Zmtimescreen, "Zm"=rep(1,length(i2)), 
                 "s_vec"=s_vec)


###########################################
### RUN STAN
###########################################

#- Number of chains
n_chains=2

#- Load the Stan model
mod <- stan_model(file = "BCSC_v1.stan", verbose = TRUE)

#- Start timer
ptm <- proc.time()

#- Fit the model
fit1 <- sampling(mod, data = data_BCSC,   # named list of data
                 init = rep(list(list(h0=1/100, h1=1/100, h2=1/100, psi=.5, w=1/10, beta=0.5)),n_chains),
                 chains = n_chains,             # number of Markov chains
                 warmup = 200,          # number of warmup iterations per chain
                 iter = 2000,            # total number of iterations per chain
                 cores = n_chains,              # number of cores
                 refresh = 1,           # no progress shown if 0
                 control = list(adapt_delta = 0.8)
)

#- Stop the timer
proc.time() - ptm

#- Save the model fit and s_vec
save(fit1, s_vec, file="Stan_BCSC_output.RData")



###########################################
### Load the Stan fit from previous run
###########################################

load(file="Stan_BCSC_output.RData")


###########################################
### Output summary and analytics
###########################################

#chain tracing
traceplot(fit1,  inc_warmup = FALSE, ncol=3,pars = c("h0", "h1", "h2", "psi", "lambda", "beta"))

## pairwise correlations
pairs(fit1, pars = c("h0", "h1", "h2", "psi", "lambda", "beta"), las = 1)


### Subfigure
pairs(fit1, pars = c("psi", "lambda"), las = 1)


# summary statistics
fit_summary <- summary(fit1)
print(fit_summary$summary)
round(fit_summary$summary[c(1:6),c(1,4,6,8) ],5)


# All chains combined
sampler_params <- get_sampler_params(fit1, inc_warmup = FALSE)
summary(do.call(rbind, sampler_params), digits = 2)

# Autocorrelation
df_of_draws <- as.data.frame(fit1)
mcmc_acf(df_of_draws, lags = 10)


#######################################################
### Prior-posterior plots
#######################################################

#-- Fraction indolent 
p.1 <- ggplot(df_of_draws, aes(x=psi) ) +
  geom_density( bw=.02,size=1) +
  stat_function(fun=dbeta, 
                args=list(shape1=1, shape2= 1, ncp = 0, log = FALSE),
                colour="blue") +
  xlim(0,.3) + 
  xlab("Fraction indolent cancers")

#-- Tumor latency
p.2<-ggplot(df_of_draws, aes(x=lambda) ) +
  geom_density(bw=.01, size=1) +
  stat_function(fun=dexp, 
                args=list(rate=.01),
                colour="blue",) +
  xlab("Progression rate")+
  ylab("")


#-- Sensitivity
p.3<-ggplot(df_of_draws, aes(x=beta) ) +
  geom_density(bw=.01,size=1) +
  stat_function(fun=dbeta, 
                args=list(shape1=38.4888, shape2= 5.7512, ncp = 0, log = FALSE),
                colour="blue",) +
  xlab("Screening sensitivity")+
  xlim(0.5,1)+
  ylab("")


#-- Onset h0
p.4<-ggplot(df_of_draws, aes(x=h0) ) +
  geom_density(size=1) +
  stat_function(fun=dexp, 
                args=list(rate=.01),
                colour="blue",) +
  xlab("Onset rate h0: 40-54y")

#-- Onset h1
p.5<-ggplot(df_of_draws, aes(x=h1) ) +
  geom_density(size=1) +
  stat_function(fun=dexp, 
                args=list(rate=.01),
                colour="blue",) +
  xlab("Onset rate h1: 55-64y")+
  ylab("")

#-- Onset h2
p.6<-ggplot(df_of_draws, aes(x=h2) ) +
  geom_density(size=1) +
  stat_function(fun=dexp, 
                args=list(rate=.01),
                colour="blue",) +
  xlab("Onset rate h2: 65y+")+
  ylab("")



multiplot(p.1, p.4, p.2, p.5,p.3, p.6, cols=3)


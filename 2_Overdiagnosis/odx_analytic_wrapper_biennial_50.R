######################################################
# Analytic overdiagnosis calculations
# Version December 22, 2021
# For revisions to Annals of Internal Medicine
######################################################


library(rstan)
library(ggridges)
library(tidyverse)
library(forcats)

#-- Load the auxiliary functions
source("odx_analytic_aux_functions.R")

#-- Load the data
load(file="Stan_BCSC_A_50_to_74_onset_40_55_65.RData")
df_of_draws <- as.data.frame(fit1)


#---------------------------------------------
### CHOICES
#---------------------------------------------

# Age at first screen
age0=50

# Number of MC draws
M=4000

# not a choice: number of biennial screens in total
ll=(74-age0)/2+1 


#---------------------------------------------
### OTHER CAUSE MORTALITY DATA
#---------------------------------------------

# #-- The survival functions and hazard functions
# cohort=max(2021-age0, 1950)
# 
# # mac version
# OCdata <- read.csv("MortalityRatesTables_17Mar2017.csv", header=TRUE)
# 
# OChaz<-subset(OCdata, Cohort==cohort)
# OChaz$Estimate<-OChaz$Estimate/10^5
# 
# save(OChaz, file="OtherCause_birthcohort_1971.RData")

load(file="OtherCause_birthcohort_1971.RData")

#---------------------------------------------
#---------------------------------------------
### PROBABILITIES OF SCREEN-DETECTION
#---------------------------------------------
#---------------------------------------------

screen.ind<-array(0,dim=c(ll,M))
screen.prog<-array(0,dim=c(ll,M))
screen.total<-array(0,dim=c(ll,M))

lambda.seq<-rep(0,M)
JART<-sample(dim(df_of_draws)[1],M, replace = FALSE)

zahler=1
for(TTT in JART){
  
  x.vec<-df_of_draws[TTT, 1:6]
  
  for(tt in (0:(ll-1)*2)){
    
    t.vec<-c(s_vec[1], seq(from=age0, to=(age0+tt), by =2))

    # number of screens
    ni=length(t.vec)-1
    
    # P(SD)
    d<-as.numeric(D_j(x.vec, t.vec, s_vec,  ni))
    # P(SD, I=1)
    d.ind<-as.numeric(D_j_IND(x.vec, t.vec, s_vec,  ni))
    # P(SD, I=0)
    d.prog<-as.numeric(D_j_PROG(x.vec, t.vec, s_vec,  ni))
    
    ## total SD probability
    screen.total[ni, zahler]<-d
    
    ## indolent SD probability
    screen.ind[ni, zahler]<-d.ind
    
    ## progressive SD probability
    screen.prog[ni,zahler]<-d.prog # *Prog.odx(age0+tt, x.vec$lambda, OChaz$Estimate)
  
    ## The lambda used
    lambda.seq[zahler]<-x.vec$lambda
  }
  zahler=zahler+1
  print(round(100*zahler/M, 1))
} 


#---------------------------------------------
#---------------------------------------------
# ADD IN OTHER CAUSE DEATH
#---------------------------------------------
#---------------------------------------------


#-- Screening sequence
tt<-seq(age0, 74, by=2)
#-- P(screen-detected cancer)
SD.stat<-rep(0,M)
#-- P(indolent screen-detected)
ODX.ind.stat<-rep(0, M)
#-- P(progressive overdiagnosed cancer)
ODX.prog.stat<-rep(0,M)


#--Probability to have a screen-detected cancer
for(z in 1:dim(screen.total)[2]){
  
  # probability to survive past each respective screening round
  hh=apply(t(tt),2, surv.fun, age0=age0, OChaz=OChaz$Estimate)
  # probability of a screen-detected cancer while still alive
  SD.stat[z]=sum(screen.total[,z]*hh)
  
}

## Probability to have an indolent screen-detected cancer

for(z in 1:dim(screen.total)[2]){
  
  # probability to survive past each respective screening round
  hh=apply(t(tt),2, surv.fun, age0=age0, OChaz=OChaz$Estimate)
  # probability of an indolent screen-detected cancer while still alive
  ODX.ind.stat[z]=sum(screen.ind[,z]*hh)
  
}

## Probability to have a progressive and overdiagnosed screen-detected cancer

for(z in 1:dim(screen.total)[2]){
  
  hh=apply(t(tt),2, Prog.odx.program, age0=age0, lambda=lambda.seq[z], OChaz=OChaz$Estimate)
  ODX.prog.stat[z]=sum(screen.prog[,z]*hh)
  
  print(z)
}

save(lambda.seq,
     screen.prog,
     screen.ind,
     screen.total,
     SD.stat,
     ODX.ind.stat,
     ODX.prog.stat,
     file=paste("OverDX_USPSTF_", age0, "_74_model_A_biennial.RData", sep=""))



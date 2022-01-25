######################################
## Biennial: starting age 50
######################################

########### Wrapper for the summary statistics

source("multiplot.R")
library(ggplot2)

wrap <- function(ODX.ind.stat, ODX.prog.stat, SD.stat){
  
  
  
  #overall
  a<-c(round(100*mean((ODX.ind.stat+ODX.prog.stat)/SD.stat),1),
       round(100*quantile((ODX.ind.stat+ODX.prog.stat)/SD.stat, probs=c(.025, .5, .975)),1))
  #indolent
  b<-c(round(100*mean((ODX.ind.stat)/SD.stat),1),
       round(100*quantile((ODX.ind.stat)/SD.stat, probs=c(.025, .5, .975)),1))
  #mortality
  c<-c(round(100*mean((ODX.prog.stat)/SD.stat),1),
       round(100*quantile((ODX.prog.stat)/SD.stat, probs=c(.025, .5, .975)),1))
  
  
  tab<-rbind(a,b,c)
  colnames(tab)[1]<-"mean"
  tab<-tab[,c(1,2,4,3)]
  rownames(tab)<-c("total", "indolent", "mortality")
  
  return(tab)
  
}

load(file="OverDX_USPSTF_50_74_model_A_biennial_r3.RData")
hist((ODX.ind.stat+ODX.prog.stat)/SD.stat, breaks=50)

#overall
tab<-wrap(ODX.ind.stat, ODX.prog.stat, SD.stat)
print(tab)


### Calculating over-diagnosis contributions


## Onset hazard function
h <- function(t, s, theta){

    out = 0;
    h0 = theta[1]; h1 = theta[2]; h2 = theta[3];
    s0 = s[1]; s1 = s[2]; s2 = s[3];
    
    if (t<s0){
      out =0;
    } else if (t>=s0 & t < s1){
      out = h0;
    } else if(t>= s1 & t<s2){
      out = h1;
    }    else {
      out = h2;
    } 
    
    return(out);
  }  


## Onset cumulative hazard function

Lambda <- function(t, s, theta){
  
    out = 0;
    h0 = theta[1]; h1 = theta[2]; h2 = theta[3];
    s0 = s[1]; s1 = s[2]; s2 = s[3];
    
    if (t<s0){
      out =0;
    } else if (t>=s0 & t < s1){
      out = h0*(t-s0);
    } else if(t>= s1 & t<s2){
      out = h0*(s1-s0)+h1*(t-s1);
    } else {
      out = h0*(s1-s0)+h1*(s2-s1)+h2*(t-s2);
    } 
    
    return(out);
    
  }

## Probability density function of preclinical onset

f <- function(t, s, theta){
  
    out = 0;
    h0 = theta[1]; h1 = theta[2]; h2 = theta[3];
    s0 = s[1]; s1 = s[2]; s2 = s[3];
    
    if (t<s0){
      out = 0;
    } else if (t>=s0 & t < s1){
      out = h0*exp(-h0*(t-s0));
    } else if(t>= s1 & t<s2){ 
      out = h1*exp(-h0*(s1-s0)-h1*(t-s1));
    } else {
      out = h2*exp(-h0*(s1-s0)-h1*(s2-s1)-h2*(t-s2));
    } 
    
    return(out);
    
  }

## Probability function of preclinical onset: P(T ≤= t)

F <- function(t, s, theta){
  
    out = 0;
    h0 = theta[1]; h1 = theta[2]; h2 = theta[3];
    s0 = s[1]; s1 = s[2]; s2 = s[3];
    
    if (t<s0) {
      out = 0;
    } else if (t>=s0 & t < s1){ 
      out = 1 - exp(-h0*(t-s0));
    } else if(t>= s1 & t<s2){
      out = 1 - exp(-h0*(s1-s0)-h1*(t-s1));
    } else {
      out = 1 - exp(-h0*(s1-s0)-h1*(s2-s1)-h2*(t-s2));
    } 
    
    return(out);
    
  }
  
##  Auxiliary function R
  
R <- function(L, U, h, s, lambda, t){
    
    out = h/(lambda-h)*(exp(U*(lambda-h))-exp(L*(lambda-h)))*exp(h*s-lambda*t);
    
    return(out);
    
  }
  
##  Auxiliary function Omega
Omega <- function(ta, tb, tc, s, theta){
    
    out=0;
    h0 = theta[1]; h1 = theta[2]; h2 = theta[3];
    s0 = s[1]; s1 = s[2]; s2 = s[3];
    lambda = theta[5];
    
    if(ta<s1 & tb>s0){
      out = out + R(max(s0,ta), min(s1, tb), h0, s0, lambda, tc);
    } 
      
    if(ta<s2 & tb>s1){
        out = out+ exp(-h0*(s1-s0))*R(max(s1, ta), min(s2, tb), h1, s1, lambda, tc);
    } 
    
    if(tb>s2){ 
      out = out+ exp(-h0*(s1-s0)-h1*(s2-s1))*R(max(s2, ta), tb, h2, s2, lambda, tc);
    } 
          
    return(out);
    
  }
  

## ALL SCREEN-DETECTED CASES likelihood contribution
  
D_j <- function(theta, t, s,  ni){
    
    l=ni+1;
    psi = theta[4];  beta = theta[6];
    A = 0; B=0;
    out = 0;
      
## sum over previous screens 
    
      for(k in 2:l){
        
        A = A +  psi*beta*'^'(1-beta, l-k)*(F(t[k], s, theta)-F(t[k-1], s, theta));
        
        B = B + (1-psi)*beta*'^'(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
        
      }
      
      out = A + B;
      
      out= out/ (psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta)));
      
      return(out);
}  


## INDOLENT SCREEN-DETECTED CASES likelihood: P(IND=1, SD)= P(SD | IND=1)P(IND=1)

D_j_IND <- function(theta, t, s,  ni){
  
  l=ni+1;
  psi = theta[4];  beta = theta[6];
  A = 0; B=0;
  out = 0;
  
  ## sum over previous screens 
  
  for(k in 2:l){
    
    # ONLY THE INDOLENT FRACTION
    A = A +  psi*beta*'^'(1-beta, l-k)*(F(t[k], s, theta)-F(t[k-1], s, theta)); 
    
    ## B = B + (1-psi)*beta*'^'(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
    
  }
  
  out = A + B;
  
  out= out / (psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta))); 
  
  return(out);
} 
  

## PROGRESSIVE SCREEN-DETECTED CASES likelihood:  P(IND=0, SD)= P(SD | IND=0)P(IND=0)

D_j_PROG <- function(theta, t, s,  ni){
  
  l=ni+1;
  psi = theta[4];  beta = theta[6];
  A = 0; B=0;
  out = 0;
  
  ## sum over previous screens 
  
  for(k in 2:l){
    
    # A = A +  psi*beta*'^'(1-beta, l-k)*(F(t[k], s, theta)-F(t[k-1], s, theta));
    
    # ONLY THE PROGRESSIVE CONTRIBUTION
    B = B + (1-psi)*beta*'^'(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
    
  }
  
  out = A + B;
  
  out= out/ (psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta)));
  
  return(out);
}  




###############################################
### Other cause death functions  ##############
###############################################

## hazard for other cause mortality

lambda_oc<-function(t, OChaz){
  
  ft<- floor(t)
  out=ifelse(ft<=119, OChaz[ft+1],OChaz[120])
  
  return(out)
  
}

lambda_oc<-Vectorize(lambda_oc, "t")


## other cause mortality survival function 

surv.fun <- function(age0, t, OChaz){
  
  a<-age0
  b<-floor(t)-1
  
  if(b>=a){ 
    
    v=apply(t(a:b),1,lambda_oc, OChaz=OChaz)
    
  } else {
    
    v=0
    
  }
  
  hh<-sum(v) + (t-floor(t))*lambda_oc(t, OChaz=OChaz)
  
  S<-exp(-hh)
  
  return(S)
}

surv.fun<-Vectorize(surv.fun,vectorize.args = "t")

######### P(TM<U)=1-P(TM>U)


Prog.odx <- function(age0, lambda, OChaz){
  
  help.f<-function(s){
    
    return(surv.fun(age0,age0+s, OChaz)*lambda*exp(-lambda*(s)))
    
  }
  out<- integrate(help.f, lower=0, upper=1000)
  
  return(1-out$value)
}


########### Prog.odx.programmatic


Prog.odx.program <- function(age0, t.k, lambda, OChaz){
  
  help.f<-function(t){
    
    return(surv.fun(age0,t.k+t, OChaz)*lambda*exp(-lambda*(t)))
    
  }
  out<- integrate(help.f, lower=0, upper=1000)
  
  return(surv.fun(age0,t.k, OChaz)-out$value)
}

Prog.odx.program<-Vectorize(Prog.odx.program,vectorize.args = "t.k")



########### Wrapper for the summary statistics

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



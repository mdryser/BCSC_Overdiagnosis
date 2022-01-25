functions {
  
  //----------------------------------
  // hazard function h
  //----------------------------------
  real h(real t, real[] s, real[] theta){
    
    real out = 0;
    real h0 = theta[1]; real h1 = theta[2]; real h2 = theta[3];
    real s0 = s[1]; real s1 = s[2]; real s2 = s[3];
    
    if (t<s0)
    out =0;
    else if (t>=s0 && t < s1)
    out = h0;
    else if(t>= s1 && t<s2)
    out = h1;
    else 
    out = h2;
    
    return(out);
  }  
  
  //----------------------------------
  // cumulative hazard function Lambda
  //----------------------------------
  real Lambda(real t, real[] s, real[] theta){
    
    real out = 0;
    real h0 = theta[1]; real h1 = theta[2]; real h2 = theta[3];
    real s0 = s[1]; real s1 = s[2]; real s2 = s[3];
    
    if (t<s0)
    out =0;
    else if (t>=s0 && t < s1)
    out = h0*(t-s0);
    else if(t>= s1 && t<s2)
    out = h0*(s1-s0)+h1*(t-s1);
    else 
    out = h0*(s1-s0)+h1*(s2-s1)+h2*(t-s2);
    
    return(out);
    
  }
  
  //---------------------------------------------------
  // PDF of preclinical onset
  //---------------------------------------------------
  real f(real t, real[] s, real[] theta){
    
    real out = 0;
    real h0 = theta[1]; real h1 = theta[2]; real h2 = theta[3];
    real s0 = s[1]; real s1 = s[2]; real s2 = s[3];
    
    if (t<s0)
    out = 0;
    else if (t>=s0 && t < s1)
    out = h0*exp(-h0*(t-s0));
    else if(t>= s1 && t<s2)
    out = h1*exp(-h0*(s1-s0)-h1*(t-s1));
    else 
    out = h2*exp(-h0*(s1-s0)-h1*(s2-s1)-h2*(t-s2));
    
    return(out);
    
  }
  
  //----------------------------------
  // CDF of preclinical onset: P(T ≤= t)
  //----------------------------------
  
  real F(real t, real[] s, real[] theta){
    
    real out = 0;
    real h0 = theta[1]; real h1 = theta[2]; real h2 = theta[3];
    real s0 = s[1]; real s1 = s[2]; real s2 = s[3];
    
    if (t<s0)
    out = 0;
    else if (t>=s0 && t < s1)
    out = 1 - exp(-h0*(t-s0));
    else if(t>= s1 && t<s2)
    out = 1 - exp(-h0*(s1-s0)-h1*(t-s1));
    else 
    out = 1 - exp(-h0*(s1-s0)-h1*(s2-s1)-h2*(t-s2));
    
    return(out);
    
  }
  
  //----------------------------------
  // Auxiliary function: R
  //----------------------------------
  real R(real L, real U, real h, real s, real lambda, real t){
    
    real out = h/(lambda-h)*(exp(U*(lambda-h))-exp(L*(lambda-h)))*exp(h*s-lambda*t);
    
    return(out);
    
  }
  
  //----------------------------------
  // Auxiliary function: Omega
  //----------------------------------
  real Omega(real ta, real tb, real tc, real[] s, real[] theta){
    
    real out=0;
    real h0 = theta[1]; real h1 = theta[2]; real h2 = theta[3];
    real s0 = s[1]; real s1 = s[2]; real s2 = s[3];
    real lambda = theta[5];
    
    if(ta<s1 && tb>s0)
    out += R(fmax(s0,ta), fmin(s1, tb), h0, s0, lambda, tc);
    
    if(ta<s2 && tb>s1)
    out += exp(-h0*(s1-s0))*R(fmax(s1, ta), fmin(s2, tb), h1, s1, lambda, tc);
    
    if(tb>s2)
    out += exp(-h0*(s1-s0)-h1*(s2-s1))*R(fmax(s2, ta), tb, h2, s2, lambda, tc);
    
    return(out);
  }
  
  //-------------------------------------
  // Likelihood contribution: NON-CASES
  //-------------------------------------
  real L_j(real[] theta, real[] t, real[] s, int ni){
    
    int l=ni+2; // number of elements in the t vector
    real psi = theta[4];  real beta = theta[6];
    real A = 0; real B=0; real C;
    real out =0;
    
    
    for(k in 2:l){
      
      A += psi*pow(1-beta,l-k)*(F(t[k], s, theta)-F(t[k-1], s, theta));
      
      B += (1-psi)*pow(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
      
    }
    
    C = 1-F(t[l], s, theta);
    
    out = A + B + C; 
    
    // normalization
    out /= psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta)); 
    
    return(out);
  }
  
  
  //-------------------------------------------
  // Likelihood contribution: SCREEN-DETECTED
  //-------------------------------------------
  real D_j(real[] theta, real[] t, real[] s, int ni){
    
    int l=ni+1;    //
    real psi = theta[4];  real beta = theta[6];
    real A = 0; real B=0;
    real out =0;
    
    
    // sum over previous screens 
    for(k in 2:l){
      
      A += psi*beta*pow(1-beta, l-k)*(F(t[k], s, theta)-F(t[k-1], s, theta));
      
      B += (1-psi)*beta*pow(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
      
    }
    
    out = A + B;
    
    // normalization
    out /= psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta)); 
    
    return(out);
  }  
  
  
  //-------------------------------------------
  // Likelihood contribution: CLINICAL CASES
  //-------------------------------------------
  real I_j(real[] theta, real[] t, real[] s, int ni){ 
    
    int l=ni+2; 
    real psi = theta[4];  real beta = theta[6];  real lambda = theta[5];
    
    real A = 0;
    real out =0;
    
    
    for(k in 2:l){
      
      A += pow(1-beta, l-k)*Omega(t[k-1], t[k], t[l], s, theta);
      
    }
    
    
    out = lambda*(1-psi)*A;
    
    // normalization
    out /= psi+(1-psi)*(1-F(t[2], s, theta)+Omega(t[1], t[2], t[2], s, theta)); 
    
    
    return(out);
    
  }
  
  
}
data {
  
  int ncol;
  
  int M0;  // Number of non-cases
  int M1;  // Number of interval cases
  int M2;  // Number of screen cases
  
  int m0[M0]; // number of screens
  int m1[M1]; // 
  int m2[M2]; //
  
  //real tXmclinical[M0]; // censoring time
  real Xmtimescreen[M0, ncol]; // the screening times for non-cases
  int Xm[M0];  // dummy for sampler
  
  // real tYmclinical[M1]; // event time of the interval case
  real Ymtimescreen[M1, ncol]; // screening times  for interval cases
  int Ym[M1];  // dummy
  
  // real tZmclinical[M2];  // event time (screen time)
  real Zmtimescreen[M2, ncol]; // screening times for screen cases
  int Zm[M2]; // dummy
  
  real s_vec[3];
  
}
parameters {
  // Piece-wise constant onset hazards
  real<lower=0> h0;
  real<lower=0> h1;
  real<lower=0> h2;
  
  // Fraction indolent
  real<lower=0, upper=1> psi;
  // Exponential rate for tumor latency
  real<lower=0> lambda; 
  //  Screening round sensitivity
  real<lower=0, upper=1> beta;           
  
  
}
transformed parameters {
  
  real theta[6];
  
  theta[1]=h0;
  theta[2]=h1;
  theta[3]=h2;
  theta[4]=psi;
  theta[5]=lambda;
  theta[6]=beta;
  
  
}
model { 
  
  // auxiliary variables for the probabilities
  real pm[M0];  
  real qm[M1];
  real rm[M2];
  
  
  
  //---------------------------------------------
  // Likelihood for NON-CASES
  //---------------------------------------------
  
  for (i in 1:M0){ 
    
    pm[i]=L_j(theta, Xmtimescreen[i,], s_vec, m0[i]);
    
  }
  
  // Update target density
  target+=bernoulli_lpmf(Xm | pm);
  
  // Make sure probabilty is properly normalized for the Bernoulli trick
  if(max(pm)>1){
    print(pm);
  }
  
  
  //---------------------------------------------
  // Likelihood for CLINICAL CASES
  //---------------------------------------------
  
  
  for (i in 1:M1){ 
    
    qm[i]=I_j(theta, Ymtimescreen[i,], s_vec, m1[i]);
    
  }
  
  // Update target density
  target+=bernoulli_lpmf(Ym | qm);
  
  // Make sure probabilty is properly normalized for the Bernoulli trick
  if(max(qm)>1){
    print(qm);
  }
  
  
  //---------------------------------------------
  // Likelihood for SCREEN-DETECTED CASES
  //---------------------------------------------
  
  for (i in 1:M2){ 
    
    rm[i]=D_j(theta, Zmtimescreen[i,], s_vec, m2[i]);
    
  }
  
  // Update target density
  target+=bernoulli_lpmf(Zm | rm);
  
  // Make sure probabilty is properly normalized for the Bernoulli trick
  if(max(rm)>1){
    print(rm);
  }
  
  
  //---------------------------------------------
  // PRIOR DISTRIBUTIONS
  //---------------------------------------------  
  
  beta ~beta(38.4888,5.7512); // screening sensitivity
  h0 ~ exponential(.01);      // onset hazards
  h1 ~ exponential(.01);
  h2 ~ exponential(.01);
  lambda ~ exponential(.01);  // tumor latency rate
  psi ~ beta(1,1);            // fraction indolent
  
}

generated quantities {
  
  // Keep track of the log likelihood
  real LKH=0;
  
  for (i in 1:M0){ 
    LKH+= bernoulli_lpmf(Xm[i] | L_j(theta, Xmtimescreen[i,], s_vec, m0[i]));
  }
  
  for (i in 1:M1){ 
    LKH+= bernoulli_lpmf(Ym[i] | I_j(theta, Ymtimescreen[i,], s_vec, m1[i]));
  }  
  
  for (i in 1:M2){ 
    LKH+= bernoulli_lpmf(Zm[i] | D_j(theta, Zmtimescreen[i,], s_vec, m2[i]));
  }    
  
  
}

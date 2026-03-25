data {
  int<lower=0> N;           //number of years
  vector[N] surv;           //logit SAS for recruits
  vector[N] lrs;            //log recruits per spawner
  vector[N] S;              //log recruits per spawner
  real mlogit_surv;         //mean logit-transformed smolt-to-adult survival
  int<lower=0> faN;         //number of Harrison fecundity measures
  int<lower=0> suN;         //number of Lower Shuswap fecundity measures
  vector[faN] ffa;          //mean Fraser Fall 0.3 Chinook fecundity
  vector[suN] fsu;          //mean Fraser Summer 0.3 Chinook fecundity
}
parameters {
  real<lower=0> Smax;          //capacity - spawners that maximizes recruitment
  real<lower=0> alpha;         //Harrison Ricker alpha parameter
  real gamma;                  //Harrison survival index parameter
  real<lower=0> sigma;         //Harrison SD within the autocorrelated process
  real<lower=-1,upper=1> rho;  //Harrison autocorrelation parameter
  real<lower=0> f_fa;          //fecundity of Fraser Fall 0.3 Chinook
  real<lower=0> f_su;          //fecundity of Fraser Summer 0.3 Chinook
  real<lower=0> sdfa;
  real<lower=0> sdsu;
}
transformed parameters {
real<lower=0> beta=1.0/Smax;   //beta - per capita density dependence parameter
vector[N] mu;                  //expectation in each year
vector[N] epsilon;             //log(R/S) residuals
real<lower=0> sigma_AR;        //sigma - corrected for autocorrelation

mu[1:N] = alpha-beta * S[1:N] + gamma * surv[1:N];//subsequent expectation with estimated survival
epsilon[1] = lrs[1] - mu[1];
for(t in 2:N){
    epsilon[t] =(lrs[t] - mu[t]);
    mu[t] = mu[t] + (rho*epsilon[t-1]); 
  }
  sigma_AR = sigma*sqrt(1-rho^2);                 //sigma corrected for autocorrelation parameter rho
}
model {
  lrs[1] ~ normal(mu[1], sigma);
  for(i in 2:N)lrs[i] ~ normal(mu[i], sigma_AR);
  Smax ~ normal(50000,25000); 
  alpha ~ normal(1,1);
  gamma ~ normal(0, 1);
  sigma ~ cauchy(0, 2);
  f_fa ~ normal(ffa, sdfa);
  f_su ~ normal(fsu, sdsu);
  sdfa ~ cauchy(0, 1000);
  sdsu ~ cauchy(0, 1000);
  rho ~ uniform(-1,1);
}
generated quantities {
  vector[N] nu_Y;
  vector[N] nu_rec;
  real alpha_fa;
  real umsy_fa;
  real umsy_su;
  real p_f_fa_su;
  
  nu_Y[1] = normal_rng(mu[1], sigma);
  nu_rec[1] = exp(nu_Y[1]) * S[1];
  for(i in 2:N){
      nu_Y[i] = normal_rng(mu[i], sigma_AR);
      nu_rec[i] = exp(nu_Y[i]) * S[i];
  }
  
  #ln(alpha) posteriors for Spring and Summer 1.3 SMUs
  alpha_fa = log(exp(alpha) * (f_fa / f_su));
  umsy_fa = 1 - lambert_w0(exp(1 - (alpha_fa)));
  umsy_su = 1 - lambert_w0(exp(1 - (alpha)));
  p_f_fa_su = f_fa / f_su;
}


data {
  int<lower=0> N;           //number of years
  vector[N] surv;           //logit SAS for recruits
  vector[N] lrs;            //log recruits per spawner
  vector[N] S;              //log recruits per spawner
  real mlogit_surv;         //mean logit-transformed smolt-to-adult survival
  int<lower=0> suN;
  int<lower=0> spN;
  int<lower=0> nicN;
  vector[suN] fsu;                //mean Fraser Summer 1.3 Chinook fecundity
  vector[spN] fsp;                //mean Fraser Spring 1.3 Chinook fecundity
  vector[nicN] fnic;               //mean Fraser Spring 1.2 (Nicola) Chinook fecundity
}
parameters {
  real<lower=0> Smax;       //capacity - spawners that maximizes recruitment
  real<lower=0> alpha;      //population-specific Ricker alpha parameter
  real gamma;               //population-specific survival index parameter
  real<lower=0> sigma;      //population-specific SD within the autocorrelated process
  real surv_est;            //estimating 1992 survival as the hatchery smolts released from the 1992 brood had a mixobacterial infection that caused high mortality rates
  real<lower=-1,upper=1> rho; //autocorrelation parameter
  real<lower=0> f_su;                  //fecundity of Fraser Summer 1.3 Chinook
  real<lower=0> f_sp;                  //fecundity of Fraser Spring 1.3 Chinook
  real<lower=0> f_nic;                 //fecundity of Fraser Spring (Nicola) 1.2 Chinook
  real<lower=0> sdsu;
  real<lower=0> sdsp;
  real<lower=0> sdnic;
}
transformed parameters {
real<lower=0> beta=1.0/Smax; //beta - per capita density dependence parameter
vector[N] mu; //expectation in each year
vector[N] epsilon; //log(R/S) residuals
real<lower=0> sigma_AR; //sigma - corrected for autocorrelation

mu[1] = alpha-beta * S[1] + gamma * surv_est; //first year with unknown survival
mu[2:N] = alpha-beta * S[2:N] + gamma * surv[2:N]; //subsequent expectation with estimated survival
epsilon[1] = lrs[1] - mu[1];
for(t in 2:N){
    epsilon[t] =(lrs[t] - mu[t]);
    mu[t] = mu[t] + (rho*epsilon[t-1]); 
  }
  sigma_AR = sigma*sqrt(1-rho^2); //sigma corrected for autocorrelation parameter rho
}
model {
  lrs[1] ~ normal(mu[1], sigma);
  for(i in 2:N)lrs[i] ~ normal(mu[i], sigma_AR);
  Smax ~ normal(10000,4000); 
  alpha ~ normal(1,1);
  gamma ~ normal(0, 1);
  surv_est ~ normal(0, 1);
  sigma ~ cauchy(0, 2);
  f_su ~ normal(fsu, sdsu);
  f_sp ~ normal(fsp, sdsp);
  f_nic ~ normal(fnic, sdnic);
  sdsu ~ cauchy(0, 1000);
  sdsp ~ cauchy(0, 1000);
  sdnic ~ cauchy(0, 1000);
  rho ~ uniform(-1,1);
}
generated quantities {
  vector[N] nu_Y;
  vector[N] nu_rec;
  real alpha_su;
  real alpha_sp;
  real alpha_prime_nic;
  real alpha_prime_su;
  real alpha_prime_sp;
  real umsy_su;
  real umsy_sp;
  real umsy_nic;
  nu_Y[1] = normal_rng(mu[1], sigma);
  nu_rec[1] = exp(nu_Y[1]) * S[1];
  for(i in 2:N){
      nu_Y[i] = normal_rng(mu[i], sigma_AR);
      nu_rec[i] = exp(nu_Y[i]) * S[i];
  }
  #ln(alpha) posteriors for Spring and Summer 1.3 SMUs
  alpha_su = log(exp(alpha) * (f_su / f_nic));
  alpha_sp = log(exp(alpha) * (f_sp / f_nic));
  alpha_prime_nic = alpha + (sigma_AR^2/2);
  alpha_prime_su = log(exp(alpha_prime_nic) * (f_su / f_nic));
  alpha_prime_sp = log(exp(alpha_prime_nic) * (f_sp / f_nic));
  umsy_su = 1 - lambert_w0(exp(1 - (alpha_su)));
  umsy_sp = 1 - lambert_w0(exp(1 - (alpha_sp)));
  umsy_nic = 1 - lambert_w0(exp(1 - (alpha)));
}


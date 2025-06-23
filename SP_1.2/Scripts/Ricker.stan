// Model with no prior on carrying capacity/srep
data {
  int<lower=0> N;           //number of years
  vector[N] surv;           //logit SAS for recruits
  vector[N] lrs;            //log recruits per spawner
  vector[N] S;              //log recruits per spawner
  real mlogit_surv;         //mean logit-transformed smolt-to-adult survival
  }
parameters {
  real<lower=0> Smax;       //capacity - spawners that maximizes recruitment
  real<lower=0> alpha;      //population-specific Ricker alpha parameter
  real gamma;               //population-specific survival index parameter
  real<lower=0> sigma;      //population-specific SD within the autocorrelated process
  real surv_est;            //estimating 1992 survival as the hatchery smolts released from the 1992 brood had a mixobacterial infection that caused high mortality rates
  real<lower=-1,upper=1> rho; //autocorrelation parameter
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
  
  Smax ~ normal(5e3,1e4); //change this to something
  alpha ~ normal(1,1);
  gamma ~ normal(0, 1);
  surv_est ~ normal(0, 1);
  sigma ~ cauchy(0, 2);

  //autocorrelation term
  rho ~ uniform(-1,1);

}

generated quantities {
  vector[N] nu_Y;
  vector[N] nu_rec;
  real alpha_prime;
  real srep;
  real smsy_85;
  real smsy;
  real umsy;
  real srep_prime;
  
  nu_Y[1] = normal_rng(mu[1], sigma);
  nu_rec[1] = exp(nu_Y[1])*S[1]; 
  for(i in 2:N){
      nu_Y[i] = normal_rng(mu[i], sigma_AR);
      nu_rec[i] = exp(nu_Y[i]*S[i]);
  }
  
  alpha_prime = alpha + (sigma_AR^2/2);
  srep = (alpha)/(beta);
  smsy_85 = 0.85*((1 - lambert_w0(exp(1 - (alpha)))) / (beta));
  smsy = ((1 - lambert_w0(exp(1 - (alpha)))) / (beta));
  umsy = 1 - lambert_w0(exp(1 - (alpha)));
  srep_prime = (alpha_prime)/(beta);
}


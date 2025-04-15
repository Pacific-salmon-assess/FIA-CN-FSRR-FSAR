// Model with no prior on carrying capacity/srep
data {
  int<lower=0> N;           //number of years
  vector[N] surv;           //logit SAS for recruits
  vector[N] lrs;            //log recruits per spawner
  vector[N] S;              //log recruits per spawner
  real mlogit_surv;         //mean logit-transformed smolt-to-adult survival
  }

parameters {
  real beta;                //population-specific Ricker beta parameter
  real<lower=0> alpha;      //population-specific Ricker alpha parameter
  real gamma;               //population-specific survival index parameter
  real<lower=0> sigma;      //population-specific SD within the autocorrelated process
}

model {
    for(i in 1:N){
      lrs[i] ~ normal(beta * S[i] + gamma * surv[i] + alpha, sigma);
}
  beta ~ normal(0, 10);
  alpha ~ cauchy(0, 5);
  gamma ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
}

generated quantities {
  vector[N] nu_Y;
  vector[N] mu_gamma_Y;
  vector[N] nu_rec;
  real srep;
  real smsy_85;
  real smsy;
  real umsy;
    for(i in 1:N){
      nu_Y[i] = normal_rng(beta * S[i] + gamma * surv[i] + alpha, sigma);
      mu_gamma_Y[i] = normal_rng(beta * S[i] + gamma * mlogit_surv + alpha, sigma);
      nu_rec[i] = exp(mu_gamma_Y[i])*S[i];
  }
  srep = (alpha + gamma*mlogit_surv)/(-1*beta);
  smsy_85 = 0.85*((1 - lambert_w0(exp(1 - (alpha + gamma*mlogit_surv)))) / (-1*beta));
  smsy = ((1 - lambert_w0(exp(1 - (alpha + gamma*mlogit_surv)))) / (-1*beta));
  umsy = 1 - lambert_w0(exp(1 - (alpha + gamma*mlogit_surv)));
}


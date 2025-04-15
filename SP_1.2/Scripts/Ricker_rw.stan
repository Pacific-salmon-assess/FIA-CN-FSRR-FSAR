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
  real<lower=0> sigma;      //population-specific SD within the autocorrelated process
}

transformed parameters {
  vector[N] rw;
  for(i in 1:N){
    rw[i] = beta * S[i] + alpha;
  }
}

model {
  lrs[1] ~ normal(beta * S[1] + alpha, sigma);
    for(i in 2:N){
      lrs[i] ~ normal(beta * S[i] + alpha + rw[i-1], sigma);
}
  beta ~ normal(0, 10);
  alpha ~ cauchy(0, 5);
  sigma ~ cauchy(0, 5);
}

generated quantities {
  vector[N] nu_Y;
  vector[N] nu_rec;
  real srep;
  real smsy_85;
  real smsy;
  real umsy;
  nu_Y[1] = normal_rng(beta * S[1] + alpha, sigma);
  for (i in 2:N){
    nu_Y[i] = normal_rng(beta * S[i] + alpha, sigma);
  }
  for(i in 1:N){
    nu_rec[i] = exp(nu_Y[i])*S[i];
  }
  srep = alpha/(-1*beta);
  smsy_85 = 0.85*((1 - lambert_w0(exp(1 - (alpha)))) / (-1*beta));
  smsy = ((1 - lambert_w0(exp(1 - (alpha)))) / (-1*beta));
  umsy = 1 - lambert_w0(exp(1 - (alpha)));
}


// Model with no prior on smax
data {
  int<lower=0> N;          //number of watersheds in original dataset
  vector[N] smax;          //smax estimated for original watersheds using alpha prime and beta from Ricker SR models
  vector[N] wa;            //corresponding watershed areas for original watersheds
  int<lower=0> K;          //number of new watersheds to predict smax for
  vector[K] nu_wa;       //vector of watershed areas from new watersheds for which we are trying to estimate beta parameters from a Ricker SR model
}

parameters {
  real slope;       //slope of relationship between smax and watershed area
  real intercept;   //intercept of relationship between smax and watershed area
  real<lower=0> sigma; //error
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N){
    mu[i] = slope * wa[i] + intercept;
  }
}

model {
  for(i in 1:N){
    smax[i] ~ normal(mu[i], sigma);
  }
  
  slope ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma ~ cauchy(0, 5);
}

generated quantities {
  vector[N] nu_Y;
  vector[K] nu_smax;
  
  for(i in 1:N){
    nu_Y[i] = normal_rng(mu[i], sigma);
  }
  for(k in 1:K){
    nu_smax[k] = slope * nu_wa[k] + intercept;
  }
}


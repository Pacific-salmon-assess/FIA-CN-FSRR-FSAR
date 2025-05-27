// Model with no prior on carrying capacity/srep
data {
  int<lower=0> N;          //number of watersheds in original dataset
  vector[N] srep;          //srep estimated for original watersheds using alpha prime and beta from Ricker SR models
  vector[N] wa;            //corresponding watershed areas for original watersheds
  int<lower=0> K;          //number of new watersheds to predict srep for
  vector[K] nu_wa;       //vector of watershed areas from new watersheds for which we are trying to estimate beta parameters from a Ricker SR model
}

parameters {
  real slope;       //slope of relationship between srep and watershed area
  real intercept;   //intercept of relationship between srep and watershed area
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
    srep[i] ~ normal(mu[i], sigma);
  }
  
  slope ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma ~ cauchy(0, 5);
}

generated quantities {
  vector[N] nu_Y;
  vector[K] nu_srep;
  
  for(i in 1:N){
    nu_Y[i] = normal_rng(mu[i], sigma);
  }
  for(k in 1:K){
    nu_srep[k] = slope * nu_wa[k] + intercept;
  }
}


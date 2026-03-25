// Model with no prior on smax
data {
  int<lower=0> o_N;          //number ocean-type alphas
  int<lower=0> s_N;          //number of stream-type alphas  
  vector[o_N] a_o;           //ocean-type alphas
  vector[s_N] a_s;           //stream-type alphas
}

parameters {
  real alpha_o;       //mean ocean-type alpha
  real alpha_s;       //mean stream-type alpha
  real sigma_o;
  real sigma_s;
}

model {
  alpha_o ~ normal(a_o, sigma_o);
  alpha_s ~ normal(a_s, sigma_s);
  sigma_o ~ cauchy(0, 1);
  sigma_s ~ cauchy(0, 1);
}

generated quantities {
  real r_sto;          //ratio of stream- to ocean-type alphas
  
  r_sto = (alpha_s / alpha_o);
}


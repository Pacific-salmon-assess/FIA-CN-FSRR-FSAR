// Non state space model comparing FMI and CYER; includes random intercepts for indicator stocks

data {
  int<lower=1> N;                          // number of observations
  int<lower=1> N_indicators;               // number of indicators
  array[N] int<lower=1, upper=N_indicators> indicator;  // indicator grouping
  vector<lower=0, upper=1>[N] can_cyer_obs;       // observed response (proportion scale)
  vector[N] logit_fmi_obs;                        // observed predictor (logit scale)
}

transformed data {
  vector[N] can_cyer_logit_obs = logit(can_cyer_obs);  // transform to logit scale
}

parameters {
  real intercept;                          // fixed intercept
  real slope;                              // fixed slope
  vector[N_indicators] z_indicator;        // standardized random effects
  real<lower=0> sigma_indicator;           // SD of random intercepts
  real<lower=0> phi;                       // Beta precision parameter (SAME AS BRMS)
}

transformed parameters {
  vector[N_indicators] indicator_effects = sigma_indicator * z_indicator;
  vector[N] mu_logit;                      // linear predictor on LOGIT scale
  vector[N] mu_prob;                       // mean on PROBABILITY scale
  
  // Construct linear predictor with random effects
  for (n in 1:N) {
    mu_logit[n] = intercept + slope * logit_fmi_obs[n] + indicator_effects[indicator[n]];
    mu_prob[n] = inv_logit(mu_logit[n]);
  }
}

model {
  // Priors
  intercept ~ normal(-1.68, 0.25);
  slope ~ normal(1, 0.25);
  z_indicator ~ std_normal();
  sigma_indicator ~ exponential(1);
  phi ~ exponential(1);
  
  // Likelihood: Beta regression 
  can_cyer_obs ~ beta(mu_prob * phi, (1 - mu_prob) .* phi);
}

generated quantities {
  vector[N] can_cyer_pred;                 // posterior predictive on proportion scale
  vector[N] log_lik;                       // log likelihood for LOO
  
  for (n in 1:N) {
    // Generate predictions using Beta distribution
    real alpha = mu_prob[n] * phi;
    real beta_param = (1 - mu_prob[n]) * phi;
    can_cyer_pred[n] = beta_rng(alpha, beta_param);
    
    // Log likelihood for the observed data
    log_lik[n] = beta_lpdf(can_cyer_obs[n] | alpha, beta_param);
  }
}

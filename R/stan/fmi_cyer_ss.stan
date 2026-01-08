// State space model comparing FMI and CYER; assumes same CV for observation error in both processes with CV (derived from PSC 2006); includes random intercepts for indicator stocks

data {
  int<lower=1> N;                          // number of observations
  int<lower=1> N_indicators;               // number of indicators
  array[N] int<lower=1, upper=N_indicators> indicator;  // indicator grouping
  vector<lower=0, upper=1>[N] can_cyer_obs;       // observed response (proportion scale)
  vector[N] logit_fmi_obs;                        // observed predictor (CENTERED, logit scale)
  vector[N] fmi_obs;                              // observed predictor (CENTERED, proportion scale)
  real<lower=0> mean_cv;                          // single CV for both variables (e.g., 0.1 for 10%)
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
  
  // Latent true values (state-space components) - CENTERED
  vector[N] logit_fmi_true_centered;       // true predictor values (CENTERED, logit scale)
  vector[N] can_cyer_logit_true;           // true response values (logit scale)
}

transformed parameters {
  vector[N_indicators] indicator_effects = sigma_indicator * z_indicator;
  vector[N] mu_logit;                      // linear predictor on LOGIT scale
  vector[N] mu_prob;                       // mean on PROBABILITY scale
  
  // Construct linear predictor with random effects
  // Uses CENTERED predictor (latent true values are centered)
  for (n in 1:N) {
    mu_logit[n] = intercept + slope * logit_fmi_true_centered[n] + indicator_effects[indicator[n]];
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
  
  // State equation: true response follows Beta regression model
  // Convert from logit scale to probability scale for Beta likelihood
  for (n in 1:N) {
    real p_true = inv_logit(can_cyer_logit_true[n]);
    real alpha = mu_prob[n] * phi;
    real beta_param = (1 - mu_prob[n]) * phi;
    // This is equivalent to: can_cyer_true ~ beta(alpha, beta_param)
    // but we're working with the logit-transformed version
    target += beta_lpdf(p_true | alpha, beta_param);
    // Jacobian adjustment for logit transformation
    target += log(p_true) + log(1 - p_true);
  }
  
  // Observation equation: observed predictor (CENTERED)
  // Delta method uses the CENTERED proportion values
  for (n in 1:N) {
    real p_fmi = fmi_obs[n];  // This is CENTERED (already provided as input)
    // Delta method: SD_logit ≈ CV / (p * (1-p))
    real sd_logit_fmi = mean_cv / (p_fmi * (1 - p_fmi));
    logit_fmi_obs[n] ~ normal(logit_fmi_true_centered[n], sd_logit_fmi);
  }
  
  // Observation equation: observed response
  for (n in 1:N) {
    real p_cyer = inv_logit(can_cyer_logit_true[n]);
    // Delta method: SD_logit ≈ CV / (p * (1-p))
    real sd_logit_cyer = mean_cv / (p_cyer * (1 - p_cyer));
    can_cyer_logit_obs[n] ~ normal(can_cyer_logit_true[n], sd_logit_cyer);
  }
}

generated quantities {
  vector[N] can_cyer_pred;                 // posterior predictive on proportion scale
  vector[N] log_lik;                       // log likelihood for LOO
  
  for (n in 1:N) {
    // Generate predictions using Beta distribution
    real alpha = mu_prob[n] * phi;
    real beta_param = (1 - mu_prob[n]) * phi;
    can_cyer_pred[n] = beta_rng(alpha, beta_param);
    
    // Log likelihood for the observed data (response only)
    real p_cyer = inv_logit(can_cyer_logit_true[n]);
    real sd_logit_cyer = mean_cv / (p_cyer * (1 - p_cyer));
    log_lik[n] = normal_lpdf(can_cyer_logit_obs[n] | can_cyer_logit_true[n], sd_logit_cyer);
  }
}

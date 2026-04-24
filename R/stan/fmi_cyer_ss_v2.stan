// State space model comparing FMI and CYER; assumes same CV for observation error in both processes
// Uses UNTRANSFORMED (proportion-scale) FMI as predictor with proper constraints

data {
  int<lower=1> N;                          // number of observations
  int<lower=1> N_indicators;               // number of indicators
  array[N] int<lower=1, upper=N_indicators> indicator;  // indicator grouping
  vector<lower=0, upper=1>[N] can_cyer_obs;       // observed response (proportion scale)
  vector[N] fmi_obs;                              // observed predictor (proportion scale)
  real<lower=0> mean_cv;                          // single CV for both variables (e.g., 0.3 for 30%)
}

transformed data {
  vector[N] can_cyer_logit_obs = logit(can_cyer_obs);  // transform to logit scale
}

parameters {
  real intercept;                          // fixed intercept
  real slope;                              // fixed slope
  vector[N_indicators] z_indicator;        // standardized random effects
  real<lower=0> sigma_indicator;           // SD of random intercepts
  real<lower=0> phi;                       // Beta precision parameter
  
  // Latent true values (state-space components)
  vector<lower=0, upper=1>[N] fmi_true;
  vector[N] can_cyer_logit_true;           // true response values (logit scale)
}

transformed parameters {
  vector[N_indicators] indicator_effects = sigma_indicator * z_indicator;
  vector[N] mu_logit;                      // linear predictor on LOGIT scale
  vector[N] mu_prob;                       // mean on PROBABILITY scale
  
  // Construct linear predictor with random effects
  for (n in 1:N) {
    mu_logit[n] = intercept + slope * fmi_true[n] + indicator_effects[indicator[n]];
    mu_prob[n] = inv_logit(mu_logit[n]);
  }
}

model {
  // Priors
  intercept ~ normal(-2.75, 1.5);
  slope ~ normal(5.5, 2);
  z_indicator ~ std_normal();
  sigma_indicator ~ exponential(1);
  phi ~ exponential(1);
  
  // State equation: true response follows Beta regression model
  for (n in 1:N) {
    real p_true = inv_logit(can_cyer_logit_true[n]);
    real alpha = mu_prob[n] * phi;
    real beta_param = (1 - mu_prob[n]) * phi;
    // Beta likelihood on probability scale
    target += beta_lpdf(p_true | alpha, beta_param);
    // Jacobian adjustment for logit transformation
    target += log(p_true) + log(1 - p_true);
  }
  
  // Observation equation: observed predictor (PROPORTION SCALE)
  for (n in 1:N) {
    real sd_fmi = mean_cv * fmi_true[n];
    fmi_obs[n] ~ normal(fmi_true[n], sd_fmi);
  }
  
  // Observation equation: observed response (LOGIT SCALE)
  for (n in 1:N) {
    real p_cyer = inv_logit(can_cyer_logit_true[n]);
    // Delta method: SD_logit ≈ CV / (p * (1-p))
    real sd_logit_cyer = mean_cv / (p_cyer * (1 - p_cyer));
    can_cyer_logit_obs[n] ~ normal(can_cyer_logit_true[n], sd_logit_cyer);
  }
}

generated quantities {
  vector[N] can_cyer_pred;                 // posterior predictive on proportion scale
  vector[N] fmi_pred;           // posterior predictive for predictor (UNCENTERED)
  vector[N] log_lik;                       // log likelihood for LOO
  
  for (n in 1:N) {
    // Generate predictions using Beta distribution
    real alpha = mu_prob[n] * phi;
    real beta_param = (1 - mu_prob[n]) * phi;
    can_cyer_pred[n] = beta_rng(alpha, beta_param);
    
    // Predictor predictions 
    fmi_pred[n] = fmi_true[n];
    
    // Log likelihood for the observed data (response only)
    real p_cyer = inv_logit(can_cyer_logit_true[n]);
    real sd_logit_cyer = mean_cv / (p_cyer * (1 - p_cyer));
    log_lik[n] = normal_lpdf(can_cyer_logit_obs[n] | can_cyer_logit_true[n], sd_logit_cyer);
  }
}

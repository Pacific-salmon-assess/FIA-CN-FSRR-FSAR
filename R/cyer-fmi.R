## FMI and CWT-based CYER Comparison
# Goal: Predict total exploitation rate as a function of FMI (based on GSI
# sampling of marine fisheries + run recronstruction catch estimates)
# 1) Use raw ERA outputs to calculate total Canadian CYER; take average from 
# marked vs unmarked fisheries
# 2) Calculate mean proportion of ER associated with US harvest since 2009
# Jan 7, 2026

# NOTE ORIGINAL VERSION USED PUBLISHED MORTALITY TABLES, replaced with raw ERA 
# outputs to accommodate missing 2024 and CKO data

library(tidyverse)
library(readxl)

ind_pal <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#54278f")
names(ind_pal) <- c("CHI", "HAR", "MSH", "SHU", "NIC", "CKO")


## CLEAN CYER DATA -------------------------------------------------------------

# CTC fishery names key
fishery_key <- read.csv(
  here::here(
    "data", "ctc", "CTCFisheryDefinitions.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  select(country, fishery_type, era_fishery_name) %>% 
  distinct()

dat_unmrk <- read.csv(
  here::here(
    "data", "ctc", "cyer_FRSR_unmarked_dec_2025.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  rename(era_fishery_name = fishery_group) %>% 
  left_join(., fishery_key, by = "era_fishery_name") %>% 
  # remove empty rows
  filter(!is.na(era_fishery_name)) %>%
  mutate(
    mark = "unmark"
  )

dat_mrk <- read.csv(
  here::here(
    "data", "ctc", "cyer_FRSR_marked_dec_2025.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  rename(era_fishery_name = fishery_group) %>% 
  left_join(., fishery_key, by = "era_fishery_name") %>% 
  # remove empty rows
  filter(!is.na(era_fishery_name)) %>%
  mutate(
    mark = "mark"
  )

cyer_dat1 <- rbind(dat_mrk, dat_unmrk) %>% 
  filter(!era_fishery_name %in% c("XCA ESC STRAY", "XUS ESC STRAY")) %>% 
  group_by(cy, stock_code, mark) %>% 
  mutate(
    total_er = sum(fishery_group_er)
  ) %>% 
  group_by(cy, stock_code, mark, country, total_er) %>% 
  summarize(
    country_er = sum(fishery_group_er)
  ) %>% 
  rename(
    year = cy, indicator = stock_code
  ) %>% 
  mutate(
    smu = case_when(
      grepl("NIC", indicator) ~ "spring_1.2",
      grepl("CKO", indicator) ~ "summer_1.3",
      grepl("SHU", indicator) |  grepl("MSH", indicator) ~ "summer_0.3",
      grepl("CHI", indicator) |  grepl("HAR", indicator) ~ "fall_0.3",
      TRUE ~ NA_character_
    )
  )

cyer_dat <- cyer_dat1 %>% 
  filter(
    country == "Canada"
  ) %>% 
  group_by(year, indicator, smu) %>% 
  summarize(
    can_cyer = mean(country_er)
  ) %>% 
  ungroup()
  

## AMONG STOCK CYER CORRELATIONS -----------------------------------------------

## FOLLOWING NEEDS TO BE REVISED WITH UPDATED CYER FORMATTING ##

dome_cyer_dat <- read.csv(
  here::here("data", "DOM_CYER.csv")
) %>% 
  janitor::clean_names() %>%
  filter(mort_type == "TM") %>% 
  mutate(
    year = as.numeric(catch_year),
    smu = "spr_1.3",
    nation = ifelse(grepl("CA", fishery_group), "can_er", "us_er"),
    indicator = stock
  ) %>% 
  group_by(
    year, smu, indicator, nation
  ) %>% 
  summarize(
    exp_rate = sum(cyer) / 100
  ) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "nation", values_from = "exp_rate") %>% 
  mutate(
    total_er = can_er + us_er
  ) %>% 
  select(
    colnames(cyer_dat)
  )


library(reshape2)
library(corrplot)
library(GGally)

cyer_wide1 <- cyer_dat %>% 
  select(year, indicator, total_er) %>% 
  pivot_wider(names_from = "indicator", values_from = "total_er") 
cyer_wide2 <- rbind(cyer_dat, dome_cyer_dat) %>% 
  filter(year %in% dome_cyer_dat$year) %>%
  select(year, indicator, total_er) %>% 
  pivot_wider(names_from = "indicator", values_from = "total_er") 

cyer_wide_can <- rbind(cyer_dat, dome_cyer_dat) %>% 
  filter(year %in% dome_cyer_dat$year) %>%
  select(year, indicator, can_er) %>% 
  pivot_wider(names_from = "indicator", values_from = "can_er") 


cor_foo <- function (x) {
  # cor_mat <- cor(x %>% select(-year),
  #                use = "pairwise.complete.obs")
  # 
  # corrplot.mixed(cor_mat,
  #                lower = "circle",   # lower triangle = circles
  #                upper = "number",   # upper triangle = numeric values
  #                tl.col = "black",
  #                tl.srt = 45,
  #                diag   = "n")
  # df2 <- 
  ggpairs(x %>% select(-year),
          upper = list(continuous = wrap("cor", size = 4)),
          lower = list(continuous = "points"),
          diag  = list(continuous = "barDiag")) +
    theme(strip.text = element_text(size = 10))
}

cor_foo(cyer_wide1)

png(here::here("figs", "cwt_indicator_total_cyer.png"), height = 6.5, 
    width = 6.5, units = "in", res = 250)
cor_foo(cyer_wide2)
dev.off()

png(here::here("figs", "cwt_indicator_can_cyer.png"), height = 6.5, 
    width = 6.5, units = "in", res = 250)
cor_foo(cyer_wide_can)
dev.off()


## CLEAN FMI DATA --------------------------------------------------------------

# # GSI based run reconstruction derived exploitation rates from a) Dobson et 
# # al. 2020 (2018 and earlier) and b) Fraser River FMI tech memo (2019-2023); 
# # assembled by A Mesmer, updated by C Freshwater to include sector specific 
# # ERs in recent years 
# rr_cyer_dat1 <-  read_xlsx(
#   here::here("data", "dobson_er_data.xlsx"),
#   sheet = 1
# ) %>% 
#   janitor::clean_names() %>% 
#   mutate(
#     year = as.character(year),
#     smu = case_when(
#       smu == "Sp. 1.2" ~ "spring_1.2",
#       smu == "Sp. 1.3" ~ "spring_1.3",
#       smu == "Su. 1.3" ~ "summer_1.3",
#       smu == "Su. 0.3" ~ "summer_0.3",
#       smu == "Fa. 0.3" ~ "fall_0.3",
#       TRUE ~ smu
#     ),
#     # convert percentages to proportional ERs
#     across(
#       where(is.numeric), ~ .x / 100
#     ),
#     source = ifelse(
#       grepl("Dobson", source), "dobson", "fmi_memo"
#     ),
#     year = as.numeric(year)
#   ) 
# 
# rr_cyer_dat_trim <- rr_cyer_dat1 %>% 
#   select(year, smu, fmi = total) 


## import updated FMI data provided by BJ
fmi_dat1 <- read.csv(
  here::here("data", "fmi_updated_dec2025.csv")
) %>% 
  janitor::clean_names() %>%
  select(-starts_with("avg")) %>% 
  pivot_longer(
    cols = starts_with("x"), names_prefix = "x", values_to = "ppn", 
    names_to = "year"
  ) %>% 
  mutate(
    year = as.numeric(year)
  )

fmi_dat <- fmi_dat1 %>% 
  # filter(
  #   grepl("sbc", strata) | grepl("nbc", strata)
  # ) %>% 
  group_by(
    smu, year
  ) %>% 
  summarize(
    fmi = sum(ppn)
  ) %>% 
  ungroup()
  

ggplot(fmi_dat) + 
  geom_point(aes(x = year, y= fmi)) + 
  facet_wrap(~smu) +
  ggsidekick::theme_sleek()


dat <- left_join(cyer_dat, fmi_dat, by = c("year", "smu")) %>% 
  filter(
    !is.na(fmi)
  ) 
saveRDS(dat, here::here("data", "cyer_fmi_dat.rds"))


## EXPLORE COMBINED DATA -------------------------------------------------------

dat <- readRDS(here::here("data", "cyer_fmi_dat.rds")) %>%
  ungroup() %>% 
  mutate(
    sd_fmi = fmi * 0.1, #converts CV to SD to pass to BRMS
    logit_fmi = qlogis(fmi),
    logit_fmi_centered = logit_fmi - mean(logit_fmi),
    sd_logit_fmi = sd_fmi / (fmi * (1 - fmi))
  )


fmi_cyer_cor <- ggplot(dat) +
  geom_point(aes(x = fmi, y = can_cyer, fill = indicator), shape = 21) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  scale_fill_manual(values = ind_pal) +
  ggsidekick::theme_sleek() 

fmi_cyer_cor_logit <- ggplot(dat) +
  geom_point(aes(x = qlogis(fmi), y = qlogis(can_cyer), fill = indicator), 
             shape = 21) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  scale_fill_manual(values = ind_pal) +
  ggsidekick::theme_sleek() 


## STAN VERSION ----------------------------------------------------------------

library(cmdstanr)
library(posterior)

# based on PSC 2008 report with Green River Chinook salmon CV of CYERs
# (reported in )Table 5-7; p 61) estimated as 0.34 (preterminal sport), 0.49 
# (terminal sport) and 0.24 (preterminal troll)
mean_cv <- 0.3 

# Prepare data
stan_data <- list(
  N = nrow(dat),
  N_indicators = length(unique(dat$indicator)),
  indicator = as.integer(factor(dat$indicator)),
  can_cyer_obs = dat$can_cyer,
  logit_fmi_obs = dat$logit_fmi_centered,
  fmi_obs = dat$fmi,  # proportion scale
  mean_cv = mean_cv # single CV for both variables
)


# state space version
mod <- cmdstan_model(
  here::here("R", "stan", "fmi_cyer_ss.stan")
)

# informative init values required for low CV, use here as well for consistency
init_fn <- function(chain_id) {
  # Add chain-specific jitter
  list(
    intercept = -1.5 + rnorm(1, 0, 0.3),
    slope = 0.7 + rnorm(1, 0, 0.2),
    sigma_indicator = abs(0.3 + rnorm(1, 0, 0.15)),
    phi = abs(8 + rnorm(1, 0, 2)),
    z_indicator = rnorm(stan_data$N_indicators, 0, 0.5),
    logit_fmi_true_centered = stan_data$logit_fmi_obs + 
      rnorm(stan_data$N, 0, 0.02),
    can_cyer_logit_true = qlogis(stan_data$can_cyer_obs) + 
      rnorm(stan_data$N, 0, 0.02)
  )
}

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.97
)


# standard version
mod2 <- cmdstan_model(
  here::here("R", "stan", "fmi_cyer.stan")
)

fit2 <- mod2$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.97
)


# state space version with reduced CV
stan_data3 <- stan_data
stan_data3$mean_cv <- 0.1

fit3 <- mod$sample(
  data = stan_data3,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.97,
  max_treedepth = 12,
  init = init_fn  # Use custom initialization
)


fit_list <- list(fit, fit3, fit2)
names_list <- c("state space - 0.3 CV", "state space - 0.1 CV", "normal")


## summarize fitted parameters
sum_dat <- purrr::map2(
  fit_list, names_list,
  function(x, y) {
    fixed_params <- c("intercept", "slope", "sigma_indicator", "phi")
    fixed_summary <- x$summary(variables = fixed_params)
    
    indicator_summary <- x$summary(variables = "indicator_effects")
    
    bind_rows(
      fixed_summary,
      indicator_summary
    ) %>%
      select(variable, mean, median, sd, q5 = `q5`, q95 = `q95`, rhat, ess_bulk, 
             ess_tail) %>% 
      mutate(
        model = y
      )
  }
) %>% 
  bind_rows()



# posterior predictions across range of fmi values
fmi_new <- seq(0.0, 0.5, by = 0.01)
logit_fmi_new <- qlogis(fmi_new)
logit_fmi_new_centered <- logit_fmi_new - mean(dat$logit_fmi)

N_new <- length(fmi_new)


pred_list <- purrr::map2(
  fit_list, names_list,
  function (x, y) {
    draws <- x$draws(format = "df")
    intercept_draws <- draws$intercept
    slope_draws <- draws$slope
    phi_draws <- draws$phi
    n_draws <- length(intercept_draws)
    
    preds <- matrix(NA, nrow = n_draws, ncol = N_new)
    
    for (i in 1:n_draws) {
      for (j in 1:N_new) {
        # Linear predictor on logit scale (average indicator, so no random effect)
        mu_ij <- intercept_draws[i] + slope_draws[i] * logit_fmi_new_centered[j]
        
        # Back-transform to proportion scale
        preds[i, j] <- plogis(mu_ij)
      }
    }
    
    # Summarize preds
    pred_dat <- tibble(
      fmi = fmi_new,
      logit_fmi = logit_fmi_new,
      logit_fmi_centered = logit_fmi_new_centered,
      mean = colMeans(preds),
      median = apply(preds, 2, median),
      sd = apply(preds, 2, sd),
      q5 = apply(preds, 2, quantile, probs = 0.05),
      q25 = apply(preds, 2, quantile, probs = 0.25),
      q75 = apply(preds, 2, quantile, probs = 0.75),
      q95 = apply(preds, 2, quantile, probs = 0.95),
      model = y
    )
    
    out_list <- list(preds, pred_dat)
    names(out_list) <- c("preds_mat", "preds_dat")
    return(out_list)
  }
)
  
pp <- purrr::map(pred_list, ~ .x$preds_dat)
pred_summary <- bind_rows(pp)


# generate observations that show assumed error distribution used in different 
# models
dat_with_cv <- purrr::map2(
  names_list, 
  c(0.3, 0.1, 0),
  function(x, y) {
    dat %>% 
      mutate(
        model = x,
        cv = y,
        fmi_lower = fmi * (1 - cv),
        fmi_upper = fmi * (1 + cv),
        cyer_lower = can_cyer * (1 - cv),
        cyer_upper = can_cyer * (1 + cv)
      )
  }
) %>% 
  bind_rows()

mod_comp_ribbon <- ggplot(pred_summary, aes(x = fmi)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, fill = "blue") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = median), color = "blue", linewidth = 1) +
  geom_errorbar(data = dat_with_cv, 
                aes(x = fmi, ymin = cyer_lower, ymax = cyer_upper),
                width = 0, alpha = 0.4, color = "gray40") +
  geom_errorbarh(data = dat_with_cv,
                 aes(y = can_cyer, xmin = fmi_lower, xmax = fmi_upper),
                 height = 0, alpha = 0.4, color = "gray40")  +
  geom_point(data = dat, aes(y = can_cyer), alpha = 0.6) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  labs(
    x = "FMI",
    y = "CYER") +
  lims(
    x = c(0.01, 0.5),
    y = c(0.01, 0.7)
  ) +
  ggsidekick::theme_sleek() +
  facet_wrap(~model)


png(here::here("figs", "model-comp-ribbon.png"), height = 4.5, 
    width = 7.5, units = "in", res = 250)
mod_comp_ribbon
dev.off()



## export posterior predictions
dim(pred_list[[1]]$preds_mat)

post_preds <- purrr::map2(
  pred_list, names_list,
  function(x, y) {
    post_draws <- x$preds_mat
    colnames(post_draws) <- fmi_new
    post_draws %>%
      as.data.frame() %>%
      mutate(iteration = row_number()) %>%
      pivot_longer(
        cols = -iteration,
        names_to = "fmi",
        values_to = "cyer"
      ) %>% 
      mutate(
        fmi = as.numeric(fmi),
        model = y
      )
  }
)

trim_seq <- seq(0, 0.5, by = 0.05)

write.csv(
  bind_rows(post_preds) %>% 
    filter(
      fmi %in% trim_seq
    ),
  row.names = FALSE,
  here::here(
    "data", "fmi_cyer_posterior_predictions.csv"
  )
)


##### (INCOMPLETE) ######
#### UPDATE TO UNCENTER BEFORE USING ######

## as above but with stock-specific predictions 

# Extract indicator-specific random effects
indicator_names <- unique(dat$indicator)
N_indicators <- length(indicator_names)

indicator_effects_draws <- draws %>%
  select(starts_with("indicator_effects[")) %>%
  as.matrix()

# Prepare observed data with indicator mapping

# Generate predictions for each indicator stock
predictions_by_indicator <- list()

for (ind in 1:N_indicators) {
  pred_matrix <- matrix(NA, nrow = n_draws, ncol = N_new)
  
  for (i in 1:n_draws) {
    # Get this indicator's random effect
    ind_effect <- indicator_effects_draws[i, ind]
    
    for (j in 1:N_new) {
      # Linear predictor on logit scale
      mu_ij <- intercept_draws[i] + slope_draws[i] * logit_fmi_new[j] + ind_effect
      
      # Add process error
      logit_pred <- rnorm(1, mu_ij, 1 / sqrt(phi_draws[i]))
      
      # Back-transform to proportion scale
      pred_matrix[i, j] <- plogis(logit_pred)
    }
  }
  
  # Summarize predictions for this indicator
  predictions_by_indicator[[ind]] <- tibble(
    indicator = indicator_names[ind],
    fmi = fmi_new,
    logit_fmi = logit_fmi_new,
    mean = colMeans(pred_matrix),
    median = apply(pred_matrix, 2, median),
    sd = apply(pred_matrix, 2, sd),
    q5 = apply(pred_matrix, 2, quantile, probs = 0.05),
    q25 = apply(pred_matrix, 2, quantile, probs = 0.25),
    q75 = apply(pred_matrix, 2, quantile, probs = 0.75),
    q95 = apply(pred_matrix, 2, quantile, probs = 0.95)
  )
}

# Combine all predictions
all_predictions <- bind_rows(predictions_by_indicator)

ggplot() +
  geom_line(data = all_predictions, 
            aes(x = fmi, y = median, color = indicator),
            linewidth = 1) +
  geom_ribbon(data = all_predictions, 
              aes(x = fmi, ymin = q25, ymax = q75, fill = indicator),
              alpha = 0.15) +
  geom_point(data = fitted_summary, 
             aes(x = fmi_prop, y = can_cyer, color = indicator),
             size = 2.5, alpha = 0.7) +
  labs(
    x = "Foreign Marine Index (FMI)",
    y = "Canadian Exploitation Rate",
    title = "All Indicator Stocks: Predictions and Observations",
    subtitle = "Lines = median predictions, ribbons = 50% credible intervals",
    color = "Indicator",
    fill = "Indicator"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )


## BRMS VERSION ----------------------------------------------------------------

# replaced with Stan version above to account for observation error in both CYER
# and FMI

# fit hiearchical slopes model to fmi data with intercept fixed at 0
library(brms)
library(DHARMa)


# random ints
fit_brms1 <- brm(
  can_cyer ~ me(logit_fmi_centered, sd_logit_fmi) + (1 | indicator),
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(normal(0, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95)
  )
fit_brms1b <- brm(
  can_cyer ~ logit_fmi_centered + (1 | indicator), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(normal(0, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)


# random slopes and ints
fit_brms2 <- brm(
  can_cyer ~ me(logit_fmi_centered, sd_logit_fmi) + (me(logit_fmi_centered, sd_logit_fmi) | indicator),
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 5), class = "b"),  # Priors for fixed effects
            prior(exponential(1), class = "sd"),
            prior(normal(-1.68, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)
fit_brms2b <- brm(
  can_cyer ~ logit_fmi_centered + (logit_fmi_centered | indicator), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(exponential(1), class = "sd"),
            prior(normal(-1.68, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)


# constrained to be nearly through 0 with strong informative prior
fit_brms3 <- brm(
  can_cyer ~ me(logit_fmi_centered, sd_logit_fmi) + (1 | indicator),
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(
    # very informative prior on the fixed‐effect intercept
    # equal to mean(dat$logit_fmi)
    prior(normal(-1.68, 0.25), class = "Intercept"),
    # weakly informative prior on the fixed slope
    prior(normal(1, 0.25), class = "b"),
    # very tight zero‐centered prior on the SD of the random intercept
    prior(exponential(1), class = "sd", group = "indicator",
          coef = "Intercept"),
    # prior on the Beta‐precision
    prior(exponential(1), class = "phi")
  ),
  chains = 4, cores = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95)
  )

fit_brms3b <- brm(
  can_cyer ~ logit_fmi_centered + (1 | indicator), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(
    # very informative prior on the fixed‐effect intercept
    # equal to mean(dat$logit_fmi)
    prior(normal(-1.68, 0.25), class = "Intercept"),
    # weakly informative prior on the fixed slope
    prior(normal(1, 0.25), class = "b"),
    # very tight zero‐centered prior on the SD of the random intercept
    prior(exponential(1), class = "sd", group = "indicator",
          coef = "Intercept"),
    # prior on the Beta‐precision
    prior(exponential(1), class = "phi")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)



pred_dat <- expand.grid(
  indicator = unique(dat$indicator),
  fmi = seq(0.01, 0.5, length.out = 30),
  sd_logit_fmi = median(dat$sd_logit_fmi)
)
pred_dat$logit_fmi <- qlogis(pred_dat$fmi)
pred_dat$logit_fmi_centered <- pred_dat$logit_fmi - mean(dat$logit_fmi)

fit_list <- list(
  fit_brms1, 
  fit_brms1b, 
  fit_brms2, 
  fit_brms2b,
  fit_brms3,
  fit_brms3b
)

mean_dat <- purrr::map2(
  fit_list,
  c("rand_i_me", 
    "rand_i", 
    "rand_s_me",
    "rand_s",
    "rand_i_informative_me", 
    "rand_i_informative"
    ),
  function (x, y) {
    pred1 <- predict(x, newdata = pred_dat)
    pred_dat2 <- cbind(pred_dat, pred1)
    
    global_pred <- pred_dat %>% 
      filter(indicator == "SHU")
    pred_fixed <- predict(x, newdata = global_pred, re.form = NA)
    global_pred2 <- cbind(global_pred, pred_fixed) %>% 
      mutate(indicator = "global")
  
    real_pred <- rbind(pred_dat2, global_pred2) %>% 
      mutate(model = y)
    
    # as above but on link scale
    pred1_link <- fitted(x, newdata = pred_dat, scale = "linear")
    pred_dat2_link <- cbind(pred_dat, pred1_link)
    
    pred_fixed_link <- fitted(x, newdata = global_pred, 
                         re.form = NA, scale = "linear")
    global_pred2_link <- cbind(global_pred, pred_fixed_link) %>% 
      mutate(indicator = "global")
    
    link_pred <- rbind(pred_dat2_link, global_pred2_link) %>% 
      mutate(model = y)
    
    real_pred %>% 
      mutate(
        link_est = link_pred$Estimate
      )
  }
) %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "rand_i", 
        "rand_i_me", 
        "rand_s", 
        "rand_s_me", 
        "rand_i_informative", 
        "rand_i_informative_me"
    ))
  )

pred_cyer_ribbon <- fmi_cyer_cor +
  geom_line(data = mean_dat %>% filter(!indicator == "global"),
            aes(x = fmi, y = Estimate, group = indicator),
            linetype = 2) +
  geom_line(data = mean_dat %>% filter(indicator == "global"),
            aes(x = fmi, y = Estimate)) +
  geom_ribbon(data = mean_dat %>% filter(indicator == "global"),
              aes(x = fmi, ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  facet_wrap(~model, ncol = 2) +
  labs(y = "Predicted CWT-based CYER", x = "FMI-based ER") +
  theme(legend.position = "top")

pred_cyer_ribbon_logit <- fmi_cyer_cor_logit +
  geom_line(data = mean_dat %>% filter(!indicator == "global"),
            aes(x = logit_fmi, y = link_est, group = indicator),
            linetype = 2) +
  geom_line(data = mean_dat %>% filter(indicator == "global"),
            aes(x = logit_fmi, y = link_est)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  facet_wrap(~model, ncol = 2) +
  labs(y = "Predicted CWT-based CYER", x = "FMI-based ER") +
  theme(legend.position = "top")
  

pred_list <- purrr::map2(
  fit_list,
  c("rand_i_me", 
    "rand_i", 
    "rand_s_me",
    "rand_s",
    "rand_i_informative_me", 
    "rand_i_informative"
  ),
  function (x, y) {
    new_data <- data.frame(fmi = c(0.05, 0.3)) %>%  # Example new fmi values
      mutate(
        logit_fmi = qlogis(fmi),
        sd_logit_fmi = median(dat$sd_logit_fmi),
        logit_fmi_centered = logit_fmi - mean(dat$logit_fmi)
      )

    pp <- posterior_predict(
      x, newdata = new_data, allow_new_levels = TRUE
    )
    data.frame(
      exp_rate = rep(new_data$fmi, each = nrow(pp)),
      est_cyer = as.numeric(pp),
      model = y
    ) %>% 
      group_by(exp_rate) %>% 
      mutate(
        median_cyer = median(est_cyer)
      )
  }
)


dd <- bind_rows(pred_list) %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "rand_i", 
        "rand_i_me", 
        "rand_s", 
        "rand_s_me", 
        "rand_i_informative", 
        "rand_i_informative_me"
      ))
  )
pred_cyer_ridges <- ggplot(dd) +
  # ggridges::geom_density_ridges() +
  ggridges::stat_density_ridges(
    aes(x = est_cyer, y = model),
    quantile_lines = TRUE,
    quantiles = 2,  # This gives you the median (50th percentile)
    alpha = 0.7
  ) +
  geom_vline(data = dd, 
             aes(xintercept = exp_rate), colour = "red") +
  facet_wrap(~exp_rate) +
  labs(x = "Predicted CWT-based CYER", y = "Model") +
  ggsidekick::theme_sleek()

pred_cyer_violin <- ggplot() +
  geom_violin(data = dd, aes(x = est_cyer, y = model), draw_quantiles = 0.5) +
  geom_vline(data = dd,
             aes(xintercept = exp_rate), colour = "red") +
  facet_wrap(~exp_rate, ncol = 1) +
  labs(x = "Predicted CWT-based CYER", y = "Model") +
  ggsidekick::theme_sleek()


# sim_dharma <- createDHARMa(
#   simulatedResponse = posterior_predict(fit_brms2, nsim = 1000) %>% t(),  # Simulate responses
#   observedResponse = dat$can_er,  # Use actual response values
#   fittedPredictedResponse = fitted(fit_brms2)[,1]  # Extract fitted values
# )
# plot(sim_dharma)  # Check residuals


calc_pit <- function(y, posterior_pred) {
  # Get the proportion of posterior samples that are less than or equal to the observed value
  n_obs <- length(y)
  pit_residuals <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    # Calculate pmin and pmax for each observation
    y_prime <- posterior_pred[i, ]  # posterior predictions for i-th observation
    pmin_i <- mean(y_prime < y[i])
    pmax_i <- mean(y_prime <= y[i])
    
    # Generate the PIT residuals as a random draw from uniform(pmin, pmax)
    pit_residuals[i] <- runif(1, pmin_i, pmax_i)
  }
  
  return(pit_residuals)
}

# check pit resid
purrr::map2(
  fit_list,
  c("rand_i_me", 
    "rand_i", 
    "rand_s_me",
    "rand_s",
    "rand_i_informative_me", 
    "rand_i_informative"),
  function (x, y) {
    pp <- posterior_predict(x, nsim = 1000) %>% t()
    pit_residuals <- calc_pit(y = dat$can_er, posterior_pred = pp)
    qqplot(qunif(ppoints(length(pit_residuals))), pit_residuals,
           main = paste("QQ-plot of PIT Residuals Model ", y))
    abline(0, 1)
  }
)

png(here::here("figs", "pred_cyer_ribbon.png"), height = 6, 
    width = 7.5, units = "in", res = 250)
pred_cyer_ribbon
dev.off()

png(here::here("figs", "pred_cyer_ribbon_logit.png"), height = 6, 
    width = 7.5, units = "in", res = 250)
pred_cyer_ribbon_logit
dev.off()

png(here::here("figs", "cyer-ts.png"), height = 4.5, 
    width = 7.5, units = "in", res = 250)
pred_cyer_ridges
dev.off()

png(here::here("figs", "cyer-ts-violin.png"), height = 4.5, 
    width = 7.5, units = "in", res = 250)
pred_cyer_violin
dev.off()


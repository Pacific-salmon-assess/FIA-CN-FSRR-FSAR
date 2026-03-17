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

# dome_cyer_dat <- read.csv(
#   here::here("data", "DOM_CYER.csv")
# ) %>% 
#   janitor::clean_names() %>%
#   filter(mort_type == "TM") %>% 
#   mutate(
#     year = as.numeric(catch_year),
#     smu = "spr_1.3",
#     nation = ifelse(grepl("CA", fishery_group), "can_er", "us_er"),
#     indicator = stock
#   ) %>% 
#   group_by(
#     year, smu, indicator, nation
#   ) %>% 
#   summarize(
#     exp_rate = sum(cyer) / 100
#   ) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = "nation", values_from = "exp_rate") %>% 
#   mutate(
#     total_er = can_er + us_er
#   ) %>% 
#   select(
#     colnames(cyer_dat)
#   )
# 
# 
# library(reshape2)
# library(corrplot)
# library(GGally)
# 
# cyer_wide1 <- cyer_dat %>% 
#   select(year, indicator, total_er) %>% 
#   pivot_wider(names_from = "indicator", values_from = "total_er") 
# cyer_wide2 <- rbind(cyer_dat, dome_cyer_dat) %>% 
#   filter(year %in% dome_cyer_dat$year) %>%
#   select(year, indicator, total_er) %>% 
#   pivot_wider(names_from = "indicator", values_from = "total_er") 
# 
# cyer_wide_can <- rbind(cyer_dat, dome_cyer_dat) %>% 
#   filter(year %in% dome_cyer_dat$year) %>%
#   select(year, indicator, can_er) %>% 
#   pivot_wider(names_from = "indicator", values_from = "can_er") 
# 
# 
# cor_foo <- function (x) {
#   # cor_mat <- cor(x %>% select(-year),
#   #                use = "pairwise.complete.obs")
#   # 
#   # corrplot.mixed(cor_mat,
#   #                lower = "circle",   # lower triangle = circles
#   #                upper = "number",   # upper triangle = numeric values
#   #                tl.col = "black",
#   #                tl.srt = 45,
#   #                diag   = "n")
#   # df2 <- 
#   ggpairs(x %>% select(-year),
#           upper = list(continuous = wrap("cor", size = 4)),
#           lower = list(continuous = "points"),
#           diag  = list(continuous = "barDiag")) +
#     theme(strip.text = element_text(size = 10))
# }
# 
# cor_foo(cyer_wide1)
# 
# png(here::here("figs", "cwt_indicator_total_cyer.png"), height = 6.5, 
#     width = 6.5, units = "in", res = 250)
# cor_foo(cyer_wide2)
# dev.off()
# 
# png(here::here("figs", "cwt_indicator_can_cyer.png"), height = 6.5, 
#     width = 6.5, units = "in", res = 250)
# cor_foo(cyer_wide_can)
# dev.off()


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
    fmi_centered = fmi - mean(fmi),
    # sd_fmi = fmi * 0.1, #converts CV to SD to pass to BRMS
    logit_fmi = qlogis(fmi),
    logit_fmi_centered = logit_fmi - mean(logit_fmi),
    # sd_logit_fmi = sd_fmi / (fmi * (1 - fmi))
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
  # logit_fmi_obs = dat$logit_fmi_centered,
  fmi_obs = dat$fmi,  # proportion scale
  mean_cv = mean_cv # single CV for both variables
)


# state space version
# switch to no logit transform on predictor Mar 2026 after reviews; also means
# predictor is not centered
mod <- cmdstan_model(
  here::here("R", "stan", "fmi_cyer_ss_v2.stan")
)

# informative init values required for low CV, use here as well for consistency
init_fn <- function(chain_id) {
  # Add chain-specific jitter
  list(
    intercept = -2.75 + rnorm(1, 0, 0.5),   # Center on prior mean ± some noise
    slope = 5.5 + rnorm(1, 0, 0.75),        # Center on prior mean ± some noise
    sigma_indicator = abs(0.3 + rnorm(1, 0, 0.15)),
    phi = abs(8 + rnorm(1, 0, 2)),
    z_indicator = rnorm(stan_data$N_indicators, 0, 0.5),
    fmi_true = stan_data$fmi_obs + 
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

fit$summary(variables = c("intercept", "slope", "sigma_indicator", "phi"))


# standard version
mod2 <- cmdstan_model(
  here::here("R", "stan", "fmi_cyer_v2.stan")
)

fit2 <- mod2$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.97
)
fit2$summary(variables = c("intercept", "slope", "sigma_indicator", "phi"))


# state space version with reduced CV (even more sensitive)
stan_data3 <- stan_data
stan_data3$mean_cv <- 0.1

init_fn_tight <- function(chain_id) {
  # Initialize latent FMI VERY close to observations (CV=0.1 is tight!)
  fmi_init <- stan_data3$fmi_obs + rnorm(stan_data3$N, 0, 0.003)  # Much smaller!
  
  # Ensure within bounds
  fmi_init <- pmax(0.001, pmin(0.999, fmi_init))
  
  list(
    intercept = -2.75,  
    slope = 5.5,        
    sigma_indicator = 0.2,
    phi = 10,
    z_indicator = rep(0, stan_data3$N_indicators),
    fmi_true = fmi_init,  # Very close to observations
    can_cyer_logit_true = qlogis(stan_data3$can_cyer_obs) + 
      rnorm(stan_data3$N, 0, 0.005)  # Much smaller!
  )
}

fit3 <- mod$sample(
  data = stan_data3,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,      # Longer warmup
  iter_sampling = 2000,
  adapt_delta = 0.98,     # Much higher (was 0.97)
  max_treedepth = 13,      # Higher (was 12)
  init = init_fn_tight,
  refresh = 100
)
fit3$summary(variables = c("intercept", "slope", "sigma_indicator", "phi"))


fit_list <- list(fit, fit3, fit2)
names_list <- c("state space - 0.3 CV", "state space - 0.1 CV", "non-state space")

saveRDS(fit_list, here::here("data", "fmi_cyer_fits_v2.rds"))

fit_list <- readRDS(here::here("data", "fmi_cyer_fits_v2.rds"))



## summarize fitted parameters
fixed_params <- c("intercept", "slope", "sigma_indicator", "phi")
sum_dat <- purrr::map2(
  fit_list, names_list,
  function(x, y) {
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
# logit_fmi_new <- qlogis(fmi_new)
# logit_fmi_new_centered <- logit_fmi_new - mean(dat$logit_fmi)

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
        # Linear predictor on logit scale
        mu_ij <- intercept_draws[i] + slope_draws[i] * fmi_new[j] 
        
        # Convert to probability scale
        mu_prob <- plogis(mu_ij)
        
        # Beta distribution parameters
        alpha <- mu_prob * phi_draws[i]
        beta_param <- (1 - mu_prob) * phi_draws[i]
        
        # Generate prediction from Beta distribution (CORRECT)
        preds[i, j] <- rbeta(1, alpha, beta_param)
      }
    }
    
    # Summarize preds
    pred_dat <- tibble(
      fmi = fmi_new,
      # logit_fmi = logit_fmi_new,
      # logit_fmi_centered = logit_fmi_new_centered,
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

ind_pal <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#54278f")
names(ind_pal) <- c("CHI", "HAR", "MSH", "SHU", "NIC", "CKO")


mod_comp_ribbon <- ggplot(pred_summary, aes(x = fmi)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, fill = "gray40") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "gray40") +
  geom_line(aes(y = median), linewidth = 1) +
  geom_vline(xintercept = c(0, 0.1, 0.2, 0.3, 0.4), linetype = 2, alpha = 0.8) +
  geom_errorbar(data = dat_with_cv, 
                aes(x = fmi, ymin = cyer_lower, ymax = cyer_upper),
                width = 0, alpha = 0.6) +
  geom_errorbarh(data = dat_with_cv,
                 aes(y = can_cyer, xmin = fmi_lower, xmax = fmi_upper),
                 height = 0, alpha = 0.6)  +
  geom_point(data = dat, aes(y = can_cyer, fill = indicator), shape = 21, 
             size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  labs(
    x = "FMI",
    y = "CYER",
    fill = "Indicator Stock") +
  lims(
    x = c(0.0, 0.5),
    y = c(0.0, 0.9)
  ) +
  scale_fill_manual(values = ind_pal) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  facet_wrap(~model)


png(here::here("figs", "model-comp-ribbon.png"), height = 4.5, 
    width = 7.5, units = "in", res = 250)
mod_comp_ribbon
dev.off()



## export posterior predictions
pp <- pred_list[[1]]$preds_mat

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

saveRDS(
  bind_rows(post_preds) ,
  here::here(
    "data", "fmi_cyer_posterior_predictions.rds"
  )
)




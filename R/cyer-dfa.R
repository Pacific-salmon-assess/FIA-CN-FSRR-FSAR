## CYER DFA
# Develop multivariate model to generate predictions of CYER for defunct stocks
# based on collinearity
# July 9, 2025

library(tidyverse)
library(MARSS)
library(readxl)
library(reshape2)
library(corrplot)
library(GGally)

## import and clean
sheet_names <- excel_sheets(
  here::here(
    "data", "ctc", 
    "2024 Fraser CWT Detailed Mortality Distribution Tables 16May2025_marked_only.xlsx"
  ))
sheet_ids <- seq_along(sheet_names)

dat <- purrr::map(
  sheet_ids, 
  function(x) {
    read_xlsx(
      here::here(
        "data", "ctc", 
        "2024 Fraser CWT Detailed Mortality Distribution Tables 16May2025_marked_only.xlsx"
      ),
      sheet = x,
      skip = 9#,
      # col_names = FALSE
    ) %>% 
      janitor::clean_names() %>% 
      select(
        stock, catch_year, can_marine_cyer = total_mortality
      ) %>% 
      filter(
        !is.na(stock)
      )
  }
) %>% 
  bind_rows() %>% 
  mutate(
    logit_cyer = boot::logit(can_marine_cyer)
  )

ggplot(dat) +
  geom_point(aes(x = catch_year, y = can_marine_cyer, color = stock))

dat_wide <- dat %>% 
  select(-can_marine_cyer) %>% 
  pivot_wider(
    names_from = "stock", values_from = "logit_cyer"
  ) 

cyer_mat <- dat_wide %>% 
  select(-catch_year) %>% 
  as.matrix() %>% 
  # center (do not scale to account for different variances)
  apply(., 2, function (x) x - mean(x, na.rm = T)) %>% 
  t()



dat %>% 
  select(-c(logit_cyer)) %>% 
  pivot_wider(
    names_from = "stock", values_from = "can_marine_cyer"
  ) %>% 
  select(-c(catch_year)) %>% 
  ggpairs(.,
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = "points"),
        diag  = list(continuous = "barDiag")) +
  theme(strip.text = element_text(size = 10))


## FIT ------------------------------------------------------------------------

## Compare state-space models that included one latent trend and either time series
# specific obs variances or a single obs variance

yrs <- dat$catch_year %>% unique()
n_ts <- nrow(cyer_mat)


# unique obs var
model.list1 <- list(
  Z = matrix("z", nrow = n_ts, ncol = 1),  # each time series loads on single latent trend
  R = "diagonal and unequal",
  A = "zero",
  x0 = "zero",  # initial state is zero
  U = "zero",   # no drift in latent trend
  Q = "unconstrained"
)

fit1 <- MARSS(cyer_mat, model = model.list1)

# shared observation variance
# model.list2 <- model.list1
# model.list2$R <- "diagonal and equal"  # constrain obs errors to be equal
# 
# fit2 <- MARSS(cyer_mat, model = model.list2)
# # fit1 has much lower AIC


# DFA version
model.list_dfa1 <- list(
  m = 1,
  R = "diagonal and unequal",   # observation error for each series
  A = "zero",                   # intercepts; set to "zero" if demeaned
  x0 = "zero",                  # initial state
  # U = "zero",                   # no drift
  Q = "diagonal and unequal",   # each trend evolves independently
  B = "identity"                # random walk structure
)
fit_dfa1 <- MARSS(cyer_mat, model = model.list_dfa1, form = "dfa")


model.list_dfa2 <- list(
  #Z = matrix(list("z1", "z2"), nrow = n_ts, ncol = 2),  # flexible loadings
  m = 2,
  R = "diagonal and unequal",   # observation error for each series
  A = "zero",                   # intercepts; set to "zero" if demeaned
  x0 = "zero",                  # initial state
  # U = "zero",                   # no drift
  Q = "diagonal and unequal",   # each trend evolves independently
  B = "identity"                # random walk structure
)
fit_dfa2 <- MARSS(cyer_mat, model = model.list_dfa2, form = "dfa")


## AIC comparison (simplest state space model favored)
fit1$AIC
fit_dfa1$AIC
fit_dfa2$AIC


## PREDICTIONS -----------------------------------------------------------------

## plot predictions 
mean_dat <- dat %>% 
  group_by(stock) %>% 
  summarize(
    mean_logit_cyer = mean(logit_cyer, na.rm = T)
  )


# predictions from one trend state-space model
pred_dat_ss <- predict(fit1, type = "ytT", interval = "confidence")$pred %>% 
  rename(
    stock = .rownames
  ) %>% 
  left_join(., mean_dat, by = "stock") %>% 
  mutate(
    year = t + 1984,
    real_obs = boot::inv.logit(y + mean_logit_cyer),
    real_estimate = boot::inv.logit(estimate + mean_logit_cyer),
    real_up = boot::inv.logit(estimate + mean_logit_cyer + (1.96 * se)),
    real_lo = boot::inv.logit(estimate + mean_logit_cyer - (1.96 * se))
  )
pred_dat_dfa <- predict(fit_dfa2, type = "ytT", interval = "confidence")$pred %>% 
  rename(
    stock = .rownames
  ) %>% 
  left_join(., mean_dat, by = "stock") %>% 
  mutate(
    year = t + 1984,
    real_obs = boot::inv.logit(y + mean_logit_cyer),
    real_estimate = boot::inv.logit(estimate + mean_logit_cyer),
    real_up = boot::inv.logit(estimate + mean_logit_cyer + (1.96 * se)),
    real_lo = boot::inv.logit(estimate + mean_logit_cyer - (1.96 * se))
  )


marss_ribbon <- ggplot() +
  geom_point(data = pred_dat_ss, aes(x = year, y = real_obs)) +
  geom_line(data = pred_dat_ss, aes(x = year, y = real_estimate), 
            color = "red") +
  geom_ribbon(data = pred_dat_ss, aes(x = year, ymin = real_lo, ymax = real_up), 
              fill = "red", alpha = 0.3) +
  geom_line(data = pred_dat_dfa, aes(x = year, y = real_estimate),
            color = "blue") +
  geom_ribbon(data = pred_dat_dfa, aes(x = year, ymin = real_lo, ymax = real_up), 
              fill = "blue", alpha = 0.3) +
  facet_wrap(
    ~ stock
  ) +
  ggsidekick::theme_sleek()

png(here::here("figs", "cyer_marss.png"), height = 4, width = 6, units = "in",
    res = 250)
marss_ribbon
dev.off()


Z.est <- coef(fit_dfa2, type = "matrix")$Z

# Smoothed latent trends
trends <- MARSSkfss(fit_dfa2)$xtT  # rows = trends, cols = time

# Plot loadings for each time series
barplot(Z.est, beside = TRUE, col = c("darkred", "darkblue"),
        legend = paste0("Trend ", 1:2),
        main = "Factor Loadings on Latent Trends")

matplot(t(trends), type = "l", lty = 1, col = c("darkred", "darkblue"),
        ylab = "Latent Trend", xlab = "Time")
legend("topright", legend = paste0("Trend ", 1:2), col = c("darkred", "darkblue"), lty = 1)

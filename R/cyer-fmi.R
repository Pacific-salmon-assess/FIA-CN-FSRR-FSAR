## FMI and CWT-based CYER Comparison
# Goal: Predict total exploitation rate as a function of FMI (based on GSI
# sampling of marine fisheries + run recronstruction catch estimates)
# 1) Calculate Canadian and total CYER exploitation rates; available as marked
# and unmarked, since FMI is not stratified take mean for CYER
# 2) Calculate mean proportion of ER associated with US harvest since 2009
# Sep. 10, 2024


library(tidyverse)
library(readxl)

ind_pal <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#54278f")
names(ind_pal) <- c("CHI", "HAR", "MSH", "SHU", "NIC", "CKO")


## CLEAN CYER DATA -------------------------------------------------------------

# sheet_names <- excel_sheets(
#   here::here(
#     "data", "ctc", 
#     "TCCHINOOK-25-01-Appendix-C-Mortality-Distribution-Tables-Detailed.xlsx"
#   ))
# 
# stocks <- c("BQR", "CHI", "COW", "HAR", "MSH", "NIC", "PHI", "QUI",
#             "EVIN", "RBT", "NWVI", "SWVI", "SHU")
# matching_sheets <- sheet_names[
#   sapply(sheet_names, function(x) any(grepl(paste(stocks, collapse = "|"), x)))
#   ]
# matching_sheets2 <- matching_sheets[
#   sapply(matching_sheets, function(x) any(grepl("TM", x)))
# ]
# sheet_ids <- which(sheet_names %in% matching_sheets2)
# 
# 
# # CWT based CYERs from Laura Tessier
# # identify sheets w/ relevant data and associated stock name
# new_col_names <- c(
#   "year", "cwt_n", "ages", "aabm_seak_t", "aabm_seak_n", "aabm_seak_s", 
#   "aabm_nbc_t", "aabm_nbc_s", "aabm_wcvi_t", "aabm_wcvi_s", "isbm_nbc_t", 
#   "isbm_nbc_n", "isbm_nbc_s", "isbm_sbc_t", "isbm_sbc_n", "isbm_sbc_s", 
#   "isbm_n_falcon_t", "isbm_n_falcon_s", "isbm_s_falcon_t", "isbm_s_falcon_s", 
#   "isbm_wac_n", "isbm_puget_n", "isbm_puget_s", "term_seak_t", "term_seak_n",
#   "term_seak_s", "term_can_n", "term_can_s", "term_sus_t", "term_sus_n", 
#   "term_sus_s", "stray", "esc", "comment"
# ) 
# 
# cwt_dat <- purrr::map2(
#   sheet_ids, matching_sheets2, 
#   function(x, y) {
#     dum <- read_xlsx(
#       here::here(
#         "data", "ctc", 
#         "TCCHINOOK-25-01-Appendix-C-Mortality-Distribution-Tables-Detailed.xlsx"
#       ),
#       sheet = x,
#       skip = 6,
#       col_names = FALSE
#     )
#     colnames(dum) <- new_col_names
#     dum %>% 
#       mutate(
#         indicator = str_split(y, " ") %>% unlist() %>% .[1],
#         mark = str_split(y, " ") %>% unlist() %>% .[2]
#       ) %>% 
#       # remove five year averages at bottom of table
#       filter(!grepl("-", year))
#   }
# ) %>% 
#   bind_rows()
# 
# # southern US harvest not available for 2023; calculate mean values 2016-22 for
# # each fishery, add to original dataset for 23 then rescale
# cwt_dat_long <- cwt_dat %>% 
#   filter(comment == "ok",
#          # focus only on unmarked given stocks of interest
#          mark == "unmarked") %>% 
#   pivot_longer(cols = c(starts_with("aabm"), starts_with("isbm"),
#                         starts_with("term"), stray, esc),
#                names_to = "strata", values_to = "percent_run") %>% 
#   mutate(
#     year = as.numeric(year),
#     southern_us = ifelse(
#       (grepl("falcon", strata) | grepl("_sus_", strata) | 
#          grepl("puget", strata) | grepl("wac", strata)),
#       TRUE,
#       FALSE
#     ),
#     canadian_er = ifelse(
#       (grepl("nbc", strata) | grepl("sbc", strata) | 
#          grepl("wcvi", strata) | grepl("term_can", strata)),
#       TRUE,
#       FALSE
#     ),
#     missing_from_fmi = ifelse(
#       strata %in% c("isbm_nbc_t", "isbm_nbc_n", "isbm_sbc_t", "isbm_sbc_n",
#                     "term_can_n"),
#       TRUE,
#       FALSE
#     )
#   ) 
# 
# # calculate mean southern US exploitation rate to use since 2023 values 
# # unavailable
# mean_sus <- cwt_dat_long %>% 
#   filter(year > 2015 & year < 2023) %>% 
#   group_by(strata, indicator) %>% 
#   summarize(mean_percent_run = mean(percent_run))
# 
# 
# cwt_dat_long2 <- left_join(cwt_dat_long, mean_sus, 
#                            by = c("indicator", "strata")) %>% 
#   mutate(
#     percent_run = ifelse(year == "2023" & southern_us == TRUE,
#                          mean_percent_run,
#                          percent_run),
#     smu = case_when(
#       indicator == "NIC" ~ "spring_1.2",
#       indicator %in% c("SHU", "MSH") ~ "summer_0.3",
#       indicator %in% c("CHI", "HAR") ~ "fall_0.3",
#       TRUE ~ NA_character_
#     )
#   ) %>% 
#   filter(
#     !is.na(smu)
#   ) %>% 
#   group_by(
#     indicator, year
#   ) %>% 
#   mutate(
#     total_percent = sum(percent_run)
#   ) %>% 
#   ungroup() %>% 
#   mutate(
#     scaled_percent = percent_run / total_percent
#   )


sheet_names <- excel_sheets(
  here::here(
    "data", "ctc",
    "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-unmarked_noTBRSEAK.xlsx"
  ))
# restrict to FR stocks
stocks <- c("CHI", "HAR", "MSH", "NIC", "SHU")

matching_sheets <- sheet_names[
  sapply(sheet_names, function(x) any(grepl(paste(stocks, collapse = "|"), x)))
]
matching_sheets2 <- matching_sheets[
  sapply(matching_sheets, function(x) any(grepl("total mort", x)))
]
sheet_ids <- which(sheet_names %in% matching_sheets2)

# identify sheets w/ relevant data and associated stock name
new_col_names <- c(
  "year", "cwt_n", "ages", "aabm_seak_t", "aabm_seak_n", "aabm_seak_s",
  "aabm_nbc_t", "aabm_nbc_s", "aabm_wcvi_t", "aabm_wcvi_s", "isbm_nbc_t",
  "isbm_nbc_n", "isbm_nbc_s", "isbm_sbc_t", "isbm_sbc_n", "isbm_sbc_s",
  "isbm_n_falcon_t", "isbm_n_falcon_s", "isbm_s_falcon_t", "isbm_s_falcon_s",
  "isbm_wac_n", "isbm_puget_n", "isbm_puget_s", "term_seak_t", "term_seak_n",
  "term_seak_s", "term_can_n", "term_can_s", "term_sus_t", "term_sus_n",
  "term_sus_s", "stray", "esc", "comment"
)
cwt_dat_unmarked <- purrr::map2(
  sheet_ids, matching_sheets2,
  function(x, y) {
    dum <- read_xlsx(
      here::here(
        "data", "ctc",
        "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-unmarked_noTBRSEAK.xlsx"
      ),
      sheet = x,
      skip = 6,
      col_names = FALSE
    )
    colnames(dum) <- new_col_names
    dum %>%
      mutate(
        indicator = str_split(y, " ") %>% unlist() %>% .[1],
        mark = "unmarked"
      ) %>%
      # remove five year averages at bottom of table
      filter(!grepl("-", year),
             ! year == "2024")
  }
) %>%
  bind_rows()

cwt_dat_marked <- purrr::map2(
  sheet_ids, matching_sheets2,
  function(x, y) {
    dum <- read_xlsx(
      here::here(
        "data", "ctc",
        "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-marked_noTBRSEAK.xlsx"
      ),
      sheet = x,
      skip = 6,
      col_names = FALSE
    )
    colnames(dum) <- new_col_names
    dum %>%
      mutate(
        indicator = str_split(y, " ") %>% unlist() %>% .[1],
        mark = "marked"
      ) %>%
      # remove five year averages at bottom of table
      filter(!grepl("-", year),
             ! year == "2024")
  }
) %>%
  bind_rows()

cwt_dat_long <- rbind(cwt_dat_unmarked, cwt_dat_marked) %>%
  filter(comment == "ok") %>%
  mutate(stock = indicator,
         indicator = paste(indicator, mark, sep = "_")) %>%
  pivot_longer(cols = c(starts_with("aabm"), starts_with("isbm"),
                        starts_with("term"), stray, esc),
               names_to = "strata", values_to = "percent_run") %>%
  mutate(
    year = as.numeric(year),
    canadian_er = ifelse(
      (grepl("nbc", strata) | grepl("sbc", strata) |
         grepl("wcvi", strata) | grepl("term_can", strata)),
      TRUE,
      FALSE
    )
  )

cwt_dat_long2 <- cwt_dat_long %>%
  group_by(
    stock, indicator, year
  ) %>%
  mutate(
    total_percent = sum(percent_run)
  ) %>%
  ungroup() %>%
  mutate(
    scaled_percent = percent_run / total_percent
  )


## import Chilko data generated by CW and NS
cko_marked <- read.csv(
  here::here(
    "data", "ctc", "cmz_CKO_marked_dec_2025.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  mutate(
    mark = "marked"
  )
cko_unmarked <- read.csv(
  here::here(
    "data", "ctc", "cmz_CKO_marked_dec_2025.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  mutate(
    mark = "unmarked"
  )
cko_long <- rbind(cko_marked, cko_unmarked) %>% 
  filter(criteria == "ok") %>%
  mutate(indicator = paste(stock, mark, sep = "_")) %>%
  pivot_longer(cols = c(starts_with("aabm"), starts_with("isbm"),
                        starts_with("term"), stray, escap),
               names_to = "strata", values_to = "percent_run") %>%
  mutate(
    year = as.numeric(catch_year),
    canadian_er = ifelse(
      grepl("canada", strata) | grepl("term", strata),
      TRUE,
      FALSE
    )
  ) 

cko_long2 <- cko_long %>%
  group_by(
    stock, indicator, year
  ) %>%
  mutate(
    total_percent = sum(percent_run)
  ) %>%
  ungroup() %>%
  mutate(
    scaled_percent = percent_run / total_percent
  )



# calculate Canadian exploitation rate
can_cyer <- rbind(
  cwt_dat_long2 %>% 
    select(indicator, year, canadian_er, scaled_percent),
  cko_long2 %>% 
    select(indicator, year, canadian_er, scaled_percent)
  ) %>% 
  filter(canadian_er == TRUE) %>% 
  group_by(
    indicator, year 
  ) %>% 
  summarize(
    can_er = sum(scaled_percent)
  )


#combine stray and escapement, then use to calculate total exploitation
cyer_dat <- rbind(
  cwt_dat_long2 %>% 
    select(stock, indicator, year, strata, canadian_er, scaled_percent),
  cko_long2 %>% 
    select(stock, indicator, year, strata, canadian_er, scaled_percent)
) %>% 
  filter(strata %in% c("esc", "stray", "escap")) %>% 
  group_by(stock, year, indicator) %>% 
  summarize(
    percent_escaped = sum(scaled_percent)
  ) %>% 
  ungroup() %>% 
  mutate(
    total_er = 1 - percent_escaped
  ) %>% 
  select(
    stock, indicator, year, total_er
  ) %>% 
  left_join(., can_cyer, by = c("year", "indicator")) %>% 
  mutate(
    us_er = total_er - can_er,
    smu = case_when(
      grepl("NIC", indicator) ~ "spring_1.2",
      grepl("CKO", indicator) ~ "summer_1.3",
      grepl("SHU", indicator) |  grepl("MSH", indicator) ~ "summer_0.3",
      grepl("CHI", indicator) |  grepl("HAR", indicator) ~ "fall_0.3",
      TRUE ~ NA_character_
    )
  ) %>% 
  group_by(
    stock, year, smu
  ) %>% 
  summarize(
    total_er = mean(total_er),
    can_er = mean(can_er),
    us_er = mean(us_er)
  ) %>% 
  ungroup
  


cyer_ts <- ggplot(cyer_dat) + 
  geom_point(aes(x = year, y= can_er)) + 
  geom_point(aes(x = year, y= total_er), color = "red") + 
  facet_wrap(~stock) +
  ggsidekick::theme_sleek()

us_er_box <- ggplot(cyer_dat) + 
  geom_boxplot(aes(x = indicator, y= us_er)) +
  ggsidekick::theme_sleek()



## AMONG STOCK CYER CORRELATIONS -----------------------------------------------

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
    sd_fmi_flat = 0.01, # median of sd_fmi
    sd_logit_fmi = sd_fmi / (fmi * (1 - fmi))
  )


fmi_cyer_cor <- ggplot(dat) +
  geom_point(aes(x = fmi, y = can_er, fill = stock), shape = 21) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  scale_fill_manual(values = ind_pal) +
  ggsidekick::theme_sleek() 

fmi_cyer_cor_logit <- ggplot(dat) +
  geom_point(aes(x = qlogis(fmi), y = qlogis(can_er), fill = stock), 
             shape = 21) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  scale_fill_manual(values = ind_pal) +
  ggsidekick::theme_sleek() 


# fit hiearchical slopes model to fmi data with intercept fixed at 0
library(brms)
library(DHARMa)


# random ints
fit_brms1 <- brm(
  can_er ~ me(logit_fmi, sd_logit_fmi) + (1 | stock),
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(normal(0, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95)
  )
fit_brms1b <- brm(
  can_er ~ logit_fmi + (1 | stock), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(normal(0, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)


# random slopes and ints
# fit_brms2 <- brm(
#   can_er ~ me(marine_fmi, sd_fmi) + (1 + me(marine_fmi, sd_fmi) | stock), 
#   data = dat,
#   family = Beta(link = "logit"),  # Beta regression
#   prior = c(prior(normal(1, 5), class = "b"),  # Priors for fixed effects
#             prior(exponential(2), class = "sd"),
#             prior(normal(0, 2), class = "Intercept"),  # Prior for intercept
#             prior(exponential(1), class = "phi")),
#   chains = 4, cores = 4, iter = 3000, warmup = 1000,
#   control = list(adapt_delta = 0.95)
# )
fit_brms2b <- brm(
  can_er ~ logit_fmi + (1 + logit_fmi | stock), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(prior(normal(1, 1), class = "b"),  # Priors for fixed effects
            prior(exponential(2), class = "sd"),
            prior(normal(0, 1), class = "Intercept"),  # Prior for intercept
            prior(exponential(1), class = "phi")),
  chains = 4, cores = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)


# constrained to be nearly through 0 with strong informative prior
# fit_brms3 <- brm(
#   can_er ~ me(marine_fmi, sd_fmi_flat) + (me(marine_fmi, sd_fmi_flat) | stock), 
#   data = dat,
#   family = Beta(link = "logit"),  # Beta regression
#   prior = c(
#     # very informative prior on the fixed‐effect intercept
#     prior(normal(-10, 0.25), class = "Intercept"),
#     # weakly informative prior on the fixed slope
#     prior(normal(1, 2.5), class = "b"),
#     # very tight zero‐centered prior on the SD of the random intercept
#     prior(exponential(50), class = "sd", group = "stock",
#           coef = "Intercept"),
#     # fairly tight zero‐centered prior on the SD of the random slope
#     prior(exponential(2), class = "sd", group = "stock", 
#           coef = "memarine_fmisd_fmi_flat"),
#     # prior on the Beta‐precision
#     prior(exponential(1), class = "phi")
#   ),
#   chains = 4, cores = 4, iter = 3000, warmup = 1000,
#   control = list(adapt_delta = 0.95)
#   )

fit_brms3b <- brm(
  can_er ~ logit_fmi + (logit_fmi | stock), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(
    # very informative prior on the fixed‐effect intercept
    prior(normal(0, 0.25), class = "Intercept"),
    # weakly informative prior on the fixed slope
    prior(normal(1, 0.25), class = "b"),
    # very tight zero‐centered prior on the SD of the random intercept
    prior(exponential(5), class = "sd", group = "stock",
          coef = "Intercept"),
    # fairly tight zero‐centered prior on the SD of the random slope
    prior(exponential(5), class = "sd", group = "stock", 
          coef = "logit_fmi"),
    # prior on the Beta‐precision
    prior(exponential(1), class = "phi")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)

fit_brms4b <- brm(
  can_er ~ logit_fmi + (logit_fmi | stock), 
  data = dat,
  family = Beta(link = "logit"),  # Beta regression
  prior = c(
    # very informative prior on the fixed‐effect intercept
    prior(normal(0, 0.15), class = "Intercept"),
    # weakly informative prior on the fixed slope
    prior(normal(1, 0.15), class = "b"),
    # very tight zero‐centered prior on the SD of the random intercept
    prior(exponential(5), class = "sd", group = "stock",
          coef = "Intercept"),
    # fairly tight zero‐centered prior on the SD of the random slope
    prior(exponential(5), class = "sd", group = "stock", 
          coef = "logit_fmi"),
    # prior on the Beta‐precision
    prior(exponential(1), class = "phi")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)



pred_dat <- expand.grid(
  stock = unique(dat$stock),
  marine_fmi = seq(0.01, 0.5, length.out = 30),
  se_fmi = 0.1
)
pred_dat$logit_fmi <- qlogis(pred_dat$marine_fmi)

fit_list <- list(
  # fit_brms1, fit_brms2, fit_brms3, 
  fit_brms1b, fit_brms2b, fit_brms3b, fit_brms4b
  # fit_brms1c, fit_brms2c, fit_brms3c
)

mean_dat <- purrr::map2(
  fit_list,
  c(#"rand_i_0.1", "rand_s_0.1", "rand_s_constrained_0.1",
    #"rand_i_0.2", "rand_s_0.2", "rand_s_constrained_0.2",
    "rand_i", "rand_s", "rand_s_constrained", "rand_s_v_constrained"),
  function (x, y) {
    pred1 <- predict(x, newdata = pred_dat)
    # pred_dat$est <- pred1[,1]
    pred_dat2 <- cbind(pred_dat, pred1)
    
    global_pred <- pred_dat %>% 
      filter(stock == "SHU")
    pred_fixed <- predict(x, newdata = global_pred, re.form = NA)
    global_pred2 <- cbind(global_pred, pred_fixed) %>% 
      mutate(stock = "global")
    
    rbind(pred_dat2, global_pred2) %>% 
      mutate(model = y)
  }
) %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(#"rand_i_0.1", "rand_s_0.1", "rand_s_constrained_0.1",
                 #"rand_i_0.2", "rand_s_0.2", "rand_s_constrained_0.2",
                 "rand_i", "rand_s", "rand_s_constrained", "rand_s_v_constrained"
    ))
  )

pred_cyer_ribbon <- fmi_cyer_cor +
  geom_line(data = mean_dat %>% filter(!stock == "global"),
            aes(x = marine_fmi, y = Estimate, group = stock),
            linetype = 2) +
  geom_line(data = mean_dat %>% filter(stock == "global"),
            aes(x = marine_fmi, y = Estimate)) +
  geom_ribbon(data = mean_dat %>% filter(stock == "global"),
              aes(x = marine_fmi, ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  facet_wrap(~model) +
  labs(y = "Predicted CWT-based CYER", x = "FMI-based ER") +
  theme(legend.position = "top")

pred_cyer_ribbon_logit <- fmi_cyer_cor_logit +
  geom_line(data = mean_dat %>% filter(!stock == "global"),
            aes(x = logit_fmi, y = qlogis(Estimate), group = stock),
            linetype = 2) +
  geom_line(data = mean_dat %>% filter(stock == "global"),
            aes(x = logit_fmi, y = qlogis(Estimate))) +
  # geom_ribbon(data = mean_dat %>% filter(stock == "global"),
  #             aes(x = logit_fmi, ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  facet_wrap(~model) +
  labs(y = "Predicted CWT-based CYER", x = "FMI-based ER") +
  theme(legend.position = "top")
  

pred_list <- purrr::map2(
  fit_list,
  c("rand_i", "rand_s", "rand_s_constrained", "rand_s_v_constrained"),
  function (x, y) {
    new_data <- data.frame(fmi = c(0.05, 0.3))  # Example new fmi values
    new_data$logit_fmi <- qlogis(new_data$fmi)
    
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


dd <- bind_rows(pred_list) 
pred_cyer_ridges <- ggplot(dd) +
  ggridges::geom_density_ridges(aes(x = est_cyer, y = model)) +
  geom_vline(data = dd, 
             aes(xintercept = exp_rate), colour = "red") +
  facet_wrap(~exp_rate, ncol = 1) +
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
  c("rand_i", "rand_s", "rand_s_constrained", "rand_s_v_constrained"),
  function (x, y) {
    pp <- posterior_predict(x, nsim = 1000) %>% t()
    pit_residuals <- calc_pit(y = dat_trim$can_er, posterior_pred = pp)
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
    width = 4.5, units = "in", res = 250)
pred_cyer_ridges
dev.off()

png(here::here("figs", "cyer-ts-violin.png"), height = 4.5, 
    width = 4.5, units = "in", res = 250)
pred_cyer_violin
dev.off()


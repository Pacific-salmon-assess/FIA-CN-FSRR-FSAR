
srep_est <- function(nu_wa, old_wa){
  #read in data to run Stan model
  wa_dat <- old_wa
  
  
  #creating the list of data to feed to Stan
  dat <- list("N" = nrow(wa_dat),                  #number of original watersheds for estimating relationship between watershed area and srep
              "K" = nrow(nu_wa),                 #number of new watersheds attempting to estimate a beta value for
              "srep" = log(wa_dat$srep), #sreps from original watersheds
              "wa" = (log(wa_dat$watershed_area) - mean(log(wa_dat$watershed_area))) / sd(log(wa_dat$watershed_area)), #standardized watershed area from the original watersheds
              "nu_wa" = nu_wa$wa) #standardized watershed areas from the watersheds for which we are trying to estimate beta 
  
  #running the Stan model
  fit_srep <- stan(file = "R/srep_wa.stan", data = dat, chains=6,
              iter=10000, cores=6, thin = 1,
              control=list("max_treedepth"=15,"adapt_delta"=0.8),
              pars=c("slope", "intercept", "sigma", "nu_srep"))
  
  #fit summary
  out_srep <- summary(fit_srep)
  nu_srep <- as.data.frame(unlist(extract(fit_srep)$nu_srep))
  print(nu_srep)
  }


#quick check that the line fits the data
# pred <- function(x){
#   (out_srep$summary[1,1] * x + out_srep$summary[2,1])
# }
# 
# 
# wa_dat %>% ggplot(aes(x = (log(watershed_area) - mean(log(watershed_area))) / sd(log(watershed_area)), y = log(srep)))+
#   geom_point()+
#   #uncertainty_ricker+
#   geom_function(fun = pred)+
#   labs(x = "log(Watershed area)", y = "log(Srep)")+
#   theme_classic()



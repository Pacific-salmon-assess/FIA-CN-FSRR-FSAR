###function to estimate correlation among spawner escapements at the individual site/population level

pop_cor <- function(dat){
  cor(dat %>% dplyr::select(pop = Pop_Name, year = Year, escape = TotalInfilled) %>% 
    pivot_wider(values_from = escape, names_from = pop, values_fill = NA) %>% dplyr::select(-year), method = "pearson", use = "pairwise.complete.obs")
}

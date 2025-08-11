library(deSolve)
library(dplyr)
library(tidyr)

set.seed(20250721)
clean_data <-  readRDS("hmpv_data_combined.rds")
RSV_data <- readRDS("nrevss_raw_scaled_RSV.rds")
 

region <- 9
my_region <- 9
fitted_parameters <- readRDS(paste0("2nd_fit_results_0725/parm_", region, ".rds"))


RSV_incidence <- RSV_data %>%  dplyr::filter(regions == region) %>% pull(scaled_cases)

RSV_average <- RSV_data %>%  
  dplyr::filter(regions == region) %>% 
  group_by(epi_week_cdc) %>% 
  summarise(rsv_ave = mean(scaled_cases)) %>% 
  pull(rsv_ave)


hmpv <- clean_data$ts_data %>% 
   filter(regions == region) %>% 
   pull(raw_cases)

hmpv2 <- clean_data$ts_data %>% 
  filter(regions == region) %>% 
  pull(scaled_cases)

age_dist <-  clean_data$age_dist %>% 
  filter(Regions == region) %>% 
  pull(prop)

my_data = list(hmpv = hmpv, 
               prop = round(age_dist * sum(hmpv)))

source("fit_hosp_data_fn.R")
parm_for_fit$data_rsv = c(rep(RSV_average, length.out = length(B_part2)),RSV_incidence)  
parm_for_fit$reporting_fraction = fitted_parameters$reporting_fraction
 



fitLL <- optim(
  par = c(-1.6, 0.1, log(9), -5),
  fn = fit_hosp_data_fn,
  dat = my_data,
  method = "L-BFGS-B",
  lower = c(-Inf, -Inf, log(8.5), -6),  # Only constrain 3rd param
  upper = c(Inf, Inf, log(10), -4),
)


parms <- c(parm_for_fit,
           Amp=plogis(fitLL$par[1]),
           phi=(2*pi*(exp(fitLL$par[2]))) / (1+exp(fitLL$par[2])),
           baseline.txn.rate = exp(fitLL$par[3]),
           interaction = -10^(fitLL$par[4]))
           

 
saveRDS(parms, paste0("parm_", region, "_viral_interference.rds"))
saveRDS(fitLL, paste0("fitLL_", region, "_viral_interference.rds"))

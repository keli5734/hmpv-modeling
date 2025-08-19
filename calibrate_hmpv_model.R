library(deSolve)
library(dplyr)
library(tidyr)

set.seed(20250721)
clean_data <-  readRDS("hmpv_data_combined.rds")
RSV_data <- readRDS("nrevss_raw_scaled_RSV.rds")
 

region <- 8
my_region <- 8
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
#parm_for_fit$baseline.txn.rate = fitted_parameters$baseline.txn.rate


start <- c(b1 = log(0.2),
           trans = 0.1,
           transmission = log(8.7),
           f = log(5e-5))

lower <- c(b1 = log(0.18),
           trans = 0,
           transmission = log(8.5),
           f = log(1e-7))

upper <- c(b1 = log(0.22),
           trans = 0.3,
           transmission = log(9),
           f = log(1e-3))

fitLL <- nlminb(
  start = start,
  objective = fit_hosp_data_fn,   # should return -logLik
  dat = my_data,
  lower = lower,
  upper = upper,
  control = list(iter.max = 1000, trace = TRUE)
)

 


parms <- c(parm_for_fit,
           Amp=exp(fitLL$par[1]),
           phi=(2*pi*(exp(fitLL$par[2]))) / (1+exp(fitLL$par[2])),
           baseline.txn.rate = exp(fitLL$par[3]),
           interaction = -exp(fitLL$par[4]))
           

 
saveRDS(parms, paste0("parm_", region, "_viral_interference.rds"))
saveRDS(fitLL, paste0("fitLL_", region, "_viral_interference.rds"))

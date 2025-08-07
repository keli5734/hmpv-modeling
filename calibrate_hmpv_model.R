library(deSolve)
library(dplyr)
library(tidyr)

set.seed(20250721)
clean_data <-  readRDS("hmpv_data_combined.rds")

region <- 1
my_region <- 1

data_pet <- readRDS("pet.rds") %>% filter(region == my_region) %>% pull(pet)
data_ppt <- readRDS("ppt.rds") %>% filter(region == my_region) %>% pull(ppt)
data_tmin <- readRDS("tmin.rds") %>% filter(region == my_region) %>% pull(tmin)
data_vap <- readRDS("vap.rds") %>% filter(region == my_region) %>% pull(vap)


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


parm_for_fit$data_pet = data_pet
parm_for_fit$data_ppt = data_ppt
parm_for_fit$data_tmin = data_tmin
parm_for_fit$data_vap = data_vap


fitLL <- optim(
  par = c(-1.6, 0.1, log(9), -1, -1.6, -1.6, -1.6, -1.6),
  fn = fit_hosp_data_fn,
  dat = my_data,
  method = "L-BFGS-B",
  lower = c(-Inf, -Inf, log(8.5), -Inf, -Inf,-Inf,-Inf,-Inf)  # Only constrain 3rd param
)


parms <- c(parm_for_fit,
           Amp=plogis(fitLL$par[1]),
           phi=(2*pi*(exp(fitLL$par[2]))) / (1+exp(fitLL$par[2])),
           baseline.txn.rate = exp(fitLL$par[3]),
           reporting_fraction = plogis(fitLL$par[4]), 
           amp_pet=plogis(fitLL$par[5]) * 0,
           amp_ppt=plogis(fitLL$par[6]) * 0,
           amp_tmin=plogis(fitLL$par[7]),
           amp_vap=plogis(fitLL$par[8]) * 0)
          
           

 
saveRDS(parms, paste0("parm_", region, "_climate.rds"))
saveRDS(fitLL, paste0("fitLL_", region, "_climate.rds"))

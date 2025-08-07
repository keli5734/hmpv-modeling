source("parameter_setting.R")
source("hmpv_transmission_model_with_climate.R")

 

fit_hosp_data_fn <- function(parameters, dat) {


  b1 <- parameters[1] # seasonal amplitude
  trans <- parameters[2] # seasonal peak offset
  transmission <- parameters[3] # transmission parameter
  f <- parameters[4] # reporting fraction

  b4 <- parameters[5] # climate effect for tmin
  
  
  
   
  Amp <- plogis(b1)
  phi <- (2*pi*(exp(trans))) / (1+exp(trans))
  baseline.txn.rate <- exp(transmission)
  reporting_fraction <- plogis(f)

  amp_pet <- 0
  amp_ppt <- 0
  amp_tmin <- plogis(b4)
  amp_vap <- 0


  results <- try(
    ode(y = yinit.vector,
        t = my_times,
        func = hmpv_transmission_model_climate,
        parms = c(parm_for_fit,
                  Amp = Amp,
                  phi = phi,
                  baseline.txn.rate = baseline.txn.rate,
                  amp_pet = amp_pet,
                  amp_ppt = amp_ppt,
                  amp_tmin = amp_tmin,
                  amp_vap = amp_vap),
        atol = 1e-5,
        rtol = 1e-5),
    silent = TRUE)

  if (inherits(results, "try-error") || any(!is.finite(results))) {
    return(1e10)
  }
    
  burnN <- t_burn_in
  results.burned <- results[-c(1:burnN),]

  
  #proportion of first infections that are LRI (by age)
  delta1=c(rep(.4,6), # < 6 mos
           rep(.3,7), # 6-24 mos
           rep(.15,3), # 2-5 years
           rep(.05,2), # 5-18 yr
           rep(.05,2), # 18-64 yr
           rep(.2,1)) # > 65 yr
 
  #proportion of second infections that are LRI
  delta2=.5*delta1
  #proportion of third infections that are LRI
  delta3=.5*delta2

  q <- 1
  beta0 <- baseline.txn.rate/(parm_for_fit$dur.days1/7)
  beta <- (beta0/100)/(sum(yinit.matrix)^(1-q))*parm_for_fit$contact
  Amp <-  Amp
  phi <-  phi

  t0 = nrow(results.burned)

  I1 <- results.burned[,grep('I1', colnames(results.burned))]
  I2 <- results.burned[,grep('I2', colnames(results.burned))]
  I3 <- results.burned[,grep('I3', colnames(results.burned))]
  I4 <- results.burned[,grep('I4', colnames(results.burned))]
  S0 <- results.burned[,grep('S0', colnames(results.burned))]
  S1 <- results.burned[,grep('S1', colnames(results.burned))]
  S2 <- results.burned[,grep('S2', colnames(results.burned))]
  S3 <- results.burned[,grep('S3', colnames(results.burned))]


  season <- 1 + Amp * cos(2 * pi * (seq_len(t0) - phi * 52.1775) / 52.1775) +
    amp_pet * data_pet[seq_len(t0)] +
    amp_ppt * data_ppt[seq_len(t0)] +
    amp_tmin * data_tmin[seq_len(t0)] +
    amp_vap * data_vap[seq_len(t0)]

  infectious <- I1 + rho1 * I2 + rho2 * I3 + rho2 * I4
  contact_effect <- infectious %*% beta
  pop_total <- rowSums(results.burned, na.rm = TRUE)
  lambda1 <- sweep(contact_effect, 1, pop_total, "/")
  lambda1 <- sweep(lambda1, 1, season, "*")



  delta1_mat <- matrix(delta1, nrow = t0, ncol = N_ages, byrow = TRUE)
  delta2_mat <- matrix(delta2, nrow = t0, ncol = N_ages, byrow = TRUE)
  delta3_mat <- matrix(delta3, nrow = t0, ncol = N_ages, byrow = TRUE)

  H1 <- delta1_mat * S0 * lambda1 +
    delta2_mat * sigma1 * S1 * lambda1 +
    delta3_mat * sigma2 * S2 * lambda1 +
    delta3_mat * sigma3 * S3 * lambda1

  H_true <- rowSums(H1,na.rm = T)
  
  
  data_fit <- dat$hmpv 
  
  LLcases <- sum(dpois(x = data_fit,
                   lambda = pmax(H_true * reporting_fraction, 1e-12),
                   log = T))
  
  
  agedist = c(sum(colSums(H1, na.rm = T)[1:12]),
              sum(colSums(H1, na.rm = T)[13:16]),
              sum(colSums(H1, na.rm = T)[17:18]),
              sum(colSums(H1, na.rm = T)[19:20]),
              sum(colSums(H1, na.rm = T)[21])) / sum(H1, na.rm = T)
              
   agedist = c(agedist[1], agedist[2], agedist[3],
              agedist[4] + agedist[5]/8,
              7*agedist[5]/8)

  

  epic_age_dist <-  dat$prop

  LLmulti <- dmultinom(x = epic_age_dist,
                       prob = agedist,
                       log = T)
  
  LL <- LLcases + LLmulti  
  
  return(-LL)
}

 

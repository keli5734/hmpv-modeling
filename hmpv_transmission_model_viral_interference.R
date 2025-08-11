###########   define RSV transmission model here #######
hmpv_transmission_model_viral_interference <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    
    
    States <- array(y, dim = dim(parms$yinit.matrix))
    dimnames(States) <- dimnames(parms$yinit.matrix)
    
    # unify the time unit of parameter inputs
    if(parms$time.step =='month'){
      period=12
      length.step=30.44 #days
    }else if(parms$time.step =='week'){
      period=52.1775
      length.step=7 #days
    }
    
    
    omega = 1/(parms$DurationMatImmunityDays/length.step) # wanning immunity from M --> S0
    
    # rate of recovery of the first, second, third and fourth infection
    gamma1 = 1/(parms$dur.days1/length.step) 
    gamma2 = 1/(parms$dur.days2/length.step)
    gamma3 = 1/(parms$dur.days3/length.step)  
    gamma4 = gamma3  
    
    
    # Relative risk of second, third and fourth infections
    sigma1 = parms$sigma1 
    sigma2 = parms$sigma2
    sigma3 = parms$sigma3
    
    
    # Relative infectiousness of second and third infections
    rho1 = parms$rho1 
    rho2 = parms$rho2
    
    
    #Pull out the states for the model as vectors
    M =  States[,'M']  
    S0 =  States[,'S0'] 
    I1 =  States[,'I1']  
    S1 =  States[,'S1']
    I2 =  States[,'I2']
    S2 =  States[,'S2']
    I3 =  States[,'I3']
    S3 =  States[,'S3']
    I4 =  States[,'I4']
    
    N_ages = length(M) # the number of age groups
    
    ## parameter related to force of infection ################
    # per capita transmission probability
    
    #baseline.txn.rate = parms$baseline.txn.rate # <<---------- estimated parameter 1
    
    # transmission probability per unit time
    b <- baseline.txn.rate / (parms$dur.days1/length.step) 
    
    # q depends on transmission type q = 1 here
    q = parms$q  
    
    # c2 is the contact matrix transmission probability per unit time in each age group
    contact = parms$contact  
    
    # transmission rate
    beta = (b/100)/(sum(yinit.matrix)^(1-q))*contact 
    
    # seasonal amplitude and phase shift (estimated)
    #Amp = parms$Amp  
    #phi = parms$phi  
    
    
    RSV_incidence <- parms$data_rsv[time]
    f_int = parms$f_int
    
    #seasonality
    seasonal.txn <- (1+Amp*cos(2*pi*(time-phi*period)/period))
    
    # seasonal transmission probability
    beta_a_i <- seasonal.txn * beta 
    
    infectiousN = (I1+rho1*I2+rho2*I3+rho2*I4) / sum(States)
    
    
    # for frequency dependent transmission
    
    lambda = infectiousN %*% beta_a_i  
    lambda = as.vector(lambda)  
    
    
    
    # create a matrix to record the changing variables
    dy <- matrix(NA, nrow = N_ages, ncol = ncol(States))
    colnames(dy) <- colnames(States)
    
    
    period.birth.rate <- log(parms$PerCapitaBirthsYear[time,]+1)/period
    
    
    # age-specified migration rates
    mu = parms$mu
    
    
    # aging rate (by time step)
    AGE = 1/parms$WidthAgeClassMonth
    if(parms$time.step=='week'){
      AGE = 1/(WidthAgeClassMonth*4.345)} 
    
    delta = c(0,AGE[1:(N_ages-1)])
    
    
    dy[,'M'] <- period.birth.rate * sum(States) - 
      omega * M - 
      AGE * M -
      mu * M[1:(N_ages)] +
      delta*c(0,M[1:(N_ages-1)]) 
    
    dy[,'S0'] <- omega * M -
      (1 + f_int * RSV_incidence)*lambda * S0 -
      AGE * S0 - 
      mu * S0[1:(N_ages)] + 
      delta*c(0,S0[1:(N_ages-1)])  
    
    dy[,'I1'] <-  (1 + f_int * RSV_incidence)*lambda*S0 - 
      gamma1 * I1 - 
      AGE * I1 - 
      mu * I1[1:(N_ages)] + 
      delta*c(0,I1[1:(N_ages-1)]) 
    
    dy[,'S1'] <- gamma1*I1 - 
      (1 + f_int * RSV_incidence)*sigma1*lambda*S1 - 
      AGE * S1 - 
      mu * S1[1:(N_ages)] + 
      delta*c(0,S1[1:(N_ages-1)])  
    
    dy[,'I2'] <- (1 + f_int * RSV_incidence)*sigma1*lambda*S1 - 
      gamma2*I2 -
      AGE * I2 - 
      mu * I2[1:(N_ages)] + 
      delta*c(0,I2[1:(N_ages-1)]) 
    
    dy[,'S2'] <- gamma2*I2 - 
      (1 + f_int * RSV_incidence)*sigma2*lambda*S2 -
      AGE * S2 - 
      mu * S2[1:(N_ages)] + 
      delta*c(0,S2[1:(N_ages-1)]) 
    
    dy[,'I3'] <- (1 + f_int * RSV_incidence)*sigma2*lambda*S2 -
      gamma3 *I3 - 
      AGE * I3 - 
      mu * I3[1:(N_ages)]  +
      delta*c(0,I3[1:(N_ages-1)]) 
    
    dy[,'S3'] <- gamma3*I3 +  
      gamma4*I4 -
      AGE * S3 -
      (1 + f_int * RSV_incidence)*sigma3*lambda*S3 -
      mu * S3[1:(N_ages)] + 
      delta*c(0,S3[1:(N_ages-1)])  
    
    dy[,'I4'] <- (1 + f_int * RSV_incidence)*sigma3*lambda*S3 - 
      gamma4*I4 - 
      AGE * I4 - 
      mu * I4[1:(N_ages)] + 
      delta*c(0,I4[1:(N_ages-1)])  
    
    derivs <- as.vector(dy)
    
    res <- list(derivs)
    
    return(res)
  })
}
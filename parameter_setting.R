 
 

N_ages <- 21



if(region %in% c(8,9)){
  
  
  B_part2 <- rep(0.01365972, each = 52*42)
} else{
  
  B_part2 <- rep(0.01365972, each = 52*41)    
}


 

region_birth_rate <- matrix(NA, nrow = 12, ncol = 10)
region_birth_rate[,1] <- c(0.01168, 0.01137, 
                           0.01099, 0.01065, 0.01056,
                           0.01037, 0.01025, 0.01023,
                           0.01013, 0.01010, 0.00990,
                           0.00970)

region_birth_rate[,2] <- c(0.01329, 0.01300,
                           0.01277, 0.01247, 0.01227,
                           0.01214, 0.01189, 0.01193,
                           0.01184, 0.01174, 0.01147,
                           0.01151)

region_birth_rate[,3] <- c(0.01302, 0.01277,
                           0.01242, 0.01211, 0.01199,
                           0.01189, 0.01171, 0.01180,
                           0.01170, 0.01157, 0.01133,
                           0.01118)

region_birth_rate[,4] <- c(0.01419, 0.01374,
                           0.01311, 0.01252, 0.01227,
                           0.01209, 0.01197, 0.01204,
                           0.01198, 0.01180, 0.01159,
                           0.01134)

region_birth_rate[,5] <- c(0.01349, 0.01320,
                           0.01281, 0.01237, 0.01221,
                           0.01213, 0.01206, 0.01214, 
                           0.01208, 0.01196, 0.01169, 
                           0.01148)


region_birth_rate[,6] <- c(0.01638, 0.01599,
                           0.01556, 0.01479, 0.01428,
                           0.01426, 0.01421, 0.01436, 
                           0.01426, 0.01389, 0.01320,
                           0.01292)

region_birth_rate[,7] <- c(0.01425, 0.01403,
                           0.01371, 0.01327, 0.01303, 
                           0.01304, 0.01290, 0.01297,
                           0.01287, 0.01273, 0.01233, 
                           0.01224)


region_birth_rate[,8] <- c(0.01618, 0.01589, 
                           0.01527, 0.01469, 0.01427, 
                           0.01421, 0.01401, 0.01398, 
                           0.01381, 0.01354, 0.01297, 
                           0.01241)

region_birth_rate[,9] <- c(0.01575, 0.01515, 
                           0.01429, 0.01368, 0.01330,
                           0.01321, 0.01290, 0.01294, 
                           0.01256, 0.01242, 0.01190, 
                           0.01149)

region_birth_rate[,10] <- c(0.01410, 0.01402,
                            0.01347, 0.01299, 0.01275,
                            0.01270, 0.01250, 0.01257, 
                            0.01242, 0.01229, 0.01170,
                            0.01124)


b <-  as.vector(region_birth_rate[,region]) 

years <- seq(from  = 2008, to = 2019, by = 1)
dates <- as.Date(paste0(years, "-07-01"))
date3_interpolate <- seq(as.Date("2008-07-01"), as.Date("2019-07-01"), by = "7 days")

B_part3 <- approx(x = dates,
                  y = b,
                  xout = date3_interpolate)$y

B_part3 <- B_part3 
B <- c(B_part2, B_part3)
B_combined <- matrix(data = 0, nrow = length(B), ncol = N_ages)
B_combined[,1] <- B


contact_mat <- readRDS("contact_mat_POLYMOND.rds")
colnames(contact_mat) <- NULL

WidthAgeClassMonth <-  c(rep(1,times=12), rep(12,times=4),  60, 120, 240, 240, 240)  #Aging rate=1/width age class (months) Vector of long N_age


Population_age_group <- readRDS(paste0("demongraphic_data/Region",region,"/population_age_group_region_",region,".rds"))


### population in 2010 ###
Population_age_group_1 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # under 1
Population_age_group_2 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 1-2
Population_age_group_3 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 2-3
Population_age_group_4 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 3-4
Population_age_group_5 <- subset(Population_age_group$Number, Population_age_group$Age_group == "under 5")[1] / 5 # 4-5
Population_age_group_6 <- subset(Population_age_group$Number, Population_age_group$Age_group == "5-19")[1] / 15 * 5 # 5-9
Population_age_group_7 <- subset(Population_age_group$Number, Population_age_group$Age_group == "5-19")[1] / 15 * 10 # 10-19
Population_age_group_8 <- subset(Population_age_group$Number, Population_age_group$Age_group == "20-39")[1] # 20-39
Population_age_group_9 <- subset(Population_age_group$Number, Population_age_group$Age_group == "40-59")[1] # 40-59
Population_age_group_10 <- subset(Population_age_group$Number, Population_age_group$Age_group == "60+")[1]  


Population_age_under_1 <- round(Population_age_group_1 / 12)


### estimated net migrations rates ###
migration_rates <-  readRDS(paste0("demongraphic_data/Region",region,"/migration_rate_region_",region,".rds"))


migration_rates_gp <- c(rep(migration_rates["mu1"], 16),
                        rep(migration_rates["mu2"], 2), 
                        migration_rates["mu3"], 
                        migration_rates["mu4"], 
                        migration_rates["mu5"])

migration_rates_gp <- as.numeric(migration_rates_gp)


### population in 2011 ###
t_burn_in <-  length(c(B_part2)) #length(B_part1) + 20 * 52   #length(c(B_part1, B_part2)) # round(52.18 * 68) 
t_burn_in2 <- t_burn_in + 52 * 3
birth.rate <- log(mean(c(B_part2))+1)/52.18

#migration_rates <- log(migration_rates+1)/52.18

N0_grp1 <- Population_age_under_1 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # under 1
N0_grp2 <- Population_age_group_2 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 1-2
N0_grp3 <- Population_age_group_3 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 2-3  
N0_grp4 <- Population_age_group_4 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 3-4  
N0_grp5 <- Population_age_group_5 * exp(-(birth.rate+migration_rates["mu1"]) * t_burn_in2)  # 4-5
N0_grp6 <- Population_age_group_6 * exp(-(birth.rate+migration_rates["mu2"]) * t_burn_in2)  # 5-9
N0_grp7 <- Population_age_group_7 * exp(-(birth.rate+migration_rates["mu2"]) * t_burn_in2)  # 10-19
N0_grp8 <- Population_age_group_8 * exp(-(birth.rate+migration_rates["mu3"]) * t_burn_in2)  # 20-39
N0_grp9 <- Population_age_group_9 * exp(-(birth.rate+migration_rates["mu4"]) * t_burn_in2)  # 40-59
N0_grp10 <- Population_age_group_10 * exp(-(birth.rate+migration_rates["mu5"]) * t_burn_in2) # above 60

Pop1 <- c(N1 = round(N0_grp1),
          N2 = round(N0_grp1), 
          N3 = round(N0_grp1),
          N4 = round(N0_grp1),
          N5 = round(N0_grp1),
          N6 = round(N0_grp1),
          N7 = round(N0_grp1),
          N8 = round(N0_grp1), 
          N9 = round(N0_grp1),
          N10 = round(N0_grp1),
          N11 = round(N0_grp1),
          N12 = round(N0_grp1), # < 1
          N13 = round(N0_grp2),
          N14 = round(N0_grp3), 
          N15 = round(N0_grp4),
          N16 = round(N0_grp5),# 1-5
          N17 = round(N0_grp6),
          N18 = round(N0_grp7), # 5-18
          N19 = round(N0_grp8),
          N20 = round(N0_grp9), # 18-60
          N21 = round(N0_grp10)) # >65

N_ages <- length(Pop1) 
agenames <- paste0('Agegrp', 1:N_ages) 
names(Pop1) <- agenames 

## Initialize the compartments (States) 
StateNames <- c("M", 
                "S0", "I1", "S1", "I2", "S2", "I3", "S3", "I4")


States <- array(NA, dim=c(N_ages, length(StateNames) )) #  N age groups x N parameters 
dimnames(States)[[1]] <- agenames
dimnames(States)[[2]] <- StateNames

yinit.matrix <- array(NA, dim=c(N_ages, length(StateNames) ))

dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

yinit.matrix[,c("S1", "I2", "S2", "I3", "S3", "I4")]  <-  0 # setting initial conditions

yinit.matrix[,'M']  <-  c(Pop1[1:6], rep(0,N_ages-6))
yinit.matrix[,'S0'] <-  c(rep(0,6),Pop1[7:N_ages]-rep(1)) 
yinit.matrix[,c('I1')]  <-  c(rep(0,6), rep(1,N_ages-6)) 


yinit.vector <- as.vector(yinit.matrix) #Vectorize the ynit matrix

# Create array that has the labels by age, State and use this to name the yinit.vector
name.array <- array(NA, dim=dim(yinit.matrix)) # dim = 21 x 25 (21 age groups x 25 model compartments) 
for(i in 1:dim(name.array)[1]){ # for 1:21 age groups 
  for(j in 1:dim(name.array)[2]){ # for 1:25 model compartnments (stages)
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i],dimnames(yinit.matrix)[[2]][j])
  }
}

name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector



start_time <- 1 # start date (in week)
tmax <- nrow(B_combined) # end_time (in week)
my_times <- seq(start_time, tmax, by = 1) # gives a sequence from start to end in increments of 1



#########################################
#Relative infectiousness for 2nd and subsequent infections
rho1 <- 0.75
rho2 <- 0.5

#########################################
# duration of infectiousness (months)
#Duration in days
dur.days1 <- 10 #days
dur.days2 <- 7 #days
dur.days3 <- 5 #days


###########################################
# 1/duration of maternal immunity (DAYS)
DurationMatImmunityDays <- 119


#############################################
PerCapitaBirthsYear <- B_combined 


#Relative risk of infection following 1st, 2nd, 3rd+ infections
sigma1 <- 0.7
sigma2 <- 0.55
sigma3 <- 0.4


q <-  1


parm_for_fit <- list(PerCapitaBirthsYear=PerCapitaBirthsYear,
                     WidthAgeClassMonth = WidthAgeClassMonth,
                     DurationMatImmunityDays = DurationMatImmunityDays,
                     mu = migration_rates_gp,
                     rho1=rho1,
                     rho2=rho2,
                     dur.days1=dur.days1,
                     dur.days2=dur.days2,
                     dur.days3=dur.days3,
                     yinit.matrix=yinit.matrix,
                     q=q,
                     contact=contact_mat,
                     sigma1=sigma1,
                     sigma2=sigma2,
                     sigma3= sigma3,
                     time.step = 'week')


 
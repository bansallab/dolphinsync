##########################################################################################
# This code was written by Melissa Collier for the paper "Breathing in sync: how a social behavior structures respiratory epidemic risk in bottlenose dolphins"
# This code takes data from the PC and SB populations and estimates a degree distribution for
# adult male, adult female and juvenile dolphins over a typical DMV infectious period
# This code covers steps 1 and 2 in Fig 1 of the paper

library(glmm)
library(ggplot2)
library(effects)
library(cowplot)
library(MASS)
library(bbmle)
library(effects)
library(sjPlot)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(car)
library(performance)
library(dplyr)
library(interactions)
library(movr)
library(boot)

###### Step 1: Cleaning the Data ####

SB <- read.csv("Clean_Degree_Data_SB.csv", header=T)
SB<- subset(SB, select = c(Dolphin_ID,Follow_length,Degree, Demo, Year, Number_Syncs))
SB$Demo_combined <- ifelse(SB$Demo == "AFNN","AF", SB$Demo) # combining mom calf pairs (AFNN and AFNNN) with Afs
SB$Demo_combined <- ifelse(SB$Demo_combined == "AFNNN","AF", SB$Demo_combined)
SB$Demo_combined <- as.factor(SB$Demo_combined)
SB$Year <- as.factor(SB$Year)
SB$Study_site <- 'SB'
SB$duration_log <- log(SB$Follow_length) #scale follow duration to encourage convergence

PC <- read.csv("Clean_Degree_Data_PCDP.csv", header=T)
PC<- subset(PC, select = c(Dolphin_ID,Follow_length,Degree, Demo, Year, Number_Syncs))
PC$Demo_combined <- ifelse(PC$Demo == "AFNN","AF", PC$Demo) # combining mom calf pairs (AFNN and AFNNN) with Afs
PC$Demo_combined <- ifelse(PC$Demo_combined == "AFNNN","AF", PC$Demo_combined)
PC$Demo_combined <- as.factor(PC$Demo_combined)
PC$Year <- as.factor(PC$Year)
PC$Study_site<- 'PC'
PC$duration_log <- log(PC$Follow_length) #scale follow duration to encourage convergence
PC<- PC[PC$Demo != "AX", ] #drop follows on adult focals where sex could not be even guessed at

#save dataset as the full version
SB_full <- SB 
PC_full<- PC

########### Step 2a: Power Law Generation ##################

#Drop all follows where syncs do not occur
SB_PL<- subset(SB_full, SB$Number_Syncs > 1) 
PC_PL<- subset(PC_full, PC$Number_Syncs > 1)

#drop any outlier data that could skew curve (only one for PC)
PC_PL<- subset(PC, PC$Number_Syncs < 30) # outlier for PC

combined_PL<- rbind(SB_PL, PC_PL)

#Fit and view power law curve
#Since power laws are known to fit data a larger size (set an x min to estimate the distribution)
#We set x min at 10 since this is the average number of syncs predicted by our glmm over our average follow length
plot(combined_PL$Number_Syncs, combined_PL$Degree, xlim = c(0,100), ylim = c(0,40), xlab = "Number of Syncs", ylab = "Degree")
fit.power.law(combined_PL$Number_Syncs,combined_PL$Degree,xmin = 10,plot = TRUE,add = TRUE, col='blue')
PL<-as.data.frame(unlist(fit.power.law(combined_PL$Number_Syncs,combined_PL$Degree,xmin = 10,plot = TRUE,add = TRUE, col='blue')))
names(PL)[1]<- "Estimate"
PL #view alpha and lamnbda parameters

####### Step 2b: Estimating the Gd (prop of time spent in a group for demographic d) ######

## A bootstrap function for proportions
boot_prop <- function(data, idx) {
  
  df <- data[ idx, ]
  #print(df)
  total = length(df$surveys)
  alone_vector<- df$group_size
  a<- length(alone_vector[! alone_vector %in% c('G')])
  
  A <- a/total
  return(A)
}

#### For SB
SB_surveys <- 1:5103 # number of surveys between 2010-2018 where the demo was known for every individual
SB_group<-rep(c("G"),each=2866) # number of surveys on a group (n > 1)
SB_alone<-rep(c("A"),each=2237) # number surveys on a lone individual
SB_group_size <- c(SB_group, SB_alone)

SB_A_df <- as.data.frame(cbind(SB_surveys, SB_group_size))
colnames(SB_A_df) = c("surveys", "group_size")

##Return an overall A with confidence interval for Shark Bay (regardless of demography)
SB_A<- boot(SB_A_df, boot_prop, R=1000)
SB_A_list <- SB_A$t
plot(SB_A) #histogram of bootstrapping results for overall A estimate
boot.ci(SB_A, type = c("norm")) #Get confidence interval

### For each demographic
SB_surveys_AM <- 1:2018 # Total surveys where an adult male is present
SB_group_AM <-rep(c("G"),each=1525) #number of surveys on a group with an adult male
SB_alone_AM <-rep(c("A"),each=493) #number of surveys on a lone adult male
SB_group_size_AM  <- c(SB_group_AM , SB_alone_AM )

SB_AM_A_df <- as.data.frame(cbind(SB_surveys_AM , SB_group_size_AM ))
colnames(SB_AM_A_df) = c("surveys", "group_size")

SB_surveys_AF <- 1:3604 #total surveys where adult female or mom calf pair is present
SB_group_AF <-rep(c("G"),each=2285) #total in a group
SB_alone_AF <-rep(c("A"),each=1319) #total alone (single female or mom calf pair)
SB_group_size_AF  <- c(SB_group_AF , SB_alone_AF )

SB_AF_A_df <- as.data.frame(cbind(SB_surveys_AF , SB_group_size_AF ))
colnames(SB_AF_A_df) = c("surveys", "group_size")

SB_surveys_JX <- 1:1818 #total surveys where a juvenile was present
SB_group_JX <-rep(c("G"),each=1393) #total in a group
SB_alone_JX <-rep(c("A"),each=425) #total alone
SB_group_size_JX  <- c(SB_group_JX , SB_alone_JX )

SB_JX_A_df <- as.data.frame(cbind(SB_surveys_JX , SB_group_size_JX ))
colnames(SB_JX_A_df) = c("surveys", "group_size")

# now combine the dataframes and bootstrap to get a distribution of Ad for each demo group 
Adf_list = list(SB_AM_A_df, SB_AF_A_df, SB_JX_A_df)
A_mega_list = list()

i = 1
for(df in Adf_list){
  
  A<- boot(df, boot_prop, R=1000)
  CI <- boot.ci(A, type = c("norm"))
  print(CI)
  A_mega_list[[i]]<-A$t
  i = i +1
  
}


### Now For PC

PC_surveys <- 1:337 #total surveys in PCDP from 2015-2022
PC_group<-rep(c("G"),each=295) #Total surveys on group 
PC_alone<-rep(c("A"),each=42) #total groups on lone individual or mom calf pair
PC_group_size <- c(PC_group, PC_alone)

PC_A_df <- as.data.frame(cbind(PC_surveys, PC_group_size))
colnames(PC_A_df) = c("surveys", "group_size")

#get distribution for A for PC overall regardless of demography
PC_A<- boot(PC_A_df, boot_prop, R=1000)
boot.ci(PC_A, type = c("norm"))
PC_A_list <- PC_A$t

SB_totals = list(2018, 3604, 1818) #number of surveys where each demo was present, AM, AF, JX

A_PC_demolist <- list()

x = 1

### We use a formula that estimates the G for each demographic group based on the trends of SB and the magnitude of PC
for( A in A_mega_list){
  i = 1
  A_d = list()
  for( a in A){
    
    pc_a = PC_A_list[[i]]
    sb_a = SB_A_list[[i]]
    
    total_pc_alone = round(pc_a * 337)
    surveys_alone_sb_d = round(a * SB_totals[[x]])
    total_sb_alone = round(sb_a * 5103)
    
    surveys_alone_pc_d = round((total_pc_alone * surveys_alone_sb_d) /total_sb_alone)
    
    total_pc = 337
    surveys_sb_d = SB_totals[[x]]
    total_sb = 5103
    
    surveys_pc_d = round((total_pc * surveys_sb_d) / total_sb   ) 
    
    A_PC_d = round(surveys_alone_pc_d / surveys_pc_d , 2)
    
    A_d[[i]] <- A_PC_d
    
    i = i +1
    
  }
  A_PC_demolist[[x]] <- A_d
  
  x = x + 1
}

## Now output the mean P(a) results with sd and min/max values as a csv
PCAM<- c("PCAM", round(mean(unlist(A_PC_demolist[[1]])), 2), round(sd(unlist(A_PC_demolist[[1]])), 2), round(min(unlist(A_PC_demolist[[1]])), 2), round(max(unlist(A_PC_demolist[[1]])), 2))
PCAF<- c("PCAF", round(mean(unlist(A_PC_demolist[[2]])), 2), round(sd(unlist(A_PC_demolist[[2]])), 2), round(min(unlist(A_PC_demolist[[2]])), 2), round(max(unlist(A_PC_demolist[[2]])), 2))
PCJX<- c("PCJX", round(mean(unlist(A_PC_demolist[[3]])), 2), round(sd(unlist(A_PC_demolist[[3]])), 2), round(min(unlist(A_PC_demolist[[3]])), 2), round(max(unlist(A_PC_demolist[[3]])), 2))


SBAM<- c("SBAM", round(mean(unlist(A_mega_list[[1]])), 2), round(sd(unlist(A_mega_list[[1]])), 2), round(min(unlist(A_mega_list[[1]])), 2), round(max(unlist(A_mega_list[[1]])), 2))
SBAF<- c("SBAF", round(mean(unlist(A_mega_list[[2]])), 2), round(sd(unlist(A_mega_list[[2]])), 2), round(min(unlist(A_mega_list[[2]])), 2), round(max(unlist(A_mega_list[[2]])), 2))
SBJX<- c("SBJX", round(mean(unlist(A_mega_list[[3]])), 2), round(sd(unlist(A_mega_list[[3]])), 2), round(min(unlist(A_mega_list[[3]])), 2), round(max(unlist(A_mega_list[[3]])), 2))

df_A_list<- list(PCAM, PCAF, PCJX, SBAM, SBAF, SBJX)
A_df<-as.data.frame(do.call(rbind, df_A_list))
colnames(A_df) <- c("Demo", "A", "SD", "Min", "Max")


### Step 2c: Estimate Sd, prob demo d will sync in a follow  ###########################################
SB_full<- SB_full %>%
  mutate(Occur = case_when(Number_Syncs ==0  ~ "N",
                           Number_Syncs > 0 ~ "Y"
  ))
SB_full$Occur<- as.factor(SB_full$Occur)
SB_full$Demo_combined <- factor(SB_full$Demo_combined, levels = c( "AF","AM", "JX"))

# main model
occurrence_model1SB<- glmer(Occur ~Demo_combined + duration_log+(1|Year)+(1|Dolphin_ID), data = SB_full, family = binomial)
simulationOutput_min <- simulateResiduals(fittedModel = occurrence_model1SB)
testDispersion(simulationOutput_min) #test fit and dispersion
plot(allEffects(occurrence_model1SB)) #check effects

#run additional models with other demos as the intercept to get at differences
SB_full$Demo_combined <- factor(SB_full$Demo_combined, levels = c( "AM","AF", "JX"))
occurrence_model2SB<- glmer(Occur ~Demo_combined + duration_log+(1|Year)+(1|Dolphin_ID), data = SB_full, family = binomial)

#view differences among demographics and study sites
sjPlot:: tab_model(occurrence_model1SB, transform = NULL, digits = 4)
sjPlot:: tab_model(occurrence_model2SB, transform = NULL, digits = 4)

#save the Sd data for SB
eff_occSB <- effects::effect("Demo_combined", occurrence_model2SB)
S_SB <- as.data.frame(eff_occSB)
S_SB

## Now for PC
PC_full<- PC_full %>%
  mutate(Occur = case_when(Number_Syncs ==0  ~ "N",
                           Number_Syncs > 0 ~ "Y"
  ))
PC_full$Occur<- as.factor(PC_full$Occur)
PC_full$Demo_combined <- factor(PC_full$Demo_combined, levels = c( "AF","AM", "JX"))

# main model
occurrence_model1PC<- glmer(Occur ~Demo_combined + duration_log+(1|Year)+(1|Dolphin_ID), data = PC_full, family = binomial)
simulationOutput_min <- simulateResiduals(fittedModel = occurrence_model1PC)
testDispersion(simulationOutput_min) #test fit and dispersion
plot(allEffects(occurrence_model1PC)) #check effects

#run additional models with other demos as the intercept to get at differences
PC_full$Demo_combined <- factor(PC_full$Demo_combined, levels = c( "AM","AF", "JX"))
occurrence_model2PC<- glmer(Occur ~Demo_combined + duration_log+(1|Year)+(1|Dolphin_ID), data = PC_full, family = binomial)

#view differences among demographics and study sites
sjPlot:: tab_model(occurrence_model1PC, transform = NULL, digits = 4)
sjPlot:: tab_model(occurrence_model2PC, transform = NULL, digits = 4)

#save the P(s) data by exporting the model estimates for P(g) by demo/study site
eff_occPC <- effects::effect("Demo_combined", occurrence_model2PC)
S_PC <- as.data.frame(eff_occPC)
S_PC

###### Step 2d: Estimate rd, the rate of syncrony for each demographic group d ####################

### For SB
SB<- SB_full[SB$Number_Syncs != 0, ] #only use follows where a sync occurred

##Run the model
SB$Demo_combined <- factor(SB$Demo_combined, levels = c( "AF", "AM", "JX"))
mod1_rSB<- glmer.nb(Number_Syncs ~Demo_combined+ duration_log+(1|Year)+(1|Dolphin_ID), data = SB)

## Test to make sure not overdispersed
simulationOutput_min <- simulateResiduals(fittedModel = mod1_rSB)
testDispersion(simulationOutput_min)
plot(allEffects(mod1_rSB))
check_convergence(mod1_rSB)

##### run additional models for interpretation with other demos at intercept
SB$Demo_combined <- factor(SB$Demo_combined, levels = c( "AM", "AF", "JX"))
mod2_rSB<- glmer.nb(Number_Syncs ~Demo_combined+ duration_log+(1|Year)+(1|Dolphin_ID), data = SB)

#view significance across demo/study site
sjPlot:: tab_model(mod1_rSB, transform = NULL, digits = 4)
sjPlot:: tab_model(mod2_rSB, transform = NULL, digits = 4)

##Save the data
eff_rSB <- effects::effect("Demo_combined", mod2_rSB)
r_SB <- as.data.frame(eff_rSB)
r_SB

#############  Now for PC

PC<- PC_full[PC$Number_Syncs != 0, ] #drop follows where no syncs occurred

##Run the model
PC$Demo_combined <- factor(PC$Demo_combined, levels = c( "AF", "AM", "JX"))
mod1_rPC<- glmer.nb(Number_Syncs ~Demo_combined+ duration_log+(1|Year)+(1|Dolphin_ID), data = PC)

## Test to make sure not overdispersed
simulationOutput_min <- simulateResiduals(fittedModel = mod1_rPC)
testDispersion(simulationOutput_min)
plot(allEffects(mod1_rPC))
check_convergence(mod1_rPC)

##### run additional models for interpretation with diff demos as intercept
PC$Demo_combined <- factor(PC$Demo_combined, levels = c( "AM", "AF", "JX"))
mod2_rPC<- glmer.nb(Number_Syncs ~Demo_combined+ duration_log+(1|Year)+(1|Dolphin_ID), data = PC)

#view significance across demo/study site
sjPlot:: tab_model(mod1_rPC, transform = NULL, digits = 4)
sjPlot:: tab_model(mod2_rPC, transform = NULL, digits = 4)

##Save the data
eff_rPC <- effects::effect("Demo_combined", mod2_rPC)
r_PC <- as.data.frame(eff_rPC)
r_PC

####### Step 2e: Get degree distributions from powerlaw, Gd, Sd, and rd for PC and SB
## The function below allows you to describe a distribution with a defined min
# max, mean and standard error
rgbeta <- function(n, mean, var, min = 0, max = 1)
{
  dmin <- mean - min
  dmax <- max - mean
  
  if (dmin <= 0 || dmax <= 0)
  {
    stop(paste("mean must be between min =", min, "and max =", max)) 
  }
  
  if (var >= dmin * dmax)
  {
    stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  }
  
  # mean and variance of the standard beta distributed variable
  mx <- (mean - min) / (max - min)
  vx <- var / (max - min)^2
  
  # find the corresponding alpha-beta parameterization
  a <- ((1 - mx) / vx - 1 / mx) * mx^2
  b <- a * (1 / mx - 1)
  
  # generate standard beta observations and transform
  x <- rbeta(n, a, b)
  y <- (max - min) * x + min
  
  return(y)
}
######################## First SB

## We will estimate degree 1000 times for each demo group
n <- 1000
average_follow_length_SB = 119

SB_m_A <- as.numeric(A_df$A[4:6])
SB_sd_A<- as.numeric(A_df$SD[4:6])
SB_min_A <- as.numeric(A_df$Min[4:6])
SB_max_A <- as.numeric(A_df$Max[4:6])

SB_m_S <- as.numeric(S_SB$fit)
SB_sd_S<- as.numeric(S_SB$se)
SB_min_S <- as.numeric(S_SB$lower)
SB_max_S<- as.numeric(S_SB$upper)

SB_m_rate <- as.numeric(r_SB$fit)
SB_sd_rate <- as.numeric(r_SB$se)
SB_min_rate <- as.numeric(r_SB$lower)
SB_max_rate <- as.numeric(r_SB$upper)

alpha = as.numeric(PL[2,])
lambda = as.numeric(PL[3,])

demos<- 1:3
samples<- 1:n
estimated_degrees_SB<- list()

IP_range<- list(5,6,7,8,9,10) #range of potential infectious periods to estimate degree for
day<- 1441 #number of minutes in a day

for(d in demos){
  
  degree_samples<- list()
  
  #### First set A parameters for the approriate demo/population
  m1 = SB_m_A[d]
  sd1 = SB_sd_A[d]
  min1 = SB_min_A[d]
  max1 = SB_max_A[d]
  
  ## make sure the sd makes sense for the function
  test1 <- (m1 - min1) * (max1 - m1)
  if ( sd1 >= test1){
    sd1 = test1-0.0001
  }
  
  #### Next set D parameters for the approriate demo/population
  m2 = SB_m_S[d]
  sd2 = SB_sd_S[d]
  min2 = SB_min_S[d]
  max2 = SB_max_S[d]
  
  ## make sure the sd makes sense for the function
  test2 <- (m2 - min2) * (max2 - m2)
  if ( sd2 >= test2){
    sd2 = test2-0.0001
  }
  
  #### Finally set rate parameters for the appropriate demo/pop
  m3 = SB_m_rate[d]
  sd3 = SB_sd_rate[d]
  min3 = SB_min_rate[d]
  max3 = SB_max_rate[d]
  
  ## make sure the sd makes sense for the function
  
  if (m3 < min3){
    min3 = m3-1
  }
  
  test3 <- (m3 - min3) * (max3 - m3)
  if ( sd3 >= test3){
    sd3 = test3-0.0001
  }
  
  #### Now we draw a P(a), P(s) and rate from distributions based on the above parameters
  for(i in samples){
    
    ## choose infectious period
    IP<- IP_range[[sample(1:length(IP_range), 1)]]
    
    #generate a random A
    A <- rgbeta(1, mean = m1, var = sd1, min = min1, max = max1)
    
    
    #generate a random S
    S <- rgbeta(1, mean = m2, var = sd2, min = min2, max = max2)
    
    #geerate a random sync rate
    rate <-rgbeta(1, mean = m3, var = sd3, min = min3, max = max3)/average_follow_length_SB
    
    # Now calcualte the estimated degree based on G (where G = 1-A), S, and r and the power law curve
    
    degree = round(alpha * ((IP*day*(1-A)*S)*rate)^-lambda)
    degree_samples[[i]]<- degree
    
    
    
  }
  
  estimated_degrees_SB[[d]]<- degree_samples
}

SB_est_deg<- as.data.frame(do.call(cbind, estimated_degrees_SB))
colnames(SB_est_deg) <- c("AM", "AF", "JX")


############### Now for PC
n <- 1000
average_follow_length_PC = 25

PC_m_A <- as.numeric(A_df$A[1:3])
PC_sd_A<- as.numeric(A_df$SD[1:3])
PC_min_A <- as.numeric(A_df$Min[1:3])
PC_max_A <- as.numeric(A_df$Max[1:3])

PC_m_S <- as.numeric(S_PC$fit)
PC_sd_S<- as.numeric(S_PC$se)
PC_min_S <- as.numeric(S_PC$lower)
PC_max_S<- as.numeric(S_PC$upper)

PC_m_rate <- as.numeric(r_PC$fit)
PC_sd_rate <- as.numeric(r_PC$se)
PC_min_rate <- as.numeric(r_PC$lower)
PC_max_rate <- as.numeric(r_PC$upper)

estimated_degrees_PC<- list()

for(d in demos){
  
  degree_samples<- list()
  
  #### First set P(a) parameters for the approriate demo/population
  m1 = PC_m_A[d]
  sd1 = PC_sd_A[d]
  min1 = PC_min_A[d]
  max1 = PC_max_A[d]
  
  ## make sure the sd makes sense for the function
  test1 <- (m1 - min1) * (max1 - m1)
  if ( sd1 >= test1){
    sd1 = test1-0.0001
  }
  
  #### Next set P(s) parameters for the approriate demo/population
  m2 = PC_m_S[d]
  sd2 = PC_sd_S[d]
  min2 = PC_min_S[d]
  max2 = PC_max_S[d]
  
  ## make sure the sd makes sense for the function
  test2 <- (m2 - min2) * (max2 - m2)
  if ( sd2 >= test2){
    sd2 = test2-0.0001
  }
  
  #### Finally set rate parameters for the appropriate demo/pop
  m3 = PC_m_rate[d]
  sd3 = PC_sd_rate[d]
  min3 = PC_min_rate[d]
  max3 = PC_max_rate[d]
  
  ## make sure the sd makes sense for the function
  
  if (m3 < min3){
    min3 = m3-1
  }
  
  test3 <- (m3 - min3) * (max3 - m3)
  if ( sd3 >= test3){
    sd3 = test3-0.0001
  }
  
  #### Now we draw a P(a), P(s) and rate from distributions based on the above parameters
  for(i in samples){
    
    ## choose infectious period
    IP<- IP_range[[sample(1:length(IP_range), 1)]]
    
    #generate a random P(a)
    A <- rgbeta(1, mean = m1, var = sd1, min = min1, max = max1)
    
    
    #generate a random P(s)
    S <- rgbeta(1, mean = m2, var = sd2, min = min2, max = max2)
    
    #geerate a random sync rate
    rate <-rgbeta(1, mean = m3, var = sd3, min = min3, max = max3)/average_follow_length_PC
    
    # Now calcualte the estimated degree based on G (where G= 1-A), S, rate and the power law curve
    
    degree = round(alpha * ((IP*day*(1-A)*S)*rate)^-lambda)
    degree_samples[[i]]<- degree
    
    
    
  }
  
  estimated_degrees_PC[[d]]<- degree_samples
}

PC_est_deg<- as.data.frame(do.call(cbind, estimated_degrees_PC))
colnames(PC_est_deg) <- c("AM", "AF", "JX")


######## Step 2f: Generate degree distributions for each demo for SB and PC over an average infectious period ####################################

combined_est_deg<- cbind(PC_est_deg, SB_est_deg)
cols <- c('PCAM', 'PCAF', 'PCJX', 'SBAM', 'SBAF', 'SBJX')
colnames(combined_est_deg)<- cols
columns<- seq(1,6, by=1)
demos<- colnames(combined_est_deg)
num_samples<- 1000
max_degree<- 300 #max the degree so we don't have issues with network generation

samples<-seq(1, num_samples, by = 1)
all_distrubutions<- list()
freq<- list()
density<- list()

for(column in columns){
  mu<- mean(unlist(combined_est_deg[column]))
  var<- var(unlist(combined_est_deg[column]))
  theta = var
  
  dist<- list()
  
  for( sample in samples ){
    dummy <- rnegbin(1, mu, theta = theta)
    while( dummy > max_degree){
      dummy <- rnegbin(1, mu, theta = theta)
    }
    dist[[sample]]<- dummy
  }
  
  freq[[column]]<- dist
  
  degrees <- seq(0, max_degree, by = 1)
  
  distribution <- list()
  
  for(deg in degrees){
    total <- 0
    for(d in dist) {
      if( d == deg){
        total <- total +1 
      }
    }
    distribution[[deg+1]] <- total/num_samples
  }
  all_distrubutions[[column]]<- distribution
}

##export final dsitrubutions for network generation
distribution_df<- as.data.frame(do.call(cbind, all_distrubutions))
colnames(distribution_df)<- demos
PC_dist<- distribution_df[c(1,2,3)]
PC_dist <- apply(PC_dist,2,as.character)
SB_dist<- distribution_df[c(4,5,6)]
SB_dist <- apply(SB_dist,2,as.character)

write.csv(PC_dist, "PC_degree_distributions.csv")
write.csv(SB_dist, "SB_degree_distributions.csv")

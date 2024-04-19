library(glmm)
library(ggplot2)
library(effects)
library(cowplot)
library(MASS)
library(bbmle)
library(effects)
library(boot)
library(broom)
library(reshape2)
library(psych)
library(expm)
library(stringr)
library(dplyr)

#bootstrapping function for the mean
boot_mean <- function(original_vector, resample_vector) {
  mean(original_vector[resample_vector])
}

### Read in and clean data for Shark Bay #######
SB <- read.csv("Clean_Degree_Data_SB.csv", header=T)
SB<- subset(SB, SB$Demo != 'AX')
SB$Demo_combined <- ifelse(SB$Demo == "AFNN","AF", SB$Demo) #combine AF and mom/calf pair follows
SB$Demo_combined <- ifelse(SB$Demo_combined == "AFNNN","AF", SB$Demo_combined)
SB$Demo_combined <- as.factor(SB$Demo_combined)
SB<- subset(SB, SB$Number_Syncs != 0) #only using follows where syncs occurred

SB_AM_data<- SB[(SB$Demo_combined=="AM"),]
SB_AF_data<- SB[(SB$Demo_combined=="AF"),]
SB_JX_data<- SB[(SB$Demo_combined=="JX"),]

SB_Data_list <- list(SB_AM_data, SB_AF_data,SB_JX_data)
Demo_list <- c('AM', 'AF', 'JX')

## For the mixing matrix
SB_df_empirical <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(SB_df_empirical) <-c("Comb", "Mean","Lower", "Upper")

h <- 1
n <- 1
# this loop takes the degree data from the data import and gets a proportion of degree 
# in each demographic for each focal individual
for (demo in Demo_list) {
  AM <- c(as.data.frame(SB_Data_list[n])[,5])
  AF <-c(as.data.frame(SB_Data_list[n])[,6])
  JX <-c(as.data.frame(SB_Data_list[n])[,7])
  
  mega_list <- list(AM, AF, JX)
  
  p <- 1
  for (i in mega_list){
    
    if (sum(i) == 0) { 
      c <- paste(c(demo,p),collapse="")
      m <- 0 
      l <- 0
      u <- 0
      SB_df_empirical[h,] <- c(c, m, l, u)
      p = p+1
      h = h +1
    } 
    else {
      results <- boot(i, boot_mean, R = 2000)
      m <- round(results$t0, 3)
      ci <- boot.ci(results)
      l <- round(ci$normal[2], 3)
      if (l < 0) {
        l = 0
      }
      u <- round(ci$normal[3], 3)
      c <- paste(c(demo,p),collapse="")
      SB_df_empirical[h,] <- c(c, m, l, u)
      p = p+1
      h = h +1
    }
    
  }
  
  n = n +1
  
}

#Once we have the proportions for all the follows, we must standardize it by
#demographic distribution of the population and the average degree, so that it makes
#sense in the network model
SB_Mean_Ratio <- c()
SB_empirical_Upper_Ratio <- c()
SB_empirical_Lower_Ratio <- c()
SB_pop_deg <- sum(c(as.numeric(as.list(SB_df_empirical$Mean))))
SB_props = c(0.42, 0.44, 0.14) #AM/AF/JX distribution
SB_degs = c(12, 7, 10) #AM/AF/JX average estimated degree
SB_totdeg = round(SB_degs[1]*SB_props[1] + SB_degs[2]*SB_props[2] + SB_degs[3]*SB_props[3]) #population average estimated degree 

#calculate the mixing matrix
r= 1
for (demo in Demo_list) {
  
  mix <- SB_df_empirical %>% filter(str_detect(Comb, demo))
  
  mean_list <- c(as.numeric(as.list(mix$Mean)))
  mean_total<- sum(unlist(mean_list))
  upper_list <- c(as.numeric(as.list(mix$Upper)))
  upper_total <- sum(unlist(upper_list))
  lower_list <- c(as.numeric(as.list(mix$Lower)))
  lower_total <- sum(unlist(lower_list))
  
  for (i in mean_list) {
    x = ((SB_degs[r] * ((i/mean_total)* SB_props[r]))/SB_totdeg)
    SB_Mean_Ratio <- append(SB_Mean_Ratio, x)
  }

  
  for (i in upper_list) {
    x = ((SB_degs[r] * ((i/upper_total)* SB_props[r]))/SB_totdeg)
    if(x > 1){
      x <-1
    }
    SB_empirical_Upper_Ratio <- append(SB_empirical_Upper_Ratio, x)
  }
  
  for (i in lower_list) {
    x = ((SB_degs[r] * ((i/lower_total)* SB_props[r]))/SB_totdeg)
    
    SB_empirical_Lower_Ratio <- append(SB_empirical_Lower_Ratio, x)
  }
  
  r = r+1 
}

SB_SD_Ratio = c()
u = 1
for (i in SB_empirical_Lower_Ratio) {
  x = abs(i -SB_Mean_Ratio[u])
  
  SB_SD_Ratio <- append(SB_SD_Ratio, x)
  u = u+1
}

SB_df_empirical_sep <- cbind(SB_df_empirical, SB_Mean_Ratio, SB_SD_Ratio, SB_empirical_Lower_Ratio, SB_empirical_Upper_Ratio)

SB_temp_data <- c(as.numeric(as.list(SB_df_empirical_sep$SB_Mean_Ratio)))
SB_M_sum <- sum(SB_temp_data)
SB_matrix_data<- SB_temp_data/SB_M_sum
SB_df_empirical_sep$SB_Mean_Ratio<-SB_matrix_data
SB_mixing_matrix <- matrix(SB_matrix_data,nrow=3,ncol=3,byrow=TRUE)

#get assortativity based on Newman's definition
SB_assort <- (tr(SB_mixing_matrix)-sum(SB_mixing_matrix %^% 2))/(1-sum(SB_mixing_matrix %^% 2))
SB_mixing_matrix
SB_assort

write.csv(SB_df_empirical_sep, "Mixing_Matrix_data_SB.csv")


### Now generate MM for PC #############
PC <- read.csv("Clean_Degree_Data_PCDP.csv", header=T)
PC<- subset(PC, PC$Demo != 'AX')
PC$Demo_combined <- ifelse(PC$Demo == "AFNN","AF", PC$Demo) #combine AF and mom/calf pair follows
PC$Demo_combined <- ifelse(PC$Demo_combined == "AFNNN","AF", PC$Demo_combined)
PC$Demo_combined <- as.factor(PC$Demo_combined)
PC<- subset(PC, PC$Number_Syncs != 0)

PC_AM_data<- PC[(PC$Demo_combined=="AM"),]
PC_AF_data<- PC[(PC$Demo_combined=="AF"),]
PC_JX_data<- PC[(PC$Demo_combined=="JX"),]

PC_Data_list <- list(PC_AM_data, PC_AF_data,PC_JX_data)
Demo_list <- c('AM', 'AF', 'JX')

## For the mixing matrix

PC_df_empirical <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(PC_df_empirical) <-c("Comb", "Mean","Lower", "Upper")

h <- 1
n <- 1
# this loop takes the degree data from the data import and gets a proportion of degree 
# in each demographic for each focal individual
for (demo in Demo_list) {
  AM <- c(as.data.frame(PC_Data_list[n])[,6])
  AF <-c(as.data.frame(PC_Data_list[n])[,7])
  JX <-c(as.data.frame(PC_Data_list[n])[,8])
  
  mega_list <- list(AM, AF, JX)
  
  p <- 1
  for (i in mega_list){
    
    if (sum(i) == 0) { 
      c <- paste(c(demo,p),collapse="")
      m <- 0 
      l <- 0
      u <- 0
      PC_df_empirical[h,] <- c(c, m, l, u)
      p = p+1
      h = h +1
    } 
    else {
      results <- boot(i, boot_mean, R = 2000)
      m <- round(results$t0, 3)
      ci <- boot.ci(results)
      l <- round(ci$normal[2], 3)
      if (l < 0) {
        l = 0
      }
      u <- round(ci$normal[3], 3)
      c <- paste(c(demo,p),collapse="")
      PC_df_empirical[h,] <- c(c, m, l, u)
      p = p+1
      h = h +1
    }
    
  }
  
  n = n +1
  
}

#Once we have the proportions for all the follows, we must standardize it by
#demographic distribution of the population and the average degree, so that it makes
#sense in the model
PC_Mean_Ratio <- c()
PC_empirical_Upper_Ratio <- c()
PC_empirical_Lower_Ratio <- c()
PC_pop_deg <- sum(c(as.numeric(as.list(PC_df_empirical$Mean))))
PC_props = c(0.42, 0.44, 0.14) #AM/AF/JX distribution
PC_degs = c(31, 23, 36) #AM/AF/JX average estimated degree
PC_totdeg = round(PC_degs[1]*PC_props[1] + PC_degs[2]*PC_props[2] + PC_degs[3]*PC_props[3]) #population average estimated degree 

#calculate the mixing matrix
r= 1
for (demo in Demo_list) {
  
  mix <- PC_df_empirical %>% filter(str_detect(Comb, demo))
  
  mean_list <- c(as.numeric(as.list(mix$Mean)))
  mean_total<- sum(unlist(mean_list))
  upper_list <- c(as.numeric(as.list(mix$Upper)))
  upper_total <- sum(unlist(upper_list))
  lower_list <- c(as.numeric(as.list(mix$Lower)))
  lower_total <- sum(unlist(lower_list))
  
  for (i in mean_list) {
    x = ((PC_degs[r] * ((i/mean_total)* PC_props[r]))/PC_totdeg)
    PC_Mean_Ratio <- append(PC_Mean_Ratio, x)
  }
  
  
  for (i in upper_list) {
    x = ((PC_degs[r] * ((i/upper_total)* PC_props[r]))/PC_totdeg)
    if(x > 1){
      x <-1
    }
    PC_empirical_Upper_Ratio <- append(PC_empirical_Upper_Ratio, x)
  }
  
  for (i in lower_list) {
    x = ((PC_degs[r] * ((i/lower_total)* PC_props[r]))/PC_totdeg)
    
    PC_empirical_Lower_Ratio <- append(PC_empirical_Lower_Ratio, x)
  }
  
  r = r+1 
}

PC_SD_Ratio = c()
u = 1
for (i in PC_empirical_Lower_Ratio) {
  x = abs(i -PC_Mean_Ratio[u])
  
  PC_SD_Ratio <- append(PC_SD_Ratio, x)
  u = u+1
}

PC_df_empirical_sep<- cbind(PC_df_empirical, PC_Mean_Ratio, PC_SD_Ratio, PC_empirical_Lower_Ratio, PC_empirical_Upper_Ratio)

PC_temp_data <- c(as.numeric(as.list(PC_df_empirical_sep$PC_Mean_Ratio)))
PC_M_sum <- sum(PC_temp_data)
PC_matrix_data<- PC_temp_data/PC_M_sum
PC_df_empirical_sep$PC_Mean_Ratio<- PC_matrix_data
PC_mixing_matrix <- matrix(PC_matrix_data,nrow=3,ncol=3,byrow=TRUE)

#get assortativity based on Newman's definition
PC_assort <- (tr(PC_mixing_matrix)-sum(PC_mixing_matrix %^% 2))/(1-sum(PC_mixing_matrix%^% 2))
PC_mixing_matrix

write.csv(PC_df_empirical_sep, "Mixing_Matrix_data_PC.csv")


# Predicted power for simulation studies
# Kendra Davis-Plourde

# INPUTS
# n: Number of clusters (I)
# l: Number of clinics within each cluster (K)
# m: Number of participants within each clinic (N)
# t: Number of periods (T)
# alpha: Vector of correlations
#         alpha_0=within clinic within period
#         rho_0=between clinic within period
#         rho_1=between clinic between period
#         alpha_1=within clinic between period
#         alpha_2=within-individual correlation
# delta: Effect size 
# bs: baseline time effect
# beta: Vector of period effects
# tot.var: total variance

#############################################################
setwd("/Users/kdavis07/Documents/GitHub/multilevel_crt_samplesize/")
source("VARd.R")
source("study_power.R")

# 1. Gaussian outcome
# define scenarios
scenarios <- read.table("Simulation_Studies/parameters_gaussian_power.txt", header=TRUE, sep="")

for(k in 1:30){
  scenario	<- subset(scenarios, scenario == k)
  n <- scenario$n                         #Number of clusters
  l <- scenario$l                         #Number of clinics
  m <- scenario$m                         #Number of observations within a clinic
  t <- scenario$t                         #Number of time-points
  delta <- scenario$delta
  alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho1,scenario$alpha1)
  #alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho0,scenario$alpha0) # ignoring between-period ICCs (CAC=1)
  tot.var <- scenario$tot.var
  
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="gaussian",alpha=alpha,tot.var=tot.var)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=0.05,df=n-2)
  
  if(k==1){power<-power0} else{power<-c(power,power0)}
}
write.table(power, file="/Users/kdavis07/Dropbox/SW-CRT Methods Development/1_Multilevel_SS/RCode/Simulations/PredictedPower_Gaussian.txt", sep="\t", row.names=F)


# 2. Binomial outcome
# define scenarios
scenarios <- read.table("Simulation_Studies/parameters_binomial_power.txt", header=TRUE, sep="")

for(k in 1:30){
  scenario	<- subset(scenarios, scenario == k)
  n <- scenario$n                         #Number of clusters
  l <- scenario$l                         #Number of clinics
  m <- scenario$m                         #Number of observations within a clinic
  t <- scenario$t                         #Number of time-points
  delta <- log(scenario$exp.delta)
  bs <-log(scenario$bs/(1-scenario$bs)) #scenarios$bs
  beta <- cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3),-0.1/(2^4),-0.1/(2^5)))[1:t]
  alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho1,scenario$alpha1)
  #alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho0,scenario$alpha0) # ignoring between-period ICCs (CAC=1)
  
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=1)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=0.05,df=n-2)
  
  if(k==1){power<-power0} else{power<-c(power,power0)}
}
write.table(power, file="/Users/kdavis07/Dropbox/SW-CRT Methods Development/1_Multilevel_SS/RCode/Simulations/PredictedPower_Binomial.txt", sep="\t", row.names=F)






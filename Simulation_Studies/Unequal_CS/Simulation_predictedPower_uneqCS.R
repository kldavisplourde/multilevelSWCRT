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
# CV.l: coefficient of variation at the subcluster level
# CV.m: coefficient of variation at the subject level

#############################################################
setwd("/Users/kdavis07/Documents/GitHub/multilevel_crt_samplesize/")
source("VARd_MC.R")
source("study_power.R")

# 1. Gaussian outcome
# define scenarios
scenarios <- read.table("Simulation_Studies/Unequal_CS/parameters_gaussian_uneqCS_power.txt", header=TRUE, sep="")

for(k in 1:80){
  scenario	<- subset(scenarios, scenario == k)
  n <- scenario$n                         #Number of clusters
  l <- scenario$l                         #Number of clinics
  m <- scenario$m                         #Number of observations within a clinic
  t <- scenario$t                         #Number of time-points
  delta <- scenario$delta
  alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho1,scenario$alpha1)
  tot.var <- scenario$tot.var
  CV.l <- scenario$CV.l
  CV.m <- scenario$CV.m
  
  vard<-VARd_MC(n=n,t=t,l=l,m=m,CV.l=CV.l,CV.m=CV.m,family="gaussian",alpha=alpha,tot.var=tot.var,
                nsims=1000,seed=9375+k)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=0.05,df=n-2)
  
  if(k==1){power<-power0} else{power<-c(power,power0)}
}
write.table(power, file="/Users/kdavis07/Dropbox/SW-CRT Methods Development/1_Multilevel_SS/RCode/Simulations/PredictedPower_uneqCS_Gaussian.txt", sep="\t", row.names=F)


# 2. Binomial outcome
# define scenarios
scenarios <- read.table("Simulation_Studies/Unequal_CS/parameters_binomial_uneqCS_power.txt", header=TRUE, sep="")

for(k in 1:80){
  scenario	<- subset(scenarios, scenario == k)
  n <- scenario$n                         #Number of clusters
  l <- scenario$l                         #Number of clinics
  m <- scenario$m                         #Number of observations within a clinic
  t <- scenario$t                         #Number of time-points
  delta <- log(scenario$exp.delta)
  bs <-log(scenario$bs/(1-scenario$bs)) #scenarios$bs
  beta <- cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3),-0.1/(2^4),-0.1/(2^5)))[1:t]
  alpha<-c(scenario$alpha0,scenario$rho0,scenario$rho1,scenario$alpha1)
  CV.l <- scenario$CV.l
  CV.m <- scenario$CV.m
  
  vard<-VARd_MC(n=n,t=t,l=l,m=m,CV.l=CV.l,CV.m=CV.m,family="binomial",alpha=alpha,delta=delta,beta=beta,phi=1,
             nsims=1000,seed=257+k)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=0.05,df=n-2)
  
  if(k==1){power<-power0} else{power<-c(power,power0)}
}
write.table(power, file="/Users/kdavis07/Dropbox/SW-CRT Methods Development/1_Multilevel_SS/RCode/Simulations/PredictedPower_uneqCS_Binomial.txt", sep="\t", row.names=F)

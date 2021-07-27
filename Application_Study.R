# Application study from "Sample size considerations for stepped wedge designs with multiple levels of clustering" 
# by Kendra Davis-Plourde, Monica Taljaard, and Fan Li
# Code created by Kendra Davis-Plourde

setwd("/Users/kdavis07/Documents/GitHub/multilevel_crt_samplesize")

####################################################################################################################################
#################################### IMPORTING FUNCTIONS ###########################################################################
source("VARd.R")                 #Importing function for generating the variance of the intervention effect 
source("study_power.R")          #Importing function for generating power
source("designEffect_SWCRT.R")   #Importing function for generating design effect under stepped wedge design (continuous outcome only)
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
#################################### LUMBAR IMAGING STUDY (LIRE) ###################################################################
## A. Outcome is continuous, we want to determine the number of participants per subcluster N using variance and power formulas
### INPUTS ###
delta<- -0.1                     #Intervention effect of interest
alpha<-c(0.046,0.04,0.02,0.023,0.1) #ICCs in the following order:
#                                             alpha_0=within subcluster within period
#                                             rho_0=between subcluster within period
#                                             rho_1=between subcluster between period         
#                                             alpha_1=within subcluster between period
#                                             alpha_2=within-individual correlation
n<-100                            #Number of clusters (I)
l<-17                             #Number of subclusters per cluster (K)
#m<-                              #Number of participants per subcluster (N)
t<-6                              #Number of periods (T)
tot.var<-2.5                      #Total variance under continuous outcome
ER1<-0.05                         #Type I error rate for t-test
df<-n-2                           #Degrees of freedom for t-test
###


# Assuming closed-cohort on subcluster and cross-sectional on participant level.
m0<-seq(50,100,by=1)                #Varying number of participants within a subcluster (N) 
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="gaussian",alpha=alpha,tot.var=tot.var)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[28] #87.5%
m0[28] #N = 77 participants per subcluster in order to achieve 87.5% power


# Assuming closed-cohort on all levels.
m0<-seq(50,100,by=1)                #Varying number of participants within a subcluster (N) 
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cohort",family="gaussian",alpha=alpha,tot.var=tot.var)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[23] #87.5%
m0[23] #N = 72 participants per subcluster in order to achieve 87.5% power


# Assuming cross-sectional on all levels.
m0<-seq(50,120,by=1)                #Varying number of participants within a subcluster (N) 
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cross-sectional",indiv="cross-sectional",family="gaussian",alpha=alpha,tot.var=tot.var)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[50] #87.5%
m0[50] #N = 99 participants per subcluster in order to achieve 87.5% power

####################################################################################################################################

## B. Using design effect
### First we need to determine the total number of participants needed under individual randomization, x.
#### INPUTS ####
delta<- -0.1                      #Intervention effect of interest
tot.var<-2.5                      #Total variance under continuous outcome
ER1<-0.05                         #Type I error rate for t-test
df<-98                            #Degrees of freedom for t-test
####

x0<-seq(9000,10000,by=10)       #Varying number of total participants under individual randomization
for(i in 1:length(x0)){
  x<-x0[i]
  vard<-tot.var*(4/x)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(x0,power)
power[87] #87.5%
x0[87] #x = 9860 total number of participants under individual randomization to achieve 87.5% power

### Next we need to generate the design effect
#### INPUTS ####
x<-9860                                     #Total number of participants under individual randomization
alpha<-c(0.046,0.04,0.02,0.023,0.1)         #ICCs (same as above)
t<-6                                        #Number of periods (T)
S<-t-1                                      #Number of steps
c<-1                                        #Number of measures taken at each step
b<-1                                        #Number of measures taken at baseline
l<-17                                       #Number of subclusters per cluster (K)
m<-77                                       #Number of participants per subcluster (N)
####

design.effect<-designEffect_SWCRT(l=l,m=m,t=t,S=S,c=c,b=b,subcluster="cohort",indiv="cross-sectional",alpha=alpha) #Approximately 12.1
x*design.effect/(l*m) #100 clusters are needed to attain 87.5% power

####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
#################################### WASHINGTON EPT TRIAL ##########################################################################

## A. Outcome is binary and we will use the canonical logit link
### INPUTS ###
delta<-log(0.7)                       #Intervention effect of interest on logit link scale
alpha<-c(0.008,0.007,0.0035,0.004,0.2) #ICCs in the following order:
#                                       alpha_0=within subcluster within period
#                                       rho_0=between subcluster within period
#                                       rho_1=between subcluster between period         
#                                       alpha_1=within subcluster between period
#                                       alpha_2=within-individual correlation
n<-24                                 #Number of clusters (I)
l<-5                                  #Number of subclusters per cluster (K)
#m<-                                  #Number of participants per subcluster (N)
t<-5                                  #Number of periods (T)
bs<-log(0.05/(1-0.05))                #Period 1 effect as a function of the prevalence. Here the prevalence is 0.05.
beta<-cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3)))[1:t] #Period effects. Here we assume a slightly decreasing effect.
phi<-1                                #Scale parameter
ER1<-0.05                             #Type I error rate for t-test
df<-n-2                               #Degrees of freedom for t-test
###

m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[42] #89.5% power
m0[42]    #N = 42 participants per subcluster in order to achieve 89.5% power


# Assuming closed-cohort on all levels.
m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cohort",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[66] #89.5% power
m0[66]    #N = 66 participants per subcluster in order to achieve 89.5% power


# Assuming cross-sectional on all levels.
m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cross-sectional",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[42] #89.5% power
m0[42]    #N = 42 participants per subcluster in order to achieve 89.5% power

####################################################################################################################################

## B. Investigating alternative time effects
### Larger decreasing time effect
#### INPUTS ####
beta<-cumsum(c(bs,-1,-1/2,-1/(2^2),-1/(2^3)))[1:t] #Assuming a larger decreasing period effect.
####

m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[139] #89.5% power
m0[139]    #N = 139 participants per subcluster in order to achieve 89.5% power

# Assuming closed-cohort on all levels.
m0<-seq(1,300,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cohort",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[218] #89.5% power
m0[218]    #N = 218 participants per subcluster in order to achieve 89.5% power

# Assuming cross-sectional on all levels.
m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cross-sectional",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[139] #89.5% power
m0[139]    #N = 139 participants per subcluster in order to achieve 89.5% power
######################################

### Smaller decreasing time effect
#### INPUTS ####
beta<-cumsum(c(bs,-0.01,-0.01/2,-0.01/(2^2),-0.01/(2^3)))[1:t] #Assuming a smaller decreasing period effect.
####

m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[37] #89.3% power
m0[37]    #N = 37 participants per subcluster in order to achieve 89.3% power


# Assuming closed-cohort on all levels.
m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cohort",indiv="cohort",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[59] #89.6% power
m0[59]    #N = 59 participants per subcluster in order to achieve 89.6% power

# Assuming cross-sectional on all levels.
m0<-seq(1,200,by=1)                   #Varying number of participants within a subcluster (N)
for(i in 1:length(m0)){
  m<-m0[i]
  vard<-VARd(n=n,l=l,m=m,t=t,subcluster="cross-sectional",indiv="cross-sectional",family="binomial",alpha=alpha,delta=delta,beta=beta,phi=phi)
  power0<-study_power(delta=delta,var.delta=vard,typeI.error=ER1,df=df)
  
  if(i==1) {
    power<-power0
  } else {
    power<-c(power,power0)
  }
}
plot(m0,power)
power[37] #89.3% power
m0[37]    #N = 37 participants per subcluster in order to achieve 89.3% power

####################################################################################################################################
####################################################################################################################################
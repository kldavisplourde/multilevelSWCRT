# Application results from "Sample size calculation for stepped wedge cluster randomized trials with multiple levels of clustering." 

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
delta<-log(0.95)                  #Intervention effect of interest
alpha<-c(0.013,0.00325,0.001625,0.0065,0.1) #ICCs in the following order:
#                                             alpha_0=within subcluster within period
#                                             rho_0=between subcluster within period
#                                             rho_1=between subcluster between period         
#                                             alpha_1=within subcluster between period
#                                             alpha_2=within-individual correlation
n<-100                            #Number of clusters (I)
l<-17                             #Number of subclusters per cluster (L)
#m<-                              #Number of participants per subcluster (N)
t<-6                              #Number of periods (T)
tot.var<-1                        #Total variance under continuous outcome
ER1<-0.05                         #Type I error rate for t-test
df<-n-2                           #Degrees of freedom for t-test
###

m0<-seq(1,10,by=1)                #Varying number of participants within a subcluster (N)
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
power[3] #81.9%
m0[3] #N = 3 participants per subcluster in order to achieve 81.9% power

####################################################################################################################################

## B. Using design effect
### First we need to determine the total number of participants needed under individual randomization, x.
#### INPUTS ####
delta<-log(0.95)                  #Intervention effect of interest
tot.var<-1                        #Total variance under continuous outcome
ER1<-0.05                         #Type I error rate for t-test
df<-98                            #Degrees of freedom for t-test
####

x0<-seq(10000,20000,by=250)       #Varying number of total participants under individual randomization
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
power[10] #80.2%
x0[10] #x = 12250 total number of participants under individual randomization to achieve 80.2% power

### Next we need to generate the design effect
#### INPUTS ####
x<-12250                                    #Total number of participants under individual randomization
alpha<-c(0.013,0.00325,0.001625,0.0065,0.1) #ICCs (same as above)
t<-6                                        #Number of periods (T)
S<-t-1                                      #Number of steps
c<-1                                        #Number of measures taken at each step
b<-1                                        #Number of measures taken at baseline
l<-17                                       #Number of subclusters per cluster (L)
m<-3                                        #Number of participants per subcluster (N)
####

design.effect<-designEffect_SWCRT(l=l,m=m,t=t,S=S,c=c,b=b,subcluster="cohort",indiv="cross-sectional",alpha=alpha) #Approximately 0.40
x*design.effect/(l*m) #96 clusters are needed to attain 80.2% power

####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
#################################### WASHINGTON EPT TRIAL ##########################################################################

## A. Outcome is binary and we will use the canonical logit link
### INPUTS ###
delta<-log(0.7)                       #Intervention effect of interest on logit link scale
alpha<-c(0.007,0.002,0.001,0.004)     #ICCs in the following order:
#                                       alpha_0=within subcluster within period
#                                       rho_0=between subcluster within period
#                                       rho_1=between subcluster between period         
#                                       alpha_1=within subcluster between period
#                                       alpha_2=within-individual correlation
n<-24                                 #Number of clusters (I)
l<-9                                  #Number of subclusters per cluster (L)
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
power[16] #81.6% power
m0[16]    #N = 16 participants per subcluster in order to achieve 81.6% power

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
power[52] #80.7% power
m0[52]    #N = 52 participants per subcluster in order to achieve 80.7% power
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
power[14] #81.1% power
m0[14]    #N = 14 participants per subcluster in order to achieve 81.1% power

####################################################################################################################################
####################################################################################################################################
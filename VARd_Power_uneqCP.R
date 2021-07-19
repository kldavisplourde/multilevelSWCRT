########################################################################################################################
# Monte Carlo function for generating variance and power estimates for multilevel SWCRT with closed-cohort at subcluster level and cross-sectional at subject level design (design variant B) 
# with unequal cluster-period sizes (between-cluster imbalance). 
# Code based on SWD_GEEPower_cp_variable_v2.R by Zibo Tian, modified by Kendra Davis-Plourde.

# INPUT
# n: Number of clusters (I)
# t: Number of periods (T)
# l: Average number of subclusters per cluster (K)
# m: Average number of participants per subcluster (N)
# CV.l: Degree of between-cluster imbalance with respect to the number of subclusters measured by coefficient of variation
# CV.m: Degree of between-subcluster imbalance measured by coefficient of variation
# family: "gaussian" or "binomial" (gaussian with identity link or binomial with logit link)
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
# delta: Effect size on the link function scale
# beta: Vector of period effect on the link function scale 
# phi: Common dispersion (only needed when outcome is binary)
# tot.var: total variance of the outcome (only needed when outcome is continuous)
# typeI.error: Type I error rate for t-test (default is 0.05)
# df: degrees of freedom for t-test
# nsims: Number of simulation runs to compute the variance of treatment effect estimator
# seed: a seed to control reproducible output

# Output: An array of power estimated under both z-test and t-test
########################################################################################################################

setwd("/Users/kdavis07/Documents/GitHub/multilevel_crt_samplesize/")
source("study_power.R")          #Importing function for generating power

VARd_Power_uneqCP <- function(n, t, l, m, CV.l=0, CV.m=0, family="gaussian", alpha, delta, beta=rep(0, t), phi=1, tot.var=1, typeI.error=0.05, df=n-2, nsims=1000, seed=2021){
  
  # elements of efficiency calculation
  scheme<-rep(n/(t-1),t-1)
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]
  if(sum(scheme) != n){stop("invalid randomization scheme")}
  
  
  # function to account for between-cluster imbalance 
  l.variable <- function(n, t, CV.l, l) {
    if (CV.l == 0) {
      return(rep(l,n))
    } else {
      Nsubc.raw <- rgamma(n, 1/CV.l^2, rate = 1/(l * CV.l^2))
      Nsubc <- as.integer(Nsubc.raw * (l/mean(Nsubc.raw)))
      Nsubc[Nsubc < 2] = 2
      return(Nsubc)
    }
  }
  
  m.variable <- function(n, t, CV.m, m) {
    if (CV.m == 0) {
      return(rep(m,n))
    } else {
      N.raw <- rgamma(n, 1/CV.m^2, rate = 1/(m * CV.m^2))
      N <- as.integer(N.raw * (m/mean(N.raw)))
      N[N < 3] = 3
      return(N)
    }
  }
  
  
  # Generate correlation matrix (for Gaussian)
  R <- function(ki, mi) {
    B <- (1 + (mi-1)*alpha[1] + mi*(ki-1)*alpha[2])/(ki*mi)
    C <- (alpha[4] + (mi-1)*alpha[4] + mi*(ki-1)*alpha[3])/(ki*mi)
    
    R <- kronecker(diag(t), B-C) + kronecker(matrix(1,nrow=t,ncol=t),C)
      
    return(R)
  }
  
  # Generate variance components using given correlations (for binomial)
  sig2_e<-pi^2/3                                                        #variance of standard logistic distribution
  sig2_b<-sig2_e*alpha[3]/(1-alpha[1])                                  #within cluster
  sig2_c<-sig2_e*(alpha[4]-alpha[3])/(1-alpha[1])                       #within subcluster
  sig2_s<-sig2_e*(alpha[2]-alpha[3])/(1-alpha[1])                       #within cluster within period
  sig2_p<-sig2_e*(alpha[1]-alpha[4]-alpha[2]+alpha[3])/(1-alpha[1])     #within subcluster within period
  sig2_g<-0                                                             #within person
  
  
  # calculate variance
  for (s in 1:nsims) {
    vardelta <- NULL
    meanCP.size <- NULL
    totalSample <- NULL
    
    set.seed(seed + s)
    #simulated baseline cluster sizes under given CV
    l_var <- l.variable(n, t, CV.l, l)
    m_var <- m.variable(n, t, CV.m, m)
    CP.size <- mean(l_var*m_var)
    
    # elements of power calculation
    Omega <- matrix(0,t+1,t+1)
    
    size <- 0
      
    for (i in 1:n){    # loop through each unique cluster
      
      # design matrix
      # put trt indicator at the first column for convenience
      period <- rep(1:t)
      X <- trtSeq[i,]
      Z <- cbind(X,diag(t))
      
      invRi <- solve(R(l_var[i], m_var[i]))
      
      if (min(eigen(R(l_var[i], m_var[i]))$values)<0){
        stop(paste0("Correlation Structure of cluster ", i, " is Not Positive Definite"))
      }
      
      size <- size + l_var[i]*m_var[i]
        
      # unique marginal means
      if (family=="gaussian"){             # gaussian with identity link
          mu <- c(Z%*%c(delta,beta))
          Omega <- Omega + (t(Z)%*%invRi%*%Z)/tot.var
       } else if (family=="binomial"){     # binomial with logit link
          gmu<-c(Z%*%c(delta,beta))
          E.Gprime<-as.vector(2+exp(0.5*(sig2_b+sig2_c+sig2_s+sig2_p+sig2_g))*(exp(-gmu)+exp(gmu)))
           
          V<-(E.Gprime/(l_var[i]*m_var[i])+sig2_p/l_var[i]+sig2_s)*diag(t)+(sig2_b+sig2_c/l_var[i]+sig2_g/(l_var[i]*m_var[i]))*matrix(1,nrow=t,ncol=t)
          invV<-solve(V)
           
          Omega<-Omega + (t(Z)%*%invV%*%Z)/phi
         }
      }
      
      # delta <- c(solve(Omega,Ry))[1]   # check whether this is close to the assumed delta
      vardelta <- c(vardelta, solve(Omega)[1,1])    # var of the trt effect estimator
      meanCP.size <- c(meanCP.size, CP.size)
      totalSample <- c(totalSample, size)
    #print(vardelta)
    
    # t-test
      t1_power <- study_power(delta=delta, var.delta=mean(vardelta), typeI.error=typeI.error, df=df)
  }
  
  return(data.frame(var.delta=mean(vardelta), power=t1_power, CP.size=meanCP.size, total.size=totalSample))
}



# From application study:
delta<- -0.1                     #Intervention effect of interest
alpha<-c(0.046,0.04,0.02,0.023,0.1) #ICCs in the following order:
#                                             alpha_0=within subcluster within period
#                                             rho_0=between subcluster within period
#                                             rho_1=between subcluster between period         
#                                             alpha_1=within subcluster between period
#                                             alpha_2=within-individual correlation
n<-100                            #Number of clusters (I)
l<-17                             #Number of subclusters per cluster (K)
m<-70                              #Number of participants per subcluster (N)
t<-6                              #Number of periods (T)
tot.var<-2.5                      #Total variance under continuous outcome
ER1<-0.05                         #Type I error rate for t-test
VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=0, CV.m=0, family="gaussian", alpha=alpha, delta=delta, tot.var=tot.var, typeI.error=ER1, df=n-2, nsims=100, seed=2021)
# Matches our application study!!! We get the same variance and power estimates!

#Testing with CVs
VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=0.25, CV.m=0.75, family="gaussian", alpha=alpha, delta=delta, tot.var=tot.var, typeI.error=ER1, df=n-2, nsims=100, seed=2021)


# From application study
delta<-log(0.7)                       #Intervention effect of interest on logit link scale
alpha<-c(0.008,0.007,0.0035,0.004,0.1) #ICCs in the following order:
#                                       alpha_0=within subcluster within period
#                                       rho_0=between subcluster within period
#                                       rho_1=between subcluster between period         
#                                       alpha_1=within subcluster between period
#                                       alpha_2=within-individual correlation
n<-24                                 #Number of clusters (I)
l<-5                                  #Number of subclusters per cluster (K)
m<-42                                  #Number of participants per subcluster (N)
t<-5                                  #Number of periods (T)
bs<-log(0.05/(1-0.05))                #Period 1 effect as a function of the prevalence. Here the prevalence is 0.05.
beta<-cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3)))[1:t] #Period effects. Here we assume a slightly decreasing effect.
phi<-1                                #Scale parameter
ER1<-0.05                             #Type I error rate for t-test
df<-n-2                               #Degrees of freedom for t-test
VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=0, CV.m=0, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=phi, typeI.error=ER1, df=n-2, nsims=100, seed=2021)
# Matches our application study!!! We get the same variance and power estimates!

#Testing with CVs
VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=0.25, CV.m=0.25, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=phi, typeI.error=ER1, df=n-2, nsims=100, seed=2021)


# SE and Power for simulation study
scenarios.Gaus <- read.table("Simulation_Studies/parameters_gaussian_power.txt", header=TRUE, sep="")

CV.l<-0;CV.m<-0;seed<-7735        # this should match our original predicted power
CV.l<-0;CV.m<-0.25;seed<-8393 
CV.l<-0;CV.m<-0.75;seed<-753 
CV.l<-0.25;CV.m<-0;seed<-237 
CV.l<-0.75;CV.m<-0;seed<-193 
CV.l<-0.25;CV.m<-0.25;seed<-6299 
CV.l<-0.25;CV.m<-0.75;seed<-2734 
CV.l<-0.75;CV.m<-0.25;seed<-3704 
CV.l<-0.75;CV.m<-0.75;seed<-53401

Gaussian.results<- matrix(0,nrow(scenarios.Gaus),4)
for(i in 1:nrow(scenarios.Gaus)){
  scenarios <- subset(scenarios.Gaus, scenario == i)
  scenario	<- i
  n <- scenarios$n                         #Number of clusters
  l <- scenarios$l                         #Number of clinics
  m <- scenarios$m                         #Number of observations within a clinic
  t <- scenarios$t                         #Number of time-points
  delta <- scenarios$delta
  bs <- scenarios$bs
  beta <- cumsum(c(bs,0.1,0.1*0.5,0.1*(0.5^2),0.1*(0.5^3),0.1*(0.5^4),0.1*(0.5^5)))[1:t]
  alpha<-c(scenarios$alpha0,scenarios$rho0,scenarios$rho1,scenarios$alpha1)
  tot.var <- scenarios$tot.var
  
  it.k<-VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=CV.l, CV.m=CV.m, family="gaussian", alpha=alpha, delta=delta, tot.var=tot.var, typeI.error=0.05, df=n-2, nsims=100, seed=seed+i)
  
  Gaussian.results[i,1] <- as.numeric(sqrt(it.k[1]))
  Gaussian.results[i,2] <- as.numeric(it.k[2])
  Gaussian.results[i,3] <- as.numeric(it.k[3])
  Gaussian.results[i,4] <- as.numeric(it.k[4])
}
Gaussian.results



scenarios.Bin <- read.table("Simulation_Studies/parameters_binomial_power.txt", header=TRUE, sep="")

CV.l<-0;CV.m<-0;seed<-7735        # this should match our original predicted power
CV.l<-0;CV.m<-0.25;seed<-8393
CV.l<-0;CV.m<-0.75;seed<-753
CV.l<-0.25;CV.m<-0;seed<-237
CV.l<-0.75;CV.m<-0;seed<-193
CV.l<-0.25;CV.m<-0.25;seed<-6299 
CV.l<-0.25;CV.m<-0.75;seed<-2734
CV.l<-0.75;CV.m<-0.25;seed<-3704
CV.l<-0.75;CV.m<-0.75;seed<-53401

Binomial.results<- matrix(0,nrow(scenarios.Bin),2)
for(i in 1:nrow(scenarios.Bin)){
  scenarios <- subset(scenarios.Bin, scenario == i)
  scenario	<- i
  n <- scenarios$n                         #Number of clusters
  l <- scenarios$l                         #Number of clinics
  m <- scenarios$m                         #Number of observations within a clinic
  t <- scenarios$t                         #Number of time-points
  delta <- log(scenarios$exp.delta)
  bs <-log(scenarios$bs/(1-scenarios$bs)) #scenarios$bs
  beta <- cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3),-0.1/(2^4),-0.1/(2^5)))[1:t]
  alpha<-c(scenarios$alpha0,scenarios$rho0,scenarios$rho1,scenarios$alpha1)
  
  it.k<-VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=CV.l, CV.m=CV.m, family="binomial", alpha=alpha, delta=delta, phi=1, typeI.error=0.05, df=n-2, nsims=1000, seed=seed+i)
  
  Binomial.results[i,1] <- as.numeric(sqrt(it.k[1]))
  Binomial.results[i,2] <- as.numeric(it.k[2])
}
Binomial.results


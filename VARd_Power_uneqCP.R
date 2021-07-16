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

n<-8                            #Number of clusters (I)
t<-3                            #Number of periods (T)
l<-4                            # Average number of subclusters per cluster (K)
m<-5
CV.l<-0.25
CV.m<-0.25
family<-"gaussian"
alpha<-c(0.046,0.04,0.02,0.023) #ICCs in the following order:
#                                             alpha_0=within subcluster within period
#                                             rho_0=between subcluster within period
#                                             rho_1=between subcluster between period         
#                                             alpha_1=within subcluster between period
delta<-0.1
tot.var<-1                      #Total variance under continuous outcome
typeI.error<-0.05                         #Type I error rate for t-test
df<-n-2                           #Degrees of freedom for t-test
nsims<-100
seed<-2764

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
  
  
  R <- function(ki, mi) {
    B <- (1 + (mi-1)*alpha[1] + mi*(ki-1)*alpha[2])/(ki*mi)
    C <- (alpha[4] + (mi-1)*alpha[4] + mi*(ki-1)*alpha[3])/(ki*mi)
    
    R <- kronecker(diag(t), B-C) + kronecker(matrix(1,nrow=t,ncol=t),C)
      
    return(R)
  }
  
  
  # calculate variance
  for (s in 1:nsims) {
    vardelta <- NULL
    
    set.seed(seed + s)
    #simulated baseline cluster sizes under given CV
    l_var <- l.variable(n, t, CV.l, l)
    m_var <- m.variable(n, t, CV.m, m)
    
    # elements of power calculation
    Omega <- matrix(0,t+1,t+1)
      
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
        
      # unique marginal means
      if (family=="gaussian"){             # gaussian with identity link
          mu <- c(Z%*%c(delta,beta))
          Omega <- Omega + (t(Z)%*%invRi%*%Z)/tot.var
       } else if (family=="binomial"){     # binomial with logit link
          gmu <- c(X%*%c(delta,beta))
          mu <- plogis(gmu)
          W <- diag(sqrt(mu*(1-mu))) %*% invRi %*% diag(sqrt(mu*(1-mu)))
          Omega <- Omega + (t(X)%*%W%*%X)/phi
        }
      }
      
      # delta <- c(solve(Omega,Ry))[1]   # check whether this is close to the assumed delta
      vardelta <- c(vardelta, solve(Omega)[1,1])    # var of the trt effect estimator
    
    #print(vardelta)
    
    # t-test
      t1_power <- study_power(delta=delta, var.delta=mean(vardelta), typeI.error=typeI.error, df=df)
  }
  
  return(data.frame(var.delta=mean(vardelta), power=t1_power))
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

VARd_Power_uneqCP(n=n, t=t, l=l, m=m, CV.l=0.25, CV.m=0.75, family="gaussian", alpha=alpha, delta=delta, tot.var=tot.var, typeI.error=ER1, df=n-2, nsims=100, seed=2021)



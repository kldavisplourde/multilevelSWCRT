########################################################################################################################
# Monte Carlo approximation of the variance of the intervention effect for multilevel stepped wedge cluster randomized trials
# with closed-cohort at subcluster level and cross-sectional at subject level design (design variant B) 
# and equal or unequal cluster-sizes (between-cluster imbalance only). 
# Code based on SWD_GEEPower_cp_variable_v2.R by Zibo Tian, modified by Kendra Davis-Plourde.

# INPUT
# n: Number of clusters (I)
# l: Average number of subclusters per cluster (K)
# m: Average number of participants per subcluster (N)
# t: Number of periods (T)
# CV.l: Degree of between-cluster imbalance with respect to the number of subclusters measured by coefficient of variation
# CV.m: Degree of between-subcluster imbalance measured by coefficient of variation
# family: "gaussian" or "binomial" (gaussian with identity link or binomial with logit link)
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
# delta: Effect size on the link function scale
# beta: Vector of period effects on the link function scale 
# phi: Common dispersion (only needed when outcome is binomial)
# tot.var: total variance of the outcome (only needed when outcome is gaussian)
# nsims: Number of simulation runs to compute the variance of treatment effect estimator
# seed: a seed to control reproducible output
########################################################################################################################

VARd_uneqCS <- function(n, l, m, t, CV.l=0, CV.m=0, family=c("gaussian","binomial"), alpha, delta, beta=rep(0, t), phi=1, tot.var=1, nsims=1000, seed=2021){
  
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
    vard <- NULL
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
      
      # delta <- c(solve(Omega,Ry))[1]
      vard <- c(vard, solve(Omega)[1,1])    # var of the trt effect estimator
      meanCP.size <- c(meanCP.size, CP.size)
      totalSample <- c(totalSample, size)
  }
  
  vardelta=mean(vard)
  #return(data.frame(vardelta=mean(vard), CP.size=meanCP.size, total.size=totalSample))
  return(vardelta)
}

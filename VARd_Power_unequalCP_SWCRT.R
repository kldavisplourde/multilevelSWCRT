########################################################################################################################
# Monte Carlo function for generating power for multilevel cross-sectional SWCRTs with unequal cluster-period sizes (between-cluster imbalance). 
# Code based on SWD_GEEPower_cp_variable_v2.R by Zibo Tian, modified by Kendra Davis-Plourde.

# INPUT
# n: Number of clusters (I)
# t: Number of periods (T)
# KN: Average number of individuals at each cluster-period (K*N)
# CV: Degree of between-cluster imbalance measured by coefficient of variation
# family: "gaussian" or "binomial" (gaussian with identity link or binomial with logit link)
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
#         alpha_2=within-individual correlation
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

swdgeepower_var <- function(n, t, KN, CV=0, family="gaussian", alpha, delta, beta=rep(0, t), phi=1, tot.var=1, typeI.error=0.05, df=n-2, nsims=1000, seed=2021){
  
  # elements of efficiency calculation
  scheme<-rep(n/(t-1),t-1)
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]
  if(sum(scheme) != n){stop("invalid randomization scheme")}
  
  
  # function to account for between-cluster imbalance 
  n.variable <- function(n, t, CV, KN) {
    if (CV == 0) {
      return(rep(KN,n))
    } else {
      N.raw <- rgamma(n, 1/CV^2, rate = 1/(KN * CV^2))
      N <- as.integer(N.raw * (KN/mean(N.raw)))
      N[N < 5] = 5
      return(rep(N,n))
    }
  }
  
  
  R <- function(ni) {
    R <- diag(t)*((1+(ni-1)*alpha0)/ni-alpha1) + matrix(1,t,t)*alpha1
      
    return(R)
  }
  
  
  # calculate variance
  for (s in 1:nsims) {
    vardelta <- NULL
    
    set.seed(seed + s)
    #simulated baseline cluster sizes under given CV
    n_var <- n.variable(n, t, CV, KN)
    
    # elements of power calculation
    Omega <- matrix(0,t+1,t+1)
      
    for (i in 1:itrt){    # loop through each unique cluster
        
      # design matrix
      # put trt indicator at the first column for convenience
      period <- rep(1:t)
      X <- rep(trtSeq[i,])
      for (j in 1:t){
        X <- cbind(X, as.numeric(period==j))
      }
        
      invRi <- solve(R(n_var[i]))
        
      if (min(eigen(R(n_var[i]))$values)<0){
        stop(paste0("Correlation Structure of cluster ", i, " is Not Positive Definite"))
      }
        
      # unique marginal means
      if (family=="gaussian"){             # gaussian with identity link
          mu <- c(X%*%c(delta,beta))
          W <- invRi
       } else if (family=="binomial"){     # binomial with logit link
            gmu <- c(X%*%c(delta,beta))
            mu <- plogis(gmu)
            W <- diag(sqrt(mu*(1-mu))) %*% invRi %*% diag(sqrt(mu*(1-mu)))
        }
        
        Omega <- Omega + (t(X)%*%W%*%X)/phi
      }
      
      # delta <- c(solve(Omega,Ry))[1]   # check whether this is close to the assumed delta
      vardelta <- c(vardelta, solve(Omega)[1,1])    # var of the trt effect estimator
    
    #print(vardelta)
    
    # t-test with df = # clusters minus 2
    df1 <- df
    t1_power <- pt(qt(typeI.error/2,df=df1) + abs(delta)/sqrt(mean(vardelta)), df=df1)
  }
  
  return(data.frame(t_test_df_cp=t1_power))
}


# ---------------------------------------------------------------------------------

# examples
n = 8; t = 9; tn = 6; KN = 100
CV=0

tnew = t+tn
design <- matrix(0,t-1,t)
design[upper.tri(design)] <- 1
design <- cbind(design, matrix(1,t-1,tnew-t))
bas = 0.25
end = 0.27   # time trend at the last period
p=seq(bas,end,length.out=t+tn)
beta = log(p/(1-p))
delta =  log(0.35/(1-0.35)) - beta[tnew]
alpha = c(0.04, 0.02)
swdgeepower_var(n=n, t=t+tn, KN=KN, CV=CV, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=1, typeI.error=0.025)
swdgeepower_var(n=n, t=t+tn, KN=KN, CV=CV, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=1, typeI.error=0.025)

swdgeepower_var(n=n, t=t+tn, KN=KN, CV=1, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=1, typeI.error=0.025)
swdgeepower_var(n=n, t=t+tn, KN=KN, CV=1, family="binomial", alpha=alpha, delta=delta, beta=beta, phi=1, typeI.error=0.025)


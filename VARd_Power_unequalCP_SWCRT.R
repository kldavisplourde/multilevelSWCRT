########################################################################################################################
# Function of power calculations for cross-sectional SWD with different responses, using cluster-period approach, 
# when between-cluster imbalance presents (Version 2: MC calculation of average variance of treatment effect 
# estimator and then a one-time calclulation of power)

# INPUT
# I: Number of clusters
# J: Number of periods
# K: Average number of individuals at each cluster-period cell
# CV: Degree of between-cluster imbalance measured by coefficient of variation
# design: Data set that describes the study design
# family: "gaussian", "binomial", or "poisson"
# link: "identity", "logit", or "log"
#       (default: link = "identity" for family = "gaussian";
#                 link = "logit" for family = "binomial";
#                 link = "log" for family = "poisson")
# windep: FALSE (default) for model-based varaince,
#         TRUE for robust sandwich varaince with an independence working correlation matrix
# corstr: "nested exchangeable" for nested exchangeable correlation structure (cross-sectional),
#         "exponential decay" for exponential decay correlation structure (cross-sectional)
# delta: Effect size on the link function scale
# beta: Vector of period effect on the link function scale
# phi: Dispersion parameter (for gaussian, it's marginal variance of the outcome)
# size: Type I error rate
# alpha: correlations
#        (within-time, between-time) for nested exchangeable,
#        (within-time, decay) for exponential decay,
# df: degrees of freedom for t-test
# nsims: Number of simulation runs to compute the variance of treatment effect estimator
# seed: a seed to control reproduceable output

# Output: An array of power estimated under both z-test and t-test
########################################################################################################################

swdgeepower_var <- function(I, J, K, CV=0, design=NULL, family="gaussian", link=NULL, windep=FALSE, corstr, delta, beta=rep(0, J), phi=1, size=0.05, alpha, df=I-2, nsims=1000, seed=2021){
  
  #check design
  if (nrow(design)==I & ncol(design)==J){
    trtSeq <- design
    itrt <- I
  } else {
    stop("Design is Not Consistent with Cluster/Period Number")
  }
  
  # link
  if (is.null(link)){
    if (family=="gaussian"){
      link <- "identity"
    } else if (family=="binomial"){
      link <- "logit"
    } else if (family=="poisson"){
      link <- "log"
    }
  }
  
  # function to generate simple AR1 correlation
  ar1 <- function(n, rho) {
    exponent <- abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE) - (1:n-1))
    res <- rho^exponent
    return(res)
  }
  
  # function to account for between-cluster imbalance 
  n.variable <- function(I, J, CV, K) {
    if (CV == 0) {
      return(rep(K,I))
    } else {
      N.raw <- rgamma(I, 1/CV^2, rate = 1/(K * CV^2))
      N <- as.integer(N.raw * (K/mean(N.raw)))
      N[N < 5] = 5
      return(rep(N,I))
    }
  }
  
  # true correlation matrix
  if (corstr=="nested exchangeable"){
    alpha0 <- alpha[1]
    alpha1 <- alpha[2]
  } else if (corstr=="exponential decay"){
    tau <- alpha[1]
    rho <- alpha[2]
  } else {
    stop("Correlation Structure is Misspecified")
  }
  
  
  
  R <- function(ni) {
    if (corstr=="nested exchangeable"){
      R <- diag(J)*((1+(ni-1)*alpha0)/ni-alpha1) + matrix(1,J,J)*alpha1
    } else if (corstr=="exponential decay"){
      R <- diag(J)*((1+(ni-1)*tau)/ni) + (ar1(J, rho)-diag(J))*tau
    }
    return(R)
  }
  
  
  # calculate variance
  for (s in 1:nsims) {
    vardelta <- NULL
    
    set.seed(seed + s)
    #simulated baseline cluster sizes under given CV
    n_var <- n.variable(I, J, CV, K)
    
    
    if (windep==FALSE){    # model-based varaince
      
      # elements of power calculation
      Omega <- matrix(0,J+1,J+1)
      # Ry <- matrix(0,J+1,1)
      
      for (i in 1:itrt){    # loop through each unique cluster
        
        # design matrix
        # put trt indicator at the first column for convenience
        period <- rep(1:J)
        X <- rep(trtSeq[i,])
        for (j in 1:J){
          X <- cbind(X, as.numeric(period==j))
        }
        
        invRi <- solve(R(n_var[i]))
        
        #if (i == 3 & s == 500) print(R(n_var[i]))
        
        if (min(eigen(R(n_var[i]))$values)<0){
          stop(paste0("Correlation Structure of cluster ", i, " is Not Positive Definite"))
        }
        
        # unique marginal means
        if (family=="gaussian"){
          if (link=="identity"){    # gaussian with identity link
            mu <- c(X%*%c(delta,beta))
            W <- invRi
          } else if (link=="log"){    # gaussian with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(mu) %*% invRi %*% diag(mu)
          }
        } else if (family=="binomial"){
          if (link=="identity"){    # binomial with identity link
            mu <- c(X%*%c(delta,beta))
            W <- diag(sqrt(1/(mu*(1-mu)))) %*% invRi %*% diag(sqrt(1/(mu*(1-mu))))
          } else if (link=="logit"){    # binomial with logit link
            gmu <- c(X%*%c(delta,beta))
            mu <- plogis(gmu)
            W <- diag(sqrt(mu*(1-mu))) %*% invRi %*% diag(sqrt(mu*(1-mu)))
          } else if (link=="log"){    # binomial with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(sqrt(mu/(1-mu))) %*% invRi %*% diag(sqrt(mu/(1-mu)))
          }
        } else if (family=="poisson"){
          if (link=="identity"){    # poisson with identity link
            mu <- c(X%*%c(delta,beta))
            W <- diag(sqrt(1/mu)) %*% invRi %*% diag(sqrt(1/mu))
          } else if (link=="log"){    # poisson with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(sqrt(mu)) %*% invRi %*% diag(sqrt(mu))
          }
        }
        
        Omega <- Omega + (t(X)%*%W%*%X)/phi
        # Ry <- Ry + ctrt*(t(X)%*%W%*%gmu)
      }
      #if (s %in% 490:500) print(solve(Omega)[1,1])
      
      # delta <- c(solve(Omega,Ry))[1]   # check whether this is close to the assumed delta
      vardelta <- c(vardelta, solve(Omega)[1,1])    # var of the trt effect estimator
    } else if (windep==TRUE){    # robust sandwich varaince with an independence working correlation matrix
      
      # elements of power calculation
      Omega0 <- matrix(0,J+1,J+1)
      Omega1 <- matrix(0,J+1,J+1)
      
      #Independence working correlation 
      invI <- function(ni) {
        return(diag(J)*ni)
      }
      
      for (i in 1:itrt){    # loop through each unique cluster
        
        invIi <- invI(n_var[i])
        Ri <- R(n_var[i])
        
        # design matrix
        # put trt indicator at the first column for convenience
        period <- rep(1:J)
        X <- rep(trtSeq[i,])
        for (j in 1:J){
          X <- cbind(X, as.numeric(period==j))
        }
        
        # unique marginal means
        if (family=="gaussian"){
          if (link=="identity"){    # gaussian with identity link
            mu <- c(X%*%c(delta,beta))
            W <- invIi
            W1 <- W
            V <- Ri
          } else if (link=="log"){    # gaussian with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(mu) %*% invIi %*% diag(mu)
            W1 <- diag(mu) %*% invIi
            V <- Ri
          }
        } else if (family=="binomial"){
          if (link=="identity"){    # binomial with identity link
            mu <- c(X%*%c(delta,beta))
            W <- diag(sqrt(1/(mu*(1-mu)))) %*% invIi %*% diag(sqrt(1/(mu*(1-mu))))
            W1 <- W
            V <- diag(sqrt(mu*(1-mu))) %*% Ri %*% diag(sqrt(mu*(1-mu)))
          } else if (link=="logit"){    # binomial with logit link
            gmu <- c(X%*%c(delta,beta))
            mu <- plogis(gmu)
            W <- diag(sqrt(mu*(1-mu))) %*% invIi %*% diag(sqrt(mu*(1-mu)))
            W1 <- diag(sqrt(mu*(1-mu))) %*% invIi %*% diag(sqrt(1/(mu*(1-mu))))
            V <- diag(sqrt(mu*(1-mu))) %*% Ri %*% diag(sqrt(mu*(1-mu)))
          } else if (link=="log"){    # binomial with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(sqrt(mu/(1-mu))) %*% invIi %*% diag(sqrt(mu/(1-mu)))
            W1 <- diag(sqrt(mu/(1-mu))) %*% invIi %*% diag(sqrt(1/(mu*(1-mu))))
            V <- diag(sqrt(mu*(1-mu))) %*% Ri %*% diag(sqrt(mu*(1-mu)))
          }
        } else if (family=="poisson"){
          if (link=="identity"){    # poisson with identity link
            mu <- c(X%*%c(delta,beta))
            W <- diag(sqrt(1/mu)) %*% invIi %*% diag(sqrt(1/mu))
            W1 <- W
            V <- diag(sqrt(mu)) %*% Ri %*% diag(sqrt(mu))
          } else if (link=="log"){    # poisson with log link
            gmu <- c(X%*%c(delta,beta))
            mu <- exp(gmu)
            W <- diag(sqrt(mu)) %*% invIi %*% diag(sqrt(mu))
            W1 <- diag(sqrt(mu)) %*% invIi %*% diag(sqrt(1/mu))
            V <- diag(sqrt(mu)) %*% Ri %*% diag(sqrt(mu))
          }
        }
        
        Omega1 <- Omega1 + (t(X)%*%W%*%X)/phi
        Omega0 <- Omega0 + (t(X)%*%W1%*%V%*%t(W1)%*%X)/phi
      }
      
      invOmega1 <- solve(Omega1)
      vardelta <- c(vardelta, (invOmega1%*%Omega0%*%invOmega1)[1,1])    # var of the trt effect estimator
    } else {
      stop("Some Parameter is Misspecified")
    }
    
    #print(vardelta)
    
    # z-test power
    z_power <- pnorm(qnorm(size/2) + abs(delta)/sqrt(mean(vardelta)))
    
    # t-test with df = # clusters minus 2
    df1 <- df
    t1_power <- pt(qt(size/2,df=df1) + abs(delta)/sqrt(mean(vardelta)), df=df1)
    
    # t-test power with df = # clusters minus p
    # df2 <- I - (J+1)
    # t2_power <- pt(qt(size/2,df=df2) + abs(delta)/sqrt(vardelta), df=df2)
    
    # final results
    # return(data.frame(z_test=z_power,t_test_df_cp=t1_power,t_test_df_cl=t2_power))
  }
  
  return(data.frame(z_test=z_power, t_test_df_cp=t1_power))
}


# ---------------------------------------------------------------------------------

# examples
I = 8; J = 9; Jn = 6; K = 100
CV=0

tnew = J+Jn
design <- matrix(0,J-1,J)
design[upper.tri(design)] <- 1
design <- cbind(design, matrix(1,J-1,tnew-J))
bas = 0.25
end = 0.27   # time trend at the last period
p=seq(bas,end,length.out=J+Jn)
beta = log(p/(1-p))
delta =  log(0.35/(1-0.35)) - beta[tnew]
alpha = c(0.04, 0.02)
swdgeepower_var(I=I, J=J+Jn, K=K, CV=CV, design=design, family="binomial", link="logit", corstr="nested exchangeable", windep=FALSE, delta=delta, beta=beta, phi=1, size=0.025, alpha=alpha)
swdgeepower_var(I=I, J=J+Jn, K=K, CV=CV, design=design, family="binomial", link="logit", corstr="nested exchangeable", windep=TRUE, delta=delta, beta=beta, phi=1, size=0.025, alpha=alpha)

swdgeepower_var(I=I, J=J+Jn, K=K, CV=1, design=design, family="binomial", link="logit", corstr="nested exchangeable", windep=FALSE, delta=delta, beta=beta, phi=1, size=0.025, alpha=alpha)
swdgeepower_var(I=I, J=J+Jn, K=K, CV=1, design=design, family="binomial", link="logit", corstr="nested exchangeable", windep=TRUE, delta=delta, beta=beta, phi=1, size=0.025, alpha=alpha)


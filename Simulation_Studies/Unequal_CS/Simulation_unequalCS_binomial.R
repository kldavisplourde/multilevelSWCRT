# Generating correlated binary outcomes under GLMM for multilevel stepped wedge trials with unequal cluster-sizes
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
# exp.delta: Effect size (in odds ratio)
# bs: baseline time effect
# beta: Vector of period effects (in log odds ratio)
# CV.l: coefficient of variation at the subcluster level
# CV.m: coefficient of variation at the subject level

#############################################################
args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

library(lme4)
library(foreach)
library(doMC)
library(doRNG)

ncores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",1))
if (is.na(ncores)) ncores<-1
registerDoMC(cores=ncores)

# define scenarios
scenarios <- read.table("parameters_binomial_power.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

######   SCENARIO X   ######
scenario	<- k
n <- scenarios$n                         #Number of clusters
l <- scenarios$l                         #Number of clinics
m <- scenarios$m                         #Number of observations within a clinic
t <- scenarios$t                         #Number of time-points
delta <- log(scenarios$exp.delta)
bs <-log(scenarios$bs/(1-scenarios$bs)) #scenarios$bs
beta <- cumsum(c(bs,-0.1,-0.1/2,-0.1/(2^2),-0.1/(2^3),-0.1/(2^4),-0.1/(2^5)))[1:t]
alpha<-c(scenarios$alpha0,scenarios$rho0,scenarios$rho1,scenarios$alpha1)
niteration 	<-scenarios$niteration
seed      <- scenarios$seed + k
CV.l <- scenarios$CV.l
CV.m <- scenarios$CV.m
#############################
paste("THIS IS THE LOG FOR SCENARIO",scenario)
paste("n=",n,"l=",l,"m=",m,"t=",t,"delta=",delta,"bs=",bs,"niteration=",niteration,"seed=",seed,
      "CV.l=",CV.l, "CV.m=",CV.m)
paste("alpha="); paste(alpha)
paste("beta="); paste(beta)
#############################

set.seed(seed)
# First we will focus on closed cohort at clinic level and cross-sectional at individual level (no sig2_g)
scheme<-rep(n/(t-1),t-1)
trtSeq<-matrix(0,t-1,t)
trtSeq[upper.tri(trtSeq)]<-1
trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]

## 1. Generate variance components using given correlations
sig2_e<-pi^2/3
sig2_b<-sig2_e*alpha[3]/(1-alpha[1])                                  #within cluster
sig2_c<-sig2_e*(alpha[4]-alpha[3])/(1-alpha[1])                       #within clinic
sig2_s<-sig2_e*(alpha[2]-alpha[3])/(1-alpha[1])                       #within cluster within period
sig2_p<-sig2_e*(alpha[1]-alpha[4]-alpha[2]+alpha[3])/(1-alpha[1])     #within clinic within period
sig2<-c(sig2_b,sig2_c,sig2_s,sig2_p)

# Functions for generating cluster-period sizes based on CV (we will assume cluster-period size does not change over time)
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

simData <- foreach(i=1:niteration, .combine=rbind) %dopar% {
  ## 2. Generate simulated dataset
  ### a. Generate cluster-period sizes
  l_var <- l.variable(n, t, CV.l, l)
  m_var <- m.variable(n, t, CV.m, m)

  ### b. Generate eta, probability, y, and subject-specific treatment
  for(i in 1:n){
    ### Generate id variables
    dt.i <- data.frame(id=(1:(t*l_var[i]*m_var[i])))
    dt.i$clust <- i
    dt.i$period <- rep(1:t,each=l_var[i]*m_var[i])
    dt.i$clinic <- rep(1:l_var[i],each=m_var[i])
    dt.i$id <- rep(1:m_var[i],l_var[i]*t) 
    
    if(i==1){dt<-dt.i} 
    else{dt<-rbind(dt,dt.i)}
    
    trt <- trtSeq[i,]
    b_i <- rnorm(1,mean=0,sd=sqrt(sig2_b))
    c_il <- rnorm(l_var[i],mean=0,sd=sqrt(sig2_c))
    s_ij <- rnorm(t,mean=0,sd=sqrt(sig2_s))
    p_ijl <- matrix(rnorm(t*l_var[i],mean=0,sd=sqrt(sig2_p)),nrow=l_var[i],ncol=t)
    
    for(j in 1:t){
      for(h in 1:l_var[i]){
        for(k in 1:m_var[i]){
        
          treatment.ijlk<-trt[j]
        
          eta.ijlk <- beta[j] + delta*treatment.ijlk + b_i + c_il[h] + s_ij[j] + p_ijl[h,j]
          p.ijlk<-exp(eta.ijlk)/(1+exp(eta.ijlk))
          y.ijlk<-rbinom(1,1,p.ijlk)
          
          if(i==1 & j==1 & h==1 & k==1){
            y<-y.ijlk
            treatment<-treatment.ijlk} else{
              y <- rbind(y,y.ijlk)
              treatment<-rbind(treatment,treatment.ijlk)}
          }
      }
    }
  }
  dt$trt<-treatment
  dt$y<-y


  ## 2. Fit the model
  dt2<-dt
  dt2$clinic<-as.numeric(paste(dt2$clust, dt2$clinic, sep = ""))
  dt2$cluster.period<-as.numeric(paste(dt2$clust, dt2$period, sep = ""))
  dt2$clinic.period<-as.numeric(paste(dt2$clust, dt2$period, dt2$clinic, sep = ""))

  m1<-glmer(y~factor(period)+treatment+(1|clust)+(1|clinic)+(1|cluster.period)+(1|clinic.period),data=dt2,family=binomial)

  fixed<-fixef(m1)
  se<-sqrt(diag(vcov(m1)))

  c(fixed,se)
}

if(t==3){
  colnames(simData)<-c("Intercept.est","Period2.est","Period3.est","Treatment.est",
                       "Intercept.se","Period2.se","Period3.se","Treatment.se")
}

if(t==4){
  colnames(simData)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Treatment.est",
                       "Intercept.se","Period2.se","Period3.se","Period4.se","Treatment.se")
}

if(t==5){
  colnames(simData)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Period5.est","Treatment.est",
                       "Intercept.se","Period2.se","Period3.se","Period4.se","Period5.se","Treatment.se")
}

if(t==6){
  colnames(simData)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Period5.est","Period6.est","Treatment.est",
                       "Intercept.se","Period2.se","Period3.se","Period4.se","Period5.se","Period6.se","Treatment.se")
}

if(t==7){
  colnames(simData)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Period5.est","Period6.est","Period7.est","Treatment.est",
                       "Intercept.se","Period2.se","Period3.se","Period4.se","Period5.se","Period6.se","Period7.se","Treatment.se")
}

if(delta==0) analysis<-"error"
if(delta != 0) analysis<-"power"

simData <- as.data.frame(simData)
write.table(simData, file=paste("results/BinResults_uneqCS_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)







# Approximation of the variance of the intervention effect for stepped wedge cluster randomized trials with subclusters
# and Gaussian or binary (with logit link) outcomes.
# Created by Kendra Davis-Plourde

# INPUTS
# n: Number of clusters (I)
# l: Number of subclusters per cluster (K)
# m: Number of participants per subcluster (N)
# t: Number of periods (T)
# subcluster: cohort or cross-sectional on subcluster level
# indiv: cohort or cross-sectional on subject level
# family: gaussian or binomial (with canonical logit link)
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
#         alpha_2=within-individual correlation
# delta: Effect of intervention in log odds ratio (only needed when outcome is binomial)
# beta: Vector of period effects (only needed when outcome is binomial)
# phi: common dispersion (only needed when outcome is binomial)
# tot.var: total variance of the outcome (only needed when outcome is Gaussian)

VARd<-function(n,l,m,t,subcluster=c("cohort","cross-sectional"),indiv=c("cohort","cross-sectional"),family=c("gaussian","binomial"),
                 alpha,delta=NA,beta=NA,phi=1,tot.var=1){
  
  # elements of efficiency calculation
  scheme<-rep(n/(t-1),t-1)
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]
  if(sum(scheme) != n){stop("invalid randomization scheme")}
  
  # loop through each cluster
  if (family=="binomial") {
    I.T<-diag(t)
    J.T<-matrix(1,nrow=t,ncol=t)
    Omega<-matrix(0,t+1,t+1)
    sig2_e<-pi^2/3
    # Generate variance components using given correlations
    if (subcluster=="cohort" & indiv=="cohort"){
      sig2_b<-sig2_e*alpha[3]/(1-alpha[1]-alpha[5]+alpha[4])                                #within cluster
      sig2_c<-sig2_e*(alpha[4]-alpha[3])/(1-alpha[1]-alpha[5]+alpha[4])                     #within subcluster
      sig2_s<-sig2_e*(alpha[2]-alpha[3])/(1-alpha[1]-alpha[5]+alpha[4])                     #within cluster within period
      sig2_p<-sig2_e*(alpha[1]-alpha[4]-alpha[2]+alpha[3])/(1-alpha[1]-alpha[5]+alpha[4])   #within subcluster within period
      sig2_g<-sig2_e*(alpha[5]-alpha[4])/(1-alpha[1]-alpha[5]+alpha[4])                     #within person
    }
    
    if (subcluster=="cohort" & indiv=="cross-sectional"){
      sig2_b<-sig2_e*alpha[3]/(1-alpha[1])                                  #within cluster
      sig2_c<-sig2_e*(alpha[4]-alpha[3])/(1-alpha[1])                       #within subcluster
      sig2_s<-sig2_e*(alpha[2]-alpha[3])/(1-alpha[1])                       #within cluster within period
      sig2_p<-sig2_e*(alpha[1]-alpha[4]-alpha[2]+alpha[3])/(1-alpha[1])     #within subcluster within period
      sig2_g<-0                                                             #within person
    }
    
    if (subcluster=="cross-sectional" & indiv=="cross-sectional") {
      sig2_b<-sig2_e*alpha[3]/(1-alpha[1])                                  #within cluster
      sig2_c<-0                                                             #within subcluster
      sig2_s<-sig2_e*(alpha[2]-alpha[3])/(1-alpha[1])                       #within cluster within period
      sig2_p<-sig2_e*(alpha[1]-alpha[2])/(1-alpha[1])                       #within subcluster within period
      sig2_g<-0                                                             #within person
    }
    
    for(i in 1:n){
      
      # design matrix
      X<-trtSeq[i,]
      Z<-cbind(X,I.T)
      
      # unique marginal means
      gmu<-c(Z%*%c(delta,beta))
      E.Gprime<-as.vector(2+exp(0.5*(sig2_b+sig2_c+sig2_s+sig2_p+sig2_g))*(exp(-gmu)+exp(gmu)))
      
      V<-(E.Gprime/(l*m)+sig2_p/l+sig2_s)*I.T+(sig2_b+sig2_c/l+sig2_g/(l*m))*J.T
      invV<-solve(V)
      
      Omega<-Omega + t(Z)%*%invV%*%Z
    }
    vardelta<-phi*solve(Omega)[1,1]
  }
 
  if (family=="gaussian"){
    U <- sum(trtSeq)
    V <- sum(rowSums(trtSeq)^2)
    W <- sum(colSums(trtSeq)^2)
    
    #Derived eigenvalues
    if (subcluster=="cohort" & indiv=="cohort"){
      eigen3 <- 1-alpha[1]-alpha[5]+alpha[4] + m*(alpha[1]-alpha[4]+(l-1)*(alpha[2]-alpha[3]))
      eigen6 <- 1-alpha[1]+(t-1)*(alpha[5]-alpha[4]) + m*(alpha[1]+(t-1)*alpha[4]+(l-1)*(alpha[2]+(t-1)*alpha[3]))
    }
    if (subcluster=="cohort" & indiv=="cross-sectional"){
      eigen3 <- 1-alpha[1]+ m*(alpha[1]-alpha[4]+(l-1)*(alpha[2]-alpha[3]))
      eigen6 <- 1-alpha[1] + m*(alpha[1]+(t-1)*alpha[4]+(l-1)*(alpha[2]+(t-1)*alpha[3]))
    }
    if (subcluster=="cross-sectional" & indiv=="cross-sectional"){
      eigen3 <- 1-alpha[1]+ m*(alpha[1]-alpha[3]+(l-1)*(alpha[2]-alpha[3]))
      eigen6 <- 1-alpha[1] + m*(alpha[1]+(t-1)*alpha[3]+(l-1)*(alpha[2]+(t-1)*alpha[3]))
    }
    vardelta<-((tot.var/(l*m))*n*t*eigen3*eigen6)/((U^2 + n*t*U - t*W - n*V)*eigen6 - (U^2 - n*V)*eigen3)
  }
  
  return(vardelta)
}
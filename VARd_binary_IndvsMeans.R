# Comparing the variance approximation of the intervention effect estimator (VARd) using individual level data versus 
# using cluster-period means approach under a binary outcome with canonical logit link and assuming common dispersion=1.

library(optimbase)

# VARd using individual level data
indVAR<-function(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g){
  # elements of efficiency calculation
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]
  if(sum(scheme) != I){stop("invalid randomization scheme")}
  I.T<-diag(t)
  J.T<-ones(t)
  I.LN<-diag(L*N)
  J.LN<-ones(L*N)
  IL.kron.JN<-kronecker(diag(L),ones(N))
  Omega<-matrix(0,t+1,t+1)
  
  # loop through each cluster
  for(i in 1:I){
    
    # design matrix
    X<-trtSeq[i,]
    Z<-cbind(X,I.T)
    
    # unique marginal means
    gmu<-c(Z%*%c(delta,beta))
    E.Gprime<-as.vector(2+exp(0.5*(sig2.b+sig2.c+sig2.s+sig2.p+sig2.g))*(exp(-gmu)+exp(gmu)))
    E.Gprime.mat<-kronecker(E.Gprime*I.T,diag(L*N))
      
    bm1<-kronecker(I.T,sig2.s*J.LN+sig2.p*IL.kron.JN)
    bm2<-kronecker(J.T,sig2.b*J.LN+sig2.c*IL.kron.JN+sig2.g*I.LN)
    kronecker(J.T,sig2.b*J.LN)+kronecker(J.T,sig2.c*IL.kron.JN)+kronecker(J.T,sig2.g*I.LN)
    V<-E.Gprime.mat+bm1+bm2
    invV<-solve(V)
    D<-kronecker(Z,transpose(matrix(rep(1,n=(L*N)),1,L*N)))
    Omega<-Omega + t(D)%*%invV%*%D
  }
  vardelta<-solve(Omega)[1,1]
  return(vardelta)
}

# VARd using a cluster-period means approach
sumVAR<-function(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g){
  # elements of efficiency calculation
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), scheme),]
  if(sum(scheme) != I){stop("invalid randomization scheme")}
  I.T<-diag(t)
  J.T<-ones(t)
  Omega<-matrix(0,t+1,t+1)
  
  # loop through each cluster
  for(i in 1:I){
    
    # design matrix
    X<-trtSeq[i,]
    Z<-cbind(X,I.T)
    
    # unique marginal means
    gmu<-c(Z%*%c(delta,beta))
    E.Gprime<-as.vector(2+exp(0.5*(sig2.b+sig2.c+sig2.s+sig2.p+sig2.g))*(exp(-gmu)+exp(gmu)))
    
    V<-(E.Gprime/(L*N)+sig2.p/L+sig2.s)*I.T+(sig2.b+sig2.c/L+sig2.g/(L*N))*J.T
    invV<-solve(V)
    
    Omega<-Omega + t(Z)%*%invV%*%Z
  }
  vardelta<-solve(Omega)[1,1]
  return(vardelta)
}

#################################################################################
#################################################################################


################ Quick Accuracy test ################
# Inputs - assuming equal number of subclusters and participants across cluster-periods
I<-4                             #Number of clusters (I)
t<-5                             #Number of periods (T)
L<-2                             #Number of subclusters per cluster (L)
N<-20                            #Number of participants per subcluster (N)
scheme<-c(1,1,1,1)               #scheme: randomization scheme
delta<-2                         #delta - effect size (log(OR))
beta<-c(0.2,0.4,0.6,0.8,1)       #beta - period effect (log(OR))
alpha<-c(0.5,0.3,0.25,0.15)      #ICCs in the following order:
#                                   alpha_0=within subcluster within period
#                                   rho_0=between subcluster within period
#                                   alpha_1=within subcluster between period
#                                   rho_1=between subcluster between period
#                                   alpha_2=within-individual correlation
sig2.b<-0.2                      #Variance of cluster random-effect
sig2.c<-0.3                      #Variance of subcluster random-effect
sig2.s<-0.25                     #Variance of cluster-time random effect
sig2.p<-0.35                     #Variance of subcluster-time random effect
sig2.g<-0.5                      #Variance of participant random effect

# Comparing strategies
indVAR(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g)
sumVAR(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g)

################ Speed test ##########################
# Inputs - assuming equal number of subclusters and participants across cluster-periods
I<-4                             #Number of clusters (I)
t<-5                             #Number of periods (T)
L<-6                             #Number of subclusters per cluster (L)
N<-50                            #Number of participants per subcluster (N)
scheme<-rep(1,4)                 #scheme: randomization scheme
delta<-2                         #delta - effect size (log(OR))
beta<-c(0.2,0.4,0.6,0.8,1)       #beta - period effect (log(OR))
sig2.b<-0.2                      #Variance of cluster random-effect
sig2.c<-0.3                      #Variance of subcluster random-effect
sig2.s<-0.25                     #Variance of cluster-time random effect
sig2.p<-0.35                     #Variance of subcluster-time random effect
sig2.g<-0.5                      #Variance of participant random effect

# Comparing strategies
indVAR(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g) #Approximately 10 seconds
sumVAR(I,t,L,N,scheme,delta,beta,alpha,sig2.b,sig2.c,sig2.s,sig2.p,sig2.g) #Less than a second

#################################################################################
#################################################################################
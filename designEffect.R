# Design effect for longitudinal cluster randomized trials with subclusters and a Gaussian outcome.
# Created by Kendra Davis-Plourde

# INPUTS
# n: Number of clusters (I)
# l: Number of subclusters per cluster (K)
# m: Number of participants per subcluster (N)
# t: Number of periods (T)
# U: Design constant
# V: Design constant
# W: Design constant
# subcluster: cohort or cross-sectional
# indiv: cohort or cross-sectional
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
#         alpha_2=within-individual correlation

designEffect<-function(n,l,m,t,U,V,W,subcluster=c("cohort","cross-sectional"),indiv=c("cohort","cross-sectional"),alpha){
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
  design_effect<-((n^2)*t*eigen3*eigen6)/(4*((U^2)+(n*t*U)-(t*W)-(n*V))*eigen6-4*((U^2)-(n*V))*eigen3)
  return(design_effect)
}
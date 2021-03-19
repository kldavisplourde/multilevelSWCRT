# Design effect for stepped wedge cluster randomized trials with a continuous outcome

# INPUTS
# l: Number of subclusters per cluster (L)
# m: Number of participants per subcluster (N)
# t: Number of periods (T)
# S: Number of steps
# c: Number of measures taken at each step
# b: Number of measures taken at baseline
# subcluster: cohort or cross-sectional
# indiv: cohort or cross-sectional
# alpha: Vector of ICCs in the following order:
#         alpha_0=within subcluster within period
#         rho_0=between subcluster within period
#         rho_1=between subcluster between period
#         alpha_1=within subcluster between period
#         alpha_2=within-individual correlation

designEffect_SWCRT<-function(l,m,t,S,c,b,subcluster=c("cohort","cross-sectional"),indiv=c("cohort","cross-sectional"),alpha){
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
  design_effect<-3/(2*c*(S-(1/S)))*(((b+S*c)*eigen3*eigen6)/((S*c/2)*eigen3+(b+S*c/2)*eigen6))
  return(design_effect)
}
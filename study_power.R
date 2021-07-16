# Generating power of clinical trial
# Created by Kendra Davis-Plourde

# INPUTS
# delta: Intervention effect on link function scale
# var.delta: Variance of intervention effect on link function scale
# typeI.error: Type I error rate for t-test (default is 0.05)
# df: degrees of freedom for t-test

study_power<-function(delta,var.delta,typeI.error=0.05,df){
  tE1<-qt(1-typeI.error/2, df=df)
  tE2<-abs(delta)/sqrt(var.delta) - tE1
  power<-pt(tE2, df=df)
  return(power)
}
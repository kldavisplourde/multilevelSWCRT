# Generating power of clinical trial
# Created by Kendra Davis-Plourde

# INPUTS
# delta: Intervention effect on link function scale
# var.delta: Variance of intervention effect on link function scale
# alpha: Type I error rate for t-test (default is 0.05)
# df: degrees of freedom for t-test

study_power<-function(delta, var.delta, alpha=0.05, df=n-2){
  tcrit<-qt(1-alpha/2, df=df)
  power<-1-pt(tcrit, df, abs(delta)/sqrt(var.delta))
  return(power)
}